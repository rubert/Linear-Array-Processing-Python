class rfClass(object):

	roiX = None
	roiY = None
	bModeImageCount = 0
		
	def __init__(self, filename, dataType):
		'''Input:  
		   filename:  The name of the file or directory the data is from 
		   dataType:  A lowercase string indicating the type of the data, choices are: 
		   		ei
				rfd
				rf
				sim
				multiSim
				 '''
		
		if not( dataType == 'ei' or dataType == 'rfd' or dataType == 'rf' or dataType =='sim' or dataType == 'multiSim' ):
			print 'Error.  Datatype must be ei,rfd,rf '
			return
		
		self.fname = filename
		self.dataType = dataType

		self.soundSpeed = 1540. #m/s
		self.data = None
												
		if dataType == 'ei':
			import os, numpy
			startdir = os.getcwd()
			os.chdir(self.fname)
			
			metaData = open('metaData.txt', 'r')
			metaList = [0,0,0,0]
			for i in range(4):
				metaList[i] = int(metaData.readline() )
		
			metaData.close()
			nFrames = len( [name for name in os.listdir('.') if os.path.isfile(name)])
			nFrames -= 1 #because of metaData.txt
			points = len( range( metaList[3]/2, metaList[2] + metaList[3]/2 ) )
			aLines = metaList[1]
			
			self.nFrames = nFrames
			self.soundSpeed = 1540.  #m/s
			self.fs = 40.*10**6  #Hz
			self.deltaX = 40./aLines  #assume a 40 mm FOV
			self.deltaY = self.soundSpeed/(2*self.fs)*10**3
			self.fovX = self.deltaX*(aLines-1)
			self.fovY = self.deltaY*(points -1 )	
			os.chdir(startdir)

			self.points = points
			self.lines = aLines	
	
		if dataType == 'rfd':
			import readrfd
			hInfo = readrfd.headerInfo(filename)
			self.fs = hInfo[1]
			self.deltaX = (1/hInfo[2])*10.
			#read a frame to find out number of points and A lines
			tempFrame = readrfd.readrfd(filename, 1)
			self.nFrames = readrfd.rfdframes(filename)
			self.deltaY = self.soundSpeed/(2*self.fs)*10**3
			self.fovX = self.deltaX*(tempFrame.shape[1] - 1)
			self.fovY = self.deltaY*(tempFrame.shape[0] -1 )	
			self.points = tempFrame.shape[0]
			self.lines = tempFrame.shape[1]	
	

		if dataType == 'rf':
			#The header is 19 32-bit integers
			import numpy 
			dataFile = open(self.fname, 'r')
			header = numpy.fromfile(dataFile, numpy.int32, 19)
			self.header = header
			self.fs = header[15]	
		        self.deltaY = self.soundSpeed/(2*self.fs)*10**3	
			self.points = header[3]
			self.lines = header[2]
			self.deltaX = 40./self.lines #parameter in header is broken, so assume a 40 mm FOV
			self.fovX = self.lines*self.deltaX
			self.fovY = self.points*self.deltaY		
			self.nFrames = header[1]
		
		if dataType == 'sim':
				
			f = open(filename, 'rb')
			
			import numpy as np
			#in Hz
			self.freqstep =float( np.fromfile(f, np.double,1) )
			self.points = int( np.fromfile(f, np.int32,1) )
			self.lines = int( np.fromfile(f, np.int32,1) )
			tempReal = np.fromfile(f,np.double, self.points*self.lines )
			tempImag = np.fromfile(f, np.double, self.points*self.lines )
			f.close()
			
			self.deltaX = .2  #beamspacing in mm
			self.fovX = self.lines*self.deltaX

		if dataType == 'multiSim':
				
			f = open(filename, 'rb')
			
			import numpy as np
			#in Hz
			self.freqstep =float( np.fromfile(f, np.double,1) )
			self.points = int( np.fromfile(f, np.int32,1) )
			self.lines = int( np.fromfile(f, np.int32,1) )
			self.nFrames = int(np.fromfile(f, np.int32,1) )
			tempReal = np.fromfile(f,np.double, self.points*self.lines )
			tempImag = np.fromfile(f, np.double, self.points*self.lines )
			f.close()
			
			self.deltaX = .2  #beamspacing in mm
			self.fovX = self.lines*self.deltaX
					
			
	def ReadFrame(self, frameNo = 0, centerFreq = 5.0E6, sigma = 1.0E6, samplingFrequency = 40.0E6):
		'''Read a single frame from the input file, the method varies depending on the file type that
		was read in.  This can handle linear array data from the Seimens S2000 recorded in either mode
		or the Ultrasonix scanner recorded using the clinical software.
		For a simulated data file, the 2nd and 3rd parameters matter, not for real data.
		
		Input::
		frameNo.  The RF data frame from a file obtained off a real ultrasound scanner.
		CenterFreq: (Hz)  The center frequency of the transmitted pulse for simulated ultrasound data.
		Sigma. (Hz)  The standard deviation of the Gaussian shaped transmit pulse for simulated data.
		samplingFrequency. (Hz)  The simulated sampling frequency for a simulated data set.  Higher
		sampling frequencies are achieved through zero-padding.'''

		if self.dataType == 'ei':
			import os, numpy
			startdir = os.getcwd()
			os.chdir(self.fname)
			
			metaData = open('metaData.txt', 'r')
			metaList = [0,0,0,0]
			for i in range(4):
				metaList[i] = int(metaData.readline() )
		
			metaData.close()
			points = len( range( metaList[3]/2, metaList[2] + metaList[3]/2 ) )
			aLines = metaList[1]
			from numpy import zeros, fromfile
			self.data = zeros( (points, aLines) )
			readThisMany = (metaList[2]+ metaList[3]/2)*metaList[1]
			
			tempIn = open(str(frameNo), 'r')
			frame = fromfile(tempIn,  numpy.int16, readThisMany)
			frame = frame.reshape( ( metaList[2]+ metaList[3]/2, metaList[1] )  , order = 'F' )
			self.data[:,:] = frame[metaList[3]/2:metaList[3]/2 + metaList[2], :].copy('C')

			os.chdir(startdir)


		if self.dataType == 'rfd':
			import readrfd
			self.data = readrfd.readrfd(self.fname, frameNo + 1)

		if self.dataType == 'rf':
			#Read in the header, then read in up to one less frame than necessary
			import numpy
			dataFile = open(self.fname, 'r')
			header = numpy.fromfile(dataFile, numpy.int32, 19)
			if frameNo == 0:
				temp = numpy.fromfile( dataFile, numpy.int16, self.points*self.lines)
				self.data = temp.reshape( (self.points, self.lines),order='F')
			else:
				dataFile.seek( 2*self.points*self.lines*(frameNo-1), 1) 
				temp = numpy.fromfile( dataFile, numpy.int16, self.points*self.lines)
				self.data = temp.reshape( (self.points, self.lines), order='F')


		if self.dataType == 'sim':
			####To get the actual sampling frequency we will zero-pad up to the
			###desired sampling frequency
					
			self.centerFreq = centerFreq #Hz
			self.sigma = sigma
			f = open(self.fname, 'rb')
			
			import numpy as np
			#in Hz
			self.freqstep =float( np.fromfile(f, np.double,1) )
			self.points = int( np.fromfile(f, np.int32,1) )
			self.lines = int( np.fromfile(f, np.int32,1) )
			tempReal = np.fromfile(f,np.double, self.points*self.lines )
			tempImag = np.fromfile(f, np.double, self.points*self.lines )
			f.close()
		
		
			pointsToDesiredFs = int(samplingFrequency/self.freqstep)
			self.fs = pointsToDesiredFs*self.freqstep
			self.freqData = (tempReal - 1j*tempImag).reshape( (self.points, self.lines), order = 'F' )
					
			import math
			f = np.arange(0,(self.points)*self.freqstep, self.freqstep)
			self.pulseSpectrum = np.exp(-(f - self.centerFreq)*(f - self.centerFreq)/(2*self.sigma**2) )
			temp = self.freqData*self.pulseSpectrum.reshape( (self.points, 1) )	

			zeroPadded = np.zeros( (pointsToDesiredFs, self.lines) ) + 1j*np.zeros( (pointsToDesiredFs, self.lines) )
			zeroPadded[0:self.points, :] = temp	
			self.data = np.fft.ifft(zeroPadded,axis=0 ).real	
			self.points = pointsToDesiredFs
			
			self.fs = self.freqstep*self.points
			self.deltaY = 1540. / (2*self.fs)*10**3
			self.fovY = self.deltaY*self.points


		if self.dataType == 'multiSim':
			####To get the actual sampling frequency we will zero-pad up to the
			###desired sampling frequency
					
			self.centerFreq = centerFreq #Hz
			self.sigma = sigma
			f = open(self.fname, 'rb')
			
			import numpy as np
			#in Hz
			self.freqstep =float( np.fromfile(f, np.double,1) )
			self.points = int( np.fromfile(f, np.int32,1) )
			self.lines = int( np.fromfile(f, np.int32,1) )
			self.nFrames = int(np.fromfile(f, np.int32, 1) )

			#now jump ahead by as many frames as necessary
			#assume 64 bit, 8 byte doubles
			f.seek(8*2*self.points*self.lines*frameNo, 1)
			#read in desired frame
			tempReal = np.fromfile(f,np.double, self.points*self.lines )
			tempImag = np.fromfile(f, np.double, self.points*self.lines )
			f.close()
		
			pointsToDesiredFs = int(samplingFrequency/self.freqstep)
			self.fs = pointsToDesiredFs*self.freqstep
			self.freqData = (tempReal - 1j*tempImag).reshape( (self.points, self.lines), order = 'F' )
					
			import math
			f = np.arange(0,(self.points)*self.freqstep, self.freqstep)
			self.pulseSpectrum = np.exp(-(f - self.centerFreq)*(f - self.centerFreq)/(2*self.sigma**2) )
			temp = self.freqData*self.pulseSpectrum.reshape( (self.points, 1) )	

			zeroPadded = np.zeros( (pointsToDesiredFs, self.lines) ) + 1j*np.zeros( (pointsToDesiredFs, self.lines) )
			zeroPadded[0:self.points, :] = temp	
			self.data = np.fft.ifft(zeroPadded,axis=0 ).real	
			self.points = pointsToDesiredFs
			
			self.fs = self.freqstep*self.points
			self.deltaY = 1540. / (2*self.fs)*10**3
			self.fovY = self.deltaY*self.points
				


	def MakeBmodeImage(self, frameNo = 0, showFig = True):
		"""Create a B-mode image and immediately display it.  This is useful for interactive viewing of particular
		frames from a data set."""
		
		rfClass.bModeImageCount +=1	
		self.ReadFrame(frameNo)
		temp = self.data


		#import signal processing modules and generate Numpy array
		from scipy.signal import hilbert
		from numpy import log10
		bMode = log10(abs(hilbert(temp, axis = 0)))
		bMode = bMode - bMode.max()

		#import matplotlib and create plot
		import matplotlib.pyplot as plt
		import matplotlib.cm as cm
		fig = plt.figure()
		fig.canvas.set_window_title("B-mode image " + str(rfClass.bModeImageCount) )
		ax = fig.add_subplot(1,1,1)
		ax.imshow(bMode, cmap = cm.gray, vmin = -3, vmax = 0, extent = [0, self.fovX,  self.fovY, 0])
		
		if showFig:
			plt.show()			

		

	def SetRoiFixedSize(self, windowX = 4, windowY = 4):
		'''The Roi size here is given in mm'''

		
		from matplotlib import pyplot
		from scipy.signal import hilbert
		from numpy import log10, arange
		
		import types
		if type(self.data) == types.NoneType:	
			self.ReadFrame()
		yExtent = self.deltaY*self.points
		xExtent = self.deltaX*self.lines
		bmode = log10(abs(hilbert(self.data, axis=0) ) )
		bmode = bmode - bmode.max()
		pyplot.imshow(bmode, extent = [0,xExtent,yExtent,0], cmap = 'gray', vmin = -3, vmax = 0 )
		x= pyplot.ginput()
		pyplot.close()
		#Compute the size of the analysis window in points and A-lines
		self.windowY = int(windowY/self.deltaY)
		self.windowX = int(windowX/self.deltaX)

		#Compute the center of the bounding box in pixels
		boxX = int( x[0][0]/self.deltaX )
		boxY = int( x[0][1]/self.deltaY )

		xLow = boxX - self.windowX/2
		if xLow < 0:
			xLow = 0
		xHigh = xLow + self.windowX
		if xHigh > self.lines:
			xHigh = self.lines
		xLow = xHigh - self.windowX

		yLow = boxY - self.windowY/2
		if yLow < 0:
			yLow = 0
		yHigh = yLow + self.windowY
		if yHigh > self.points+1:
			yHigh = self.points+1
		yLow = yHigh - self.windowY

		rfClass.roiX = [xLow, xHigh]
		rfClass.roiY = [yLow, yHigh]

		#Draw Analysis Window
		import matplotlib.cm as cm
		palette = cm.gray
		palette.set_bad('r')
		from numpy import ma, zeros
		bmode = ma.array(bmode)
		bmode.mask = zeros(bmode.shape)


		bmode.mask[yLow:yHigh, xLow] = True
		bmode.mask[yLow:yHigh, xHigh] = True
		#since the pixels are so much thinner in the axial direction we use ten of them to
		#draw the ROI
		lineThickY = 10
		if yLow < lineThickY:
			yLow = lineThickY 
		bmode.mask[yLow-lineThickY:yLow, xLow:xHigh] = True

		if yHigh > bmode.mask.shape[0] - lineThickY:
			yHigh = bmode.mask.shape[0] - lineThickY 
		bmode.mask[yHigh-lineThickY:yHigh, xLow:xHigh] = True


		pyplot.imshow(bmode, cmap = palette, extent = [0, xExtent, yExtent, 0 ],vmin = -3, vmax = 0)
		pyplot.show()


	def SetRoiBoxSelect(self):
		"""
		Do a mouseclick somewhere, move the mouse to some destination, release
		the button.  This class gives click- and release-events and also draws
		a line or a box from the click-point to the actual mouseposition
		(within the same axes) until the button is released.  Within the
		method 'self.ignore()' it is checked wether the button from eventpress
		and eventrelease are the same.

		"""
		#Check that a Data set has been loaded	
		if self.fname == None:
			return	

		from matplotlib.widgets import RectangleSelector
		import numpy as np
		import matplotlib.pyplot as plt
		current_ax = plt.subplot(111)               # make a new plotingrangej


		def on_select(eclick, erelease):
			rfClass.roiX = [0,0]
			rfClass.roiY = [0,0]
			#lower boundary, make sure it doesn't violate block matching constraints
			tempX	= int(erelease.xdata/self.deltaX)
			if tempX < self.minimumROIBoundaryX:
				tempX = self.minimumROIBoundaryX
			if tempX > self.data.shape[1] - self.minimumROIBoundaryX - 1:
				tempX = self.data.shape[1] -self.minimumROIBoundaryX - 1			
			rfClass.roiX[0] = tempX
			
			tempX	= int(eclick.xdata/self.deltaX)
			if tempX < self.minimumROIBoundaryX:
				tempX = self.minimumROIBoundaryX
			if tempX > self.data.shape[1] - self.minimumROIBoundaryX - 1:
				tempX = self.data.shape[1] -self.minimumROIBoundaryX - 1			
			rfClass.roiX[1] = tempX
			
			tempY	= int(eclick.ydata/self.deltaY)
			if tempY < self.minimumROIBoundaryY:
				tempY = self.minimumROIBoundaryY
			if tempY > self.data.shape[0] - self.minimumROIBoundaryY - 1:
				tempY = self.data.shape[0] -self.minimumROIBoundaryY - 1			
			rfClass.roiY[0] = tempY

			tempY	= int(erelease.ydata/self.deltaY)
			if tempY < self.minimumROIBoundaryY:
				tempY = self.minimumROIBoundaryY
			if tempY > self.data.shape[0] - self.minimumROIBoundaryY - 1:
				tempY = self.data.shape[0] -self.minimumROIBoundaryY - 1			
			rfClass.roiY[1] = tempY


		# drawtype is 'box' or 'line' or 'none'
		rs = RectangleSelector(current_ax, on_select,
						       drawtype='box', useblit=True,
						       button=[1,3], # don't use middle button
						       minspanx=5, minspany=5,
						       spancoords='data')
			
		
		#could be image sequence or just a 2-D image
		import types
		if type(self.data) == types.NoneType:	
			self.ReadFrame(0)
		temp = self.data


		from scipy.signal import hilbert
		from numpy import log10
		bMode = log10(abs(hilbert(temp, axis = 0)))
		bMode = bMode - bMode.max()
		bMode[bMode < -3] = -3

		#import matplotlib and create plot
		import matplotlib.cm as cm

		plt.imshow(bMode, cmap = cm.gray,  extent = [0, self.fovX,  self.fovY, 0])
		plt.show()


	def CreateParametricImage(self, paramImage, origin, spacing, colormap = 'jet'):
		'''Input:
		   paramImage: The image to place within the B-mode image
		   origin:  The B-mode pixel at which the upper left hand corner of the parametric image is located.
		   		[row, column]
		   spacing:  The number of B-mode pixels separating paramImage's pixels
		   		[rows, columns]
		   colormap:  The colormap of the parametric image
		   
		 '''	

		from numpy import arange,zeros		
		from scipy import interpolate

		self.paramValMax = paramImage.max()
		self.paramValMin = paramImage.min()	
		#work out size of region in B-mode image
		bModeSizeY = spacing[0]*(paramImage.shape[0] - 1) + 1
		bModeSizeX = spacing[1]*(paramImage.shape[1] - 1) + 1

		paramImageUpInY = zeros( (bModeSizeY, paramImage.shape[1] ) )
		paramImageUp = zeros( (bModeSizeY, bModeSizeX) )
					 
		#Upsample strain image to same spacing as B-mode image.
		paramRfYIndexes = arange(origin[0], origin[0] + spacing[0]*paramImage.shape[0], spacing[0] )
		paramRfYIndexesNew = arange(origin[0], origin[0] + spacing[0]*(paramImage.shape[0]-1) + 1 )
		paramRfXIndexes = arange(origin[1], origin[1] + spacing[1]*paramImage.shape[1], spacing[1] )
		paramRfXIndexesNew = arange(origin[1], origin[1] + spacing[1]*(paramImage.shape[1]-1) + 1 )

		#use old image to interpolate
		for x in range(paramImage.shape[1]):	
			interp = interpolate.interp1d(paramRfYIndexes, paramImage[:,x] )
			paramImageUpInY[:, x] = interp(paramRfYIndexesNew )

		
		for y in range(paramImageUp.shape[0]):	
			interp = interpolate.interp1d(paramRfXIndexes, paramImageUpInY[y,:] )
			paramImageUp[y, :] = interp( paramRfXIndexesNew )



		'''Convert array containing param values to RGBALpha array'''
		from matplotlib import cm

		palette = cm.ScalarMappable()
		palette.set_cmap(colormap)
		tempIm = palette.to_rgba(paramImageUp)

		'''Create B-mode image'''	
		from scipy.signal import hilbert
		from numpy import log10
		bMode = log10(abs(hilbert(self.data, axis = 0)))
		bMode = bMode - bMode.max()
		bMode[bMode < -3]  = -3
		palette = cm.ScalarMappable()
		palette.set_cmap('gray')
		bMode = palette.to_rgba(bMode)

		bMode[origin[0]:origin[0] + tempIm.shape[0],origin[1]:origin[1] + tempIm.shape[1], :] = tempIm
	
	
		return bMode
		
				
	def CreateRoiArray(self):
#		#Check that a Data set has been loaded	
		#could be image sequence or just a 2-D image
		import types
		if type(self.data) == types.NoneType:	
			self.ReadFrame(0)
		temp = self.data


		#Create masked array with B-mode data, mask entries for box
		from scipy.signal import hilbert
		from numpy import log10
		bMode = log10(abs(hilbert(temp, axis = 0)))
		bMode = bMode - bMode.max()
		bMode[bMode < -3] = -3

		from numpy import ma, zeros
		bMode = ma.array(bMode)
		bMode.mask = zeros(bMode.shape)
		#replace some B-mode data with sentinels
		xLow = min(rfClass.roiX); xHigh = max(rfClass.roiX)
		yLow = min(rfClass.roiY); yHigh = max(rfClass.roiY)
		
		bMode.mask[yLow:yHigh, xLow] = True
		bMode.mask[yLow:yHigh, xHigh] = True
		#since the pixels are so much thinner in the axial direction we use ten of them to
		#draw the ROI
		lineThickY = 10
		if yLow < lineThickY:
			yLow = lineThickY 
		bMode.mask[yLow-lineThickY:yLow, xLow:xHigh] = True

		if yHigh > bMode.mask.shape[0] - lineThickY:
			yHigh = bMode.mask.shape[0] - lineThickY 
		bMode.mask[yHigh-lineThickY:yHigh, xLow:xHigh] = True

		return bMode

	def ShowRoiImage(self):
		bMode = self.createRoiArray()	
		#import matplotlib and create plot
		import matplotlib.pyplot as plt
		import matplotlib.cm as cm
		palette = cm.gray
		palette.set_bad('r')

		im = plt.imshow(bMode, cmap = palette, extent = [0, self.fovX, self.fovY, 0])
		plt.show()



	def SaveRoiImage(self, fname):
		bMode = self.createRoiArray()	
		#import matplotlib and create plot
		import matplotlib.pyplot as plt
		import matplotlib.cm as cm
		palette = cm.gray
		palette.set_bad('r')
	
		im = plt.imshow(bMode, cmap = palette, extent = [0, self.fovX, self.fovY, 0])
		plt.savefig(fname)
		plt.close()

	
