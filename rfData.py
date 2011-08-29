class rfClass:

	roiX = None
	roiY = None
		
	def __init__(self, filename, dataType):
		'''Input:  filename:  The name of the file or directory the data is from \n
		   dataType:  A lowercase string indicating the type of the data, choices are: \n
		   		ei
				rfd
				 '''
		
		if not( dataType == 'ei' or dataType == 'rfd' or dataType == 'rf' ):
			print 'Error.  Datatype must be ei,rfd,rf '
			return

		self.soundSpeed = 1540. #m/s
		self.fs = 40.*10**6 #Default sampling frequency in Hz
		self.deltaY = self.soundSpeed/(2*self.fs)*10**3 #Pixel spacing in mm
		self.data = None
		self.isSectorData = False
		self.minimumROIBoundaryY = 0  #for processing strain data
		self.minimumROIBoundaryX = 0

		#set parameters for block matching
		from blockMatch import blockMatchingParameters
		self.params = blockMatchingParameters()
		self.strainLimit = .02
		self.strain = None

								
		self.fname = filename
		self.dataType = dataType

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
		
		
	def readFrame(self, frameNo = 0):
		'''Read a single frame from the input file, the method varies depending on the file type that
		was read in.  This can handle linear array data from the Seimens S2000 recorded in either mode
		or the Ultrasonix scanner recorded using the clinical software.'''
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

		
	def makeBmodeImage(self, frameNo = 0):
		"""Create a B-mode image and immediately display it.  This is useful for interactive viewing of particular
		frames from a data set."""
			
		#Check that a Data set has been loaded	
		self.readFrame(frameNo)
		temp = self.data


		#import signal processing modules and generate Numpy array
		from scipy.signal import hilbert
		from numpy import log10
		bMode = log10(abs(hilbert(temp, axis = 0)))
		bMode = bMode - bMode.max()

		#import matplotlib and create plot
		import matplotlib.pyplot as plt

		if not self.isSectorData:
			import matplotlib.cm as cm
			im = plt.imshow(bMode, cmap = cm.gray, vmin = -3, vmax = 0, extent = [0, self.fovX,  self.fovY, 0])
			plt.show()

		else:
			fig = plt.figure()
			ax = fig.add_subplot(111)
			#ax.invert_yaxis()
			bMode[bMode < -3] = -3	
			ax.pcolormesh(self.X, self.Y, bMode, cmap = 'bone')
			ax.invert_yaxis()
			plt.show()
			

	def setRoiFixedSize(self, windowX = 4, windowY = 4):
		'''The Roi size here is given in mm'''

		
		from matplotlib import pyplot
		from scipy.signal import hilbert
		from numpy import log10, arange
		
		self.readFrame()
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


	def setRoi(self):
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
		self.readFrame(0)
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

	def createRoiArray(self):
#		#Check that a Data set has been loaded	
		#could be image sequence or just a 2-D image
		self.readFrame(0)
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

	def showRoiImage(self):
		bMode = self.createRoiArray()	
		#import matplotlib and create plot
		import matplotlib.pyplot as plt
		import matplotlib.cm as cm
		palette = cm.gray
		palette.set_bad('r')

		im = plt.imshow(bMode, cmap = palette, extent = [0, self.fovX, self.fovY, 0])
		plt.show()



	def saveRoiImage(self, fname):
		bMode = self.createRoiArray()	
		#import matplotlib and create plot
		import matplotlib.pyplot as plt
		import matplotlib.cm as cm
		palette = cm.gray
		palette.set_bad('r')
	
		im = plt.imshow(bMode, cmap = palette, extent = [0, self.fovX, self.fovY, 0])
		plt.savefig(fname)
		plt.close()

	def createStrainImage(self, preFrameNo, postFrameNo):
		'''Takes two frame numbers as input and produces displacement quality and strain images for plotting'''
		self.preFrame = preFrameNo
		from numpy import arange, meshgrid, zeros
		from scipy import interpolate
		from blockMatch import blockMatchClass
	
		self.readFrame(preFrameNo)
		pre = self.data.copy()	
		self.readFrame(postFrameNo)
		post = self.data.copy()	
		blockMatchInstance = blockMatchClass(pre, post, self.params)
		
		blockMatchInstance.trackSeeds()
		blockMatchInstance.trackNonSeeds()
		blockMatchInstance.displacementToStrain()

		self.startY = blockMatchInstance.startY
		self.stopY = blockMatchInstance.stopY
		self.startX = blockMatchInstance.startX
		self.stopX = blockMatchInstance.stopX

		self.startYStrain = blockMatchInstance.startYstrain 
		self.stopYStrain = blockMatchInstance.stopYstrain
		self.stepY = blockMatchInstance.stepY

		strainRfIndexes = arange(self.startYStrain, self.stopYStrain, self.stepY)
		strainRfIndexesNew = arange(self.startYStrain, self.stopYStrain - self.stepY)
		
		self.strain = zeros( (self.data.shape[0], self.data.shape[1] ) )

		'''Upsample strain image to same spacing as B-mode image.'''
		for x in range(self.startX, self.stopX):	
			interp = interpolate.interp1d(strainRfIndexes, blockMatchInstance.strain[:,x - self.startX] )
			self.strain[strainRfIndexesNew, x] = interp(strainRfIndexesNew )

	
		'''Convert array containing strain values to RGBALpha array'''
		from matplotlib import cm

		palette = cm.ScalarMappable()
		palette.set_cmap('gray')
		self.strain = abs(self.strain)
		self.strain[self.strain > self.strainLimit ] = self.strainLimit
		self.strain = palette.to_rgba(self.strain)

		del blockMatchInstance, post			
		'''Create B-mode image'''	
		from scipy.signal import hilbert
		from numpy import log10
		bMode = log10(abs(hilbert(pre, axis = 0)))
		bMode = bMode - bMode.max()
		bMode[bMode < -3]  = -3
		palette = cm.ScalarMappable()
		palette.set_cmap('gray')

		bMode = palette.to_rgba(bMode)

		if not rfClass.roiX == None:		
			self.strain[0:min(rfClass.roiY), :, :] = bMode[0:min(rfClass.roiY), :, :]
			self.strain[max(rfClass.roiY):-1, :, :] = bMode[max(rfClass.roiY):-1, :, :]
			self.strain[:, 0:min(rfClass.roiX), :] = bMode[:, 0:min(rfClass.roiX), :]
			self.strain[:, max(rfClass.roiX):-1, :] = bMode[:, max(rfClass.roiX):-1, :]
					
		self.fovStrainX = self.deltaX*self.strain.shape[1]
		self.fovStrainY = self.deltaY*self.strain.shape[0] 


	def deleteStrainImage(self):
		if not self.strain == None:
			self.strain = None


	def plotBmodeStrain(self):
		'''Plot b-mode and strain images together, saves them to output directory.'''

		self.readFrame(self.preFrame)	
		temp = self.data
		#import signal processing modules and generate Numpy array
		from scipy.signal import hilbert
		from numpy import log10
		bMode = log10(abs(hilbert(temp, axis = 0)))
		bMode = bMode - bMode.max()

		#import matplotlib and create plot
		import matplotlib.pyplot as plt

		fig = plt.figure()
		ax = fig.add_subplot(121)
		import matplotlib.cm as cm
		plt.imshow(bMode, cmap = cm.gray, vmin = -3, vmax = 0, extent = [0, self.fovX,  self.fovY, 0])
		ax2 = fig.add_subplot(122)	
	
		
		plt.imshow(self.strain, extent = [0, self.fovStrainX, self.fovStrainY, 0] )
		plt.show()

	def writeBmodeStrain(self, fname):
		self.readFrame(self.preFrame)
		temp = self.data

		#import signal processing modules and generate Numpy array
		from scipy.signal import hilbert
		from numpy import log10
		bMode = log10(abs(hilbert(temp, axis = 0)))
		bMode = bMode - bMode.max()

		#import matplotlib and create plot
		import matplotlib.pyplot as plt
		fig = plt.figure()
		ax = fig.add_subplot(121)
		import matplotlib.cm as cm
		plt.imshow(bMode, cmap = cm.gray, vmin = -3, vmax = 0, extent = [0, self.fovX,  self.fovY, 0])
		ax2 = fig.add_subplot(122)	
		
		plt.imshow(self.strain, extent = [0, self.fovStrainX, self.fovStrainY, 0] )
		plt.savefig(fname)
		plt.close()

	def createWriteStrainItk(self, preFrameNo, fname):
		'''This function calculates strain for a fiven frame number and writes the file to fname + .mhd'''
		import itk, numpy
	        converter = itk.PyBuffer[itk.Image.F2]	
		writerType = itk.ImageFileWriter[itk.Image.F2]
		writer = writerType.New()
		from numpy import arange, meshgrid, zeros
		from scipy import interpolate
		from blockMatch import blockMatchClass
	
		self.readFrame(preFrameNo)
		pre = self.data.copy()	
		postFrameNo = preFrameNo + 1
		self.readFrame(postFrameNo)
		post = self.data.copy()	
		blockMatchInstance = blockMatchClass(pre, post, self.params)
		
		blockMatchInstance.trackSeeds()
		blockMatchInstance.trackNonSeeds()
		blockMatchInstance.displacementToStrain()

		self.startY = blockMatchInstance.startY
		self.stopY = blockMatchInstance.stopY
		self.startX = blockMatchInstance.startX
		self.stopX = blockMatchInstance.stopX

		self.startYStrain = blockMatchInstance.startYstrain 
		self.stopYStrain = blockMatchInstance.stopYstrain
		self.stepY = blockMatchInstance.stepY

		strainRfIndexes = arange(self.startYStrain, self.stopYStrain, self.stepY)
		strainRfIndexesNew = arange(self.startYStrain, self.stopYStrain - self.stepY)
		
		self.strain = zeros( (self.data.shape[0], self.data.shape[1] ) )

		'''Upsample strain image to same spacing as B-mode image.'''
		for x in range(self.startX, self.stopX):	
			interp = interpolate.interp1d(strainRfIndexes, blockMatchInstance.strain[:,x - self.startX] )
			self.strain[strainRfIndexesNew, x] = interp(strainRfIndexesNew )

	
		if not rfClass.roiX == None:		
			strain = numpy.float32( abs(self.strain[ min(rfClass.roiY):max(rfClass.roiY), min(rfClass.roiX):max(rfClass.roiX) ] ) )	
		else:	
			strain = numpy.float32(abs(self.strain)	)		
	
		strain[strain > self.strainLimit] = self.strainLimit	
		im = converter.GetImageFromArray(strain)	
		im.SetSpacing( (self.deltaX, self.deltaY) )
		writer.SetFileName(fname)
		writer.SetInput( im )
		writer.Update()
