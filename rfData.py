class rfClass:
	
	def __init__(self):
		"""  A class for managing B-mode and strain imaging of RF data """
		self.soundSpeed = 1540. #m/s
		self.fs = 40.*10**6 #Default sampling frequency in Hz
		self.deltaY = self.soundSpeed/(2*self.fs)*10**3 #Pixel spacing in mm
		self.data = None
		self.roiX = None
		self.roiY = None
		self.isSectorData = False
		self.minimumROIBoundaryY = 0  #for processing strain data
		self.minimumROIBoundaryX = 0

		#set parameters for block matching
		from blockMatch import blockMatchingParameters
		self.params = blockMatchingParameters()
		self.strainLimit = .02
		self.strain = None

		#parameters for file reading
		self.fname = None
		self.ftype = None

	def setInputData(self, filename, dataType):
		'''Input:  filename:  The name of the file or directory the data is from \n
		   dataType:  A lowercase string indicating the type of the data, choices are: \n
		   		ei
				rfd
				'''

		if not( dataType == 'ei' or dataType == 'rfd' ):
			print 'Error.  Datatype must be ei or rfd'
			return

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
				
	
		
			
	def readFrame(self, frameNo):
		'''When calling this from python it will get frames numbered from 0 to numFrames -1'''
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
			

	def readUltrasonixData(self, filename):
		"""Import data from the wobbler transducer on the Ultrasonix SonixTouch ultrasound system, output from a GUI created by 
		myself."""
		from readUltrasonixFiles import readUltrasonixData
		from numpy import zeros, sin, cos, pi
		(self.data, header) = readUltrasonixData(filename)
		self.fs = 40.*10**6 #Hz, always 40 Mhz
		self.deltaTheta = .6088  #degrees
		self.rt = 42.9  #distance from array foci to image
		self.rm = 27.25 #distance to motor rotation axis, mm
		self.deltaR = self.soundSpeed/(2.*self.fs)*10**3
		self.isSectorData = True

		self.X = zeros((self.data.shape[0],self.data.shape[1]) )
		self.Y = zeros((self.data.shape[0],self.data.shape[1]) )

		centerLine = self.data.shape[1]/2


		for x in range(self.data.shape[1]):
				startX = self.rt*sin( (self.deltaTheta * (x - centerLine) )*pi/180.0 )
				startY = self.rt*cos( (self.deltaTheta * (x - centerLine) )*pi/180.0 )
				deltaX = self.deltaR*sin( (self.deltaTheta * (x -centerLine) )*pi/180.0 )
				deltaY = self.deltaR*cos( (self.deltaTheta * (x - centerLine) )*pi/180.0 )


				for y in range(self.data.shape[0]): 
					self.X[y,x] = startX + deltaX*y 
					self.Y[y,x] = startY + deltaY*y 
				


	
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
			self.roiX = [0,0]
			self.roiY = [0,0]
			#lower boundary, make sure it doesn't violate block matching constraints
			tempX	= int(erelease.xdata/self.deltaX)
			if tempX < self.minimumROIBoundaryX:
				tempX = self.minimumROIBoundaryX
			if tempX > self.data.shape[1] - self.minimumROIBoundaryX - 1:
				tempX = self.data.shape[1] -self.minimumROIBoundaryX - 1			
			self.roiX[0] = tempX
			
			tempX	= int(eclick.xdata/self.deltaX)
			if tempX < self.minimumROIBoundaryX:
				tempX = self.minimumROIBoundaryX
			if tempX > self.data.shape[1] - self.minimumROIBoundaryX - 1:
				tempX = self.data.shape[1] -self.minimumROIBoundaryX - 1			
			self.roiX[1] = tempX
			
			tempY	= int(eclick.ydata/self.deltaY)
			if tempY < self.minimumROIBoundaryY:
				tempY = self.minimumROIBoundaryY
			if tempY > self.data.shape[0] - self.minimumROIBoundaryY - 1:
				tempY = self.data.shape[0] -self.minimumROIBoundaryY - 1			
			self.roiY[0] = tempY

			tempY	= int(erelease.ydata/self.deltaY)
			if tempY < self.minimumROIBoundaryY:
				tempY = self.minimumROIBoundaryY
			if tempY > self.data.shape[0] - self.minimumROIBoundaryY - 1:
				tempY = self.data.shape[0] -self.minimumROIBoundaryY - 1			
			self.roiY[1] = tempY


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
		xLow = min(self.roiX); xHigh = max(self.roiX)
		yLow = min(self.roiY); yHigh = max(self.roiY)
		
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

		if not self.roiX == None:		
			self.strain[0:min(self.roiY), :, :] = bMode[0:min(self.roiY), :, :]
			self.strain[max(self.roiY):-1, :, :] = bMode[max(self.roiY):-1, :, :]
			self.strain[:, 0:min(self.roiX), :] = bMode[:, 0:min(self.roiX), :]
			self.strain[:, max(self.roiX):-1, :] = bMode[:, max(self.roiX):-1, :]
					
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

	
		if not self.roiX == None:		
			strain = numpy.float32( abs(self.strain[ min(self.roiY):max(self.roiY), min(self.roiX):max(self.roiX) ] ) )	
		else:	
			strain = numpy.float32(abs(self.strain)	)		
	
		strain[strain > self.strainLimit] = self.strainLimit	
		im = converter.GetImageFromArray(strain)	
		im.SetSpacing( (self.deltaX, self.deltaY) )
		writer.SetFileName(fname)
		writer.SetInput( im )
		writer.Update()
