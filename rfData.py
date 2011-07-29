
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

	def readSeimensRFD(self, filename):
		"""Load functions from a C library for reading files from the Seimens Antares, created by Matt McCormick."""
		import readrfd
		self.data = readrfd.readrfd(filename)
		hInfo = readrfd.headerInfo(filename)
		self.fs = hInfo[1]
		self.deltaX = (1/hInfo[2])*10.
		self.deltaY = self.soundSpeed/(2*self.fs)*10**3
		self.fovX = self.deltaX*(self.data.shape[1] - 1)
		self.fovY = self.deltaY*(self.data.shape[0] -1 )	

	def readSeimensEIFiles(self, dirname):
		import os, numpy
		startdir = os.getcwd()
		os.chdir(dirname)
		
		metaData = open('metaData.txt', 'r')
		metaList = [0,0,0,0]
		for i in range(4):
			metaList[i] = int(metaData.readline() )
	
		metaData.close()
		nFrames = len( [name for name in os.listdir('.') if os.path.isfile(name)])
		nFrames -= 1 #because of metaData.txt
		points = len( range( metaList[3]/2, metaList[2] + metaList[3]/2 ) )
		aLines = metaList[1]
		from numpy import zeros, fromfile
		self.data = zeros( (points, aLines, nFrames) )
		readThisMany = (metaList[2]+ metaList[3]/2)*metaList[1]
		for x in range(nFrames):
			tempIn = open(str(x), 'r')
			frame = fromfile(tempIn,  numpy.int16, readThisMany)
			frame = frame.reshape( ( metaList[2]+ metaList[3]/2, metaList[1] )  , order = 'F' )
			self.data[:,:,x] = frame[metaList[3]/2:metaList[3]/2 + metaList[2], :].copy('C')

		self.soundSpeed = 1540.  #m/s
		self.fs = 40.*10**6  #Hz
		self.deltaX = 40./aLines  #assume a 40 mm FOV
		self.deltaY = self.soundSpeed/(2*self.fs)*10**3
		self.fovX = self.deltaX*(self.data.shape[1] - 1)
		self.fovY = self.deltaY*(self.data.shape[0] -1 )	

		os.chdir(startdir)

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
		if self.data == None:
			return	
		
		
		#could be image sequence or just a 2-D image
		if len(self.data.shape) > 2:
			temp = self.data[:,:,frameNo]

		else:
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
		if self.data == None:
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
		if len(self.data.shape) > 2:
			temp = self.data[:,:,0]

		else:
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
		if len(self.data.shape) > 2:
			temp = self.data[:,:,0]

		else:
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



	def saveROIImage(self, fname):
		bMode = self.createRoiArray()	
		#import matplotlib and create plot
		import matplotlib.pyplot as plt
		import matplotlib.cm as cm
		palette = cm.gray
		palette.set_bad('r')
	
		im = plt.imshow(bMode, cmap = palette, extent = [0, self.fovX, self.fovY, 0])
		plt.savefig(fname)

