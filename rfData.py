class rfClass:
	
	def __init__(self):
		"""  A class for managing B-mode and strain imaging of RF data """
		self.soundSpeed = 1540. #m/s
		self.fs = 40.*10**6 #Default sampling frequency in Hz
		self.deltaY = self.soundSpeed/(2*self.fs)*10**3 #Pixel spacing in mm
		self.data = None
		self.isSectorData = False


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
		import types
		if self.data.__class__ == types.NoneType:
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
			


	def setROI(self, xLow, xHigh, yLow, yHigh):
		"""Input boundaries in array Indexes.  This function will highlight the ROI on the image."""
		
		#Check that a Data set has been loaded	
		import types
		if self.data.__class__ == types.NoneType:
			return	
		

		#Assume it is a square data set for the moment.
		self.roi_xLow = xLow 
		self.roi_xHigh = xHigh 
		self.roi_yLow = yLow
		self.roi_yHigh = yHigh
	
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

		#import matplotlib and create plot
		import matplotlib.pyplot as plt

		import matplotlib.cm as cm
		palette = cm.gray
		palette.set_bad('r')	
	
		from numpy import ma, zeros
		bMode = ma.array(bMode)
		bMode.mask = zeros(bMode.shape)
		#replace some B-mode data with sentinels
		
		bMode.mask[yLow:yHigh, xLow] = True
		bMode.mask[yLow:yHigh, xHigh] = True
		if yLow < 10:
			yLow = 10		
		bMode.mask[yLow-10:yLow, xLow:xHigh] = True
		
		if yHigh > bMode.mask.shape[0] - 10:
			yHigh = bMode.mask.shape[0] - 10
		bMode.mask[yHigh-10:yHigh, xLow:xHigh] = True

		im = plt.imshow(bMode, cmap = palette,  extent = [0, self.fovX,  self.fovY, 0])
		plt.show()


