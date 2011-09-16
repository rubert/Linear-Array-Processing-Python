from rfData import rfClass
from faranScattering import faranBsc

class attenuation(rfClass):

	def __init__(self, sampleName, refName, sampleType, refType = None, refAttenuation = .5, freqLow = 2., freqHigh = 8., attenuationKernelSizeYmm = 15, blockYmm = 8, blockXmm = 8, overlapY = .85, overlapX = .85, spectralShiftRangeMhz):
		'''Description:
			This class implements the reference phantom method of Yao et al.  It inherits from the RF data class
			defined for working with simulations and Seimens rfd files.
		
		Input:
			sampleName:  Filename of sample RF data
			refName:  Filename of reference RF data
			sampleType:  The file type of the sample RFdata
			refTYpe:  The file type of the reference RF data.  Defaults to the sampleType
			refAttenuation:  Assuming a linear dependence of attenuation on frequency,
					the attenuation slope in dB/(cm MHz)
			freqLow:  The low frequency for the CZT in MHz
			freqHigh: The high frequency for the CZT in MHz
			attenuationKernelSizeMM:  The size of the data segment used to do the least squares fitting
			to find center frequency shift with depth.
			blockSizeY,Xmm:
			blockOverlap:

		  Throughout the code I'll call the 1-D segment where a single FFT is performed a window

		  A block will refer to a 2-D region spanning multiple windows axially and several A-lines
		  where a power spectrum is calculated by averaging FFTs 
					'''
		import chirpz
		import numpy
		
		if not refType:
			refType = sampleType

		super(attenuation, self).__init__(sampleName, sampleType)
		self.refRf = rfClass(refName, refType)
			
		#Check to see that reference data and sample data contain
		#the same number of points
		if self.points != self.refRf.points or self.lines != self.refRf.lines:
			print "Error.  Sample and reference images must be the same size. \n\ "
			return

		#Attenuation estimation parameters
		self.betaRef = refAttenuation
		#get window sizes and overlap
		self.blockYmm = blockYmm
		self.blockXmm = blockXmm
		self.overlapY = overlapY
		self.overlapX = overlapX
		self.spectralShiftRangeMhz = spectralShiftRangeMHz 
		
		self.blockX =int( self.blockXmm/self.deltaX)
		self.blockY =int( self.blockYmm/self.deltaY)
		
		#make the block sizes in pixels odd numbers for the sake of calculating their centers
		if not self.blockY%2:
			self.blockY +=1

		if not self.blockX%2:
			self.blockX +=1

		self.halfY = self.blockY//2
		self.halfX = self.blockX//2
			
		#overlap the blocks axially by self.overlapY%
		stepY = int(  (1-self.overlapY)*self.blockY )
		startY = self.halfY
		stopY = self.points - self.halfY - 1
		self.blockCenterY = range(startY, stopY, stepY)
	
		#Set the attenuation kernel size in mm
		#Work it out in points
		#Make kernel size an odd number of poitns
		#Work out kernel size in mm
		self.attenuationKernelSizeYmm = attenuationKernelSizeYmm #attenuation estimation size used in least squares fit
		self.lsqFitPoints = int(self.attenuationKernelSizeYmm/(stepY*self.deltaY) ) #make this number odd
		if not self.lsqFitPoints%2:
			self.lsqFitPoints += 1
		self.halfLsq = self.lsqFitPoints//2
		self.attenuationKernelSizeYmm = self.lsqFitPoints*stepY*self.deltaY
		
			
		#cutoff some more points because of least squares fitting, cross correlation	
		self.attenCenterY= self.blockCenterY[self.halfLsq + 1:-(self.halfLsq + 1)]	
	
		stepX = int( (1-self.overlapX)*self.blockX )
		if stepX < 1:
			stepX = 1	
		startX =self.halfX
		stopX = self.lines - self.halfX
		self.blockCenterX = range(startX, stopX, stepX)

		##Within each block a Welch-Bartlett style spectrum will be estimated		
		##Figure out the number of points used in an individual FFT based on
		##a 50% overlap and rounding the block size to be divisible by 4
		self.blockY -= self.blockY%4
		self.bartlettY = self.blockY//2
		self.spectrumFreqStep = (freqHigh - freqLow)/self.bartlettY
		self.spectrumFreq = numpy.arange(0, self.bartlettY)*self.spectrumFreqStep + freqLow

		#set-up the chirpZ transform
		self.czt = chirpz.ZoomFFT(self.bartlettY, freqLow, freqHigh, self.bartlettY, self.fs/10.**6)
		self.spectralShiftRangePoints = int(self.spectralShiftRangeMhz/self.spectrumFreqStep)
			
	
	
	def CalculateAttenuationImage(self,convertToRgb = True ):
		'''Estimate the center frequency by fitting to a Gaussian'''
		'''Loop through the image and calculate the spectral shift at each depth.
		Perform the operation 1 A-line at a time to avoid repeating calculations.
		Input:
		convertToRgb:  A switch to make the output an RGB image that I can plot directly, but
		I'll lose the attenuation slope values.
		'''

		import numpy
		import types
		if type(self.data) == types.NoneType:
			self.ReadFrame()
			self.refRf.ReadFrame()
		
		self.FindCenterFrequency()	
		numY = len(self.attenCenterY)
		numX = len(self.blockCenterX)
		self.attenuationImage = numpy.zeros( (numY, numX) )	
		startY = self.attenCenterY[0] 
		startX = self.blockCenterX[0]
		stepY = self.blockCenterY[1] - self.blockCenterY[0]
		stepX = self.blockCenterX[1] - self.blockCenterX[0]

		#first compute the power spectrum at each depth for the reference phantom and
		#average over all the blocks
		self.computeReferenceSpectrum()
		self.computeSampleSpectrum()
	
	
		#compute the log ratio
		#fit the log ratio at each depth to a line
		#to get derivative with respect to frequency
		for countY in range(len(refSpectrumList)):
			for countX in range(numX):
				logRatio = numpy.log( self.sampleSpectrum[:,countY,countX]/self.refSpectrum[:,countY, countX] )
				dFreqLogRatio[countY, countX] = self.lsqFit(logRatio)/self.spectrumFreqStep
				
		attenKernelCm = self.attenuationKernelSizeYmm/10.
		#now compute the derivative with respect to depth
		for countY in range(numY):
			for countX in range(numX):
				self.attenuationImage = self.lsqFit(dFreqLogRatio[countY:countY+self., countX] ) /attenKernelCm
					
					
		#convert slope value to attenuation value
		self.attenuationImage *= -8.686/4.
		self.attenuationImage += self.betaRef
		
		print "Mean attenuation value of: " + str( self.attenuationImage.mean() )
				
		if convertToRgb:
			self.attenuationImage = self.CreateParametricImage(self.attenuationImage,[startY, startX], [stepY, stepX] )

	def ComputeReferenceSpectrum(self):

		'''Calculate the spectrum of the reference region over every beamline at every
		depth, average over all the beamlines

		'''
		import numpy
	
		#list comprehension, avoids all list entries pointing to same array.  A little ugly	
		self.refSpectrum = numpy.zeros(self.bartlettY, len(self.blockCenterY))
		
		for countX,x in enumerate(self.blockCenterX):
			for countY,y in enumerate(self.blockCenterY):
				
				#######REFERENCE REGION#######
				maxDataWindow = refRegion[y - self.halfY:y + self.halfY+1, :]
				fftRef = self.CalculateSpectrumBlock(maxDataWindow)
				self.refSpectrum[:,y] += fftRef

		#normalize
		for y in range(len(refSpectrumList)):
			self.refSpectrumList /= self.refSpectrumList[:,y].max()
		
	
			
	def ComputeSampleSpectrum(self):
		'''Calculate the spectra for the sample. 
		'''
		import numpy
		
		self.sampleSpectrum = numpy.zeros(self.bartlettY, len(self.blockCenterY), len(self.blockCenterX))
		
		for countX, x in enumerate(self.blockCenterX):
			for countY,y in enumerate(self.blockCenterY):
				maxDataWindow = sampleRegion[y - self.halfY:y + self.halfY+1, x - self.halfX: x + self.halfX + 1]
				fftSample = self.CalculateSpectrumBlock(maxDataWindow)
				self.sampleSpectrum[:,countY,countX] = fftSample
					

		
	
	def lsqFit(self, inputArray):
		'''
		Input:  
		inputArray: An array containing frequency shifts over depth.  The units are
		MHz
	
		Output:
		slope: A scalar. The value of the slope of the least squares line fit.
		slope is in units of Mhz/cm
			
		#perform linear least squares fit to line
		#want to solve equation y = mx + b for m
		#[x1   1      [m     = [y1
		# x2   1       b ]	y2
		# x3   1]               y3]
		#
		'''
		import numpy
		pointStep = self.blockCenterY[1] - self.blockCenterY[0]
		#in cm
		deltaZ = 1540./(2*self.fs)*10**2*pointStep
		
		A = numpy.ones( (self.lsqFitPoints,2) )
		A[:,0] = numpy.arange(0, self.lsqFitPoints)*deltaZ
		b = inputArray
		out = numpy.linalg.lstsq(A, b)
	
		return out[0][0]

	def CalculateSpectrumBlock(self, region):
		'''Return the power spectrum of a region based on a Welch-Bartlett method.
		The block used in each FFT is half the length of the total window.
		The step size is half the size of the FFT window.
		Average over A-lines.

		This function assumes the size of the region is divisible by 4.

		It uses a zoomed in FFT to compute the power spectrum.  The zoomed in FFT is given by the
		chirpz transform.
		'''
		from scipy.signal import hamming
		import numpy
		points = region.shape[0]
		points -= points%4
		points /= 2
		#######SAMPLE REGION#############
		maxDataWindow = region[0:2*points, :]
		
		#compute 3 fourier transforms and average them	
		windowFunc = hamming(points).reshape(points,1)
		fftList = []
		for f in range(3):
			dataWindow = maxDataWindow[(points/2)*f:(points/2)*f + points, :]*windowFunc	
			fourierData = self.czt(dataWindow, axis = 0).copy()
			fftList.append( fourierData.copy() )
	
		fftSample = numpy.zeros( fourierData.shape)
		for f in range(3):
			fftSample += abs(fftList[f])

		fftSample = fftSample.mean(axis = 1)
		return fftSample	
