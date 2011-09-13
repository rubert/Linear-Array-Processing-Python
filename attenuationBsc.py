from rfData import rfClass
from faranScattering import faranBsc

class attenuation(rfClass):

	def __init__(self, sampleName, refName, sampleType, refType = None, refAttenuation = .5, freqLow = 2., freqHigh = 8.):
		'''Input:
			sampleName:  Filename of sample RF data
			refName:  Filename of reference RF data
			sampleType:  The file type of the sample RFdata
			refTYpe:  The file type of the reference RF data.  Defaults to the sampleType
			refAttenuation:  Assuming a linear dependence of attenuation on frequency,
					the attenuation slope in dB/(cm MHz)
			freqLow:  The low frequency for the CZT in MHz
			freqHigh: The high frequency for the CZT in MHz

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
		self.blockYmm = 8
		self.blockXmm = 10
		self.overlapY = .9
		self.overlapX = .95
		self.spectralShiftRangeMhz = 2.
		
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
		self.attenuationKernelSizeYmm = 12 #attenuation estimation size used in least squares fit
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
			
	def FindCenterFrequency(self, centerFreq = None):
		'''Show a B-mode image of the reference phantom
		and pick the focus by visual inspection.  Then fit the transmit pulse
		to a Gaussian through a least-squares fit with 2 parameters,
		center frequency and variance.'''

		import numpy
		from scipy import optimize

		self.refRf.SetRoiFixedSize(self.blockXmm, self.blockYmm)

		#get the depth of the ROI
		focalDepth = (self.refRf.roiY[0] + self.refRf.roiY[1])//2
		
		
		#ask for an initial guess on the center frequency
		from matplotlib import pyplot
		tempRoi =self.refRf.data[focalDepth - self.halfY:focalDepth + self.halfY + 1, self.blockCenterX[0] - self.halfX:self.blockCenterX[0] + self.halfX + 1]
		if not centerFreq:
			spec = self.CalculateSpectrumBlock(tempRoi)
			pyplot.plot(self.spectrumFreq, spec)
			pyplot.show()
			centerFreq = input("Enter a center frequency (MHz) " )

		#create functions for LSQ gaussian fit
		fitfunc = lambda p,f: numpy.exp( - (f - p[0])**2 / (2*p[1]**2) )
		errfunc = lambda p,f, y: (fitfunc(p,f) - y )**2
		p0 = numpy.array( [float(centerFreq), 3.0] )


		#loop through bealines, calculating power spectrum, performing Gaussian fit
		self.centerFreq = 0.
		self.sigma = 0.
		count = 0
		for x in self.blockCenterX:
			spec = self.CalculateSpectrumBlock(self.refRf.data[focalDepth - self.halfY:focalDepth + self.halfY + 1, x - self.halfX:x + self.halfX + 1] )
			spec /= spec.max()	
			#fit a Gaussian to the power spectrum
			[p1, success] = optimize.leastsq(errfunc, p0, args = (self.spectrumFreq,spec) )
			self.centerFreq += p1[0]
			self.sigma += p1[1]
			count += 1

		self.sigma /= count
		self.centerFreq /= count
		#show fitted spectrum
		self.gaussianFilter = fitfunc([self.centerFreq, self.sigma],self.spectrumFreq)
		pyplot.plot(self.spectrumFreq, self.gaussianFilter )
		pyplot.show()
	
	
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

		for x in range(numX):
			if not x:
				from time import time
				t1 = time()
			tempRegionSample = self.data[:, self.blockCenterX[x] - self.halfX:self.blockCenterX[x] + self.halfX + 1]
			tempRegionRef = self.refRf.data[:, self.blockCenterX[x] - self.halfX:self.blockCenterX[x] + self.halfX + 1]
			self.attenuationImage[:, x] = self.CalculateAttenuationAlongBeamLine(tempRegionSample, tempRegionRef)
			if not x:
				t2 = time()
				print "Time for a beamline was" + str(t2 - t1) + " seconds"
				print "Time for all lines is: " + str( (t2-t1)*numX/60 ) + "minutes"
					
		#convert slope value to attenuation value
		self.attenuationImage *= -8.686/(4*self.sigma**2)
		self.attenuationImage += self.betaRef
		
		print "Mean attenuation value of: " + str( self.attenuationImage.mean() )
				
		if convertToRgb:
			self.attenuationImage = self.CreateParametricImage(self.attenuationImage,[startY, startX], [stepY, stepX] )

	def CalculateAttenuationAlongBeamLine(self, sampleRegion, refRegion):
		'''Calculate the attenuation along a single beam line using the hybrid method.
		Input:
		sampleRegion:  RF data from the sample, this is the entire column in the axial
		direction, and the block size in the lateral direction
		refRegion:  RF data from the reference phantom

		'''
		import numpy
		
		spectrumList = []
		for count,y in enumerate(self.blockCenterY):
			#get first spectrum for x-correlation
			maxDataWindow = sampleRegion[y - self.halfY:y + self.halfY+1, :]
			fftSample = self.CalculateSpectrumBlock(maxDataWindow)
				
			#######REFERENCE REGION#######
			maxDataWindow = refRegion[y - self.halfY:y + self.halfY+1, :]
			fftRef = self.CalculateSpectrumBlock(maxDataWindow)

			###Divide sample spectrum by reference spectrum to perform deconvolution
			#gfr = gaussian filtered ratio
			gfr = fftSample/fftRef*self.gaussianFilter	
			spectrumList.append(gfr.copy())

		

		#Compute centroid by fitting spectrum to 
		#a gaussian
		centroids = numpy.zeros(self.lsqFitPoints)	
		slope = numpy.zeros(len(self.attenCenterY))
		
		shift = [0]*self.lsqFitPoints
	
		
		def runningSum(iterable):
			sum = 0
			for x in iterable:
				sum += x
				yield sum
	
		#include contribution from attenuation terms				
		for y,centerY in enumerate(self.attenCenterY):
			if y == 0:
				betaDiffNepers = 0.0
			else:
				betaDiff= -8.686*slope[0:y].mean()/(4*self.sigma**2)
				betaDiffNepers = betaDiff/8.686
			
			depthBlockOne = self.blockCenterY[y]*self.deltaY/10.
			factorOne = numpy.exp(-4*betaDiffNepers*depthBlockOne*self.spectrumFreq)
				
			for w in range(self.lsqFitPoints):
				depthBlockTwo = (self.blockCenterY[y+w])*self.deltaY/10.
				factorTwo = numpy.exp(-4*betaDiffNepers*depthBlockTwo*self.spectrumFreq)
				shift[w] =self.findShift(spectrumList[y]*factorOne, spectrumList[y + w]*factorTwo )

			slope[y] = self.lsqFit(shift)	
		
		
		return slope

	def findShift(self, ps1, ps2):
		'''Use the openCV library to compute a cross-correlation between
		two power spectra and return the frequency shift'''
		import cv
		import numpy as np
		#allocate array to hold results of cross correlation
		resultNp = np.float32( np.zeros( (2*self.spectralShiftRangePoints + 1,1)  ) )
		resultCv = cv.fromarray(resultNp)

		#convert template and image to openCV friendly arrays
		template = cv.fromarray( np.float32( ps1.reshape(len(ps1), 1) ) )
		imageNp = np.zeros( (len(ps2) + 2*self.spectralShiftRangePoints , 1) )
		imageNp[self.spectralShiftRangePoints:self.spectralShiftRangePoints + len(ps2), 0] = ps2
		image = cv.fromarray( np.float32(imageNp))

		cv.MatchTemplate(template, image, resultCv, cv.CV_TM_CCORR_NORMED )
		resultNp = np.asarray(resultCv)

		#FIND MAXIMUM, PERFORM SUB-SAMPLE FITTING
		maxInd = resultNp.argmax()

		delta = 0.0
		if maxInd > 0 and maxInd < len(resultNp) - 1:
			c = resultNp[maxInd]
			a = (resultNp[maxInd-1] + resultNp[maxInd + 1] )/2 - c
			b = (resultNp[maxInd-1] - resultNp[maxInd + 1])/2
			delta = float( - b/ (2*a) )
		if abs(delta) > 1:
			delta = 0.0

		return (maxInd - self.spectralShiftRangePoints + delta)*self.spectrumFreqStep
	
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
