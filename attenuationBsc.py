from rfData import rfClass
from faranScattering import faranBsc

class attenuation(rfClass):

	def __init__(self, sampleName, refName, sampleType, refType = None, refAttenuation = .5):

		if refType == None:
			refType = sampleType

		super(attenuation, self).__init__(sampleName, sampleType)
		self.refRf = rfClass(refName, refType)
			
		#Check to see that reference data and sample data contain
		#the same number of points
		if self.points != self.refRf.points or self.lines != self.refRf.lines:
			print "Error.  Sample and reference images must be the same size. \n\ "
			return

		#make this number odd
		self.lsqFitPoints = 13
		self.halfLsq = self.lsqFitPoints//2
		self.betaRef = refAttenuation
		#get window sizes and overlap
		self.windowYmm = 4
		self.windowXmm = 4
		self.overlapY = .75
		self.overlapX = 1
		
		self.windowX =int( self.windowXmm/self.deltaX)
		self.windowY =int( self.windowYmm/self.deltaY)
		self.windowY -= self.windowY%4
		
		#make the windows odd numbers
		if not self.windowY%2:
			self.windowY +=1

		if not self.windowX%2:
			self.windowX +=1

		self.halfY = self.windowY//2
		self.halfX = self.windowX//2
			
		#overlap the windows axially by 50%
		stepY = int(  (1-self.overlapY)*self.windowY )
		startY = self.halfY
		stopY = self.points - self.halfY - 1
		self.winCenterYpriorLsq = range(startY, stopY, stepY)
		
		#cutoff some more points because of least squares fitting	
		self.winCenterY= self.winCenterYpriorLsq[self.halfLsq:-self.halfLsq]	
	
		stepX = int( (1-self.overlapX)*self.windowX )
		if stepX < 1:
			stepX = 1	
		startX =self.halfX
		stopX = self.lines - self.halfX - 1
		self.winCenterX = range(startX, stopX, stepX)

		##Within each window a Welch-Bartlett style spectrum will be estimated		
		##Figure out the number of points used in an individual FFT based on
		##a 50% overlap and rounding the window size to be divisible by 4
		self.bartlettY = self.windowY//2
		self.spectrumFreqStep = self.fs/self.barlettY

	def CalculateSpectrumWelch(self, region):
		'''Return the power spectrum of a region based on a Welch-Bartlett method.
		The window used in each FFT is half the length of the total window.
		The step size is half the size of the FFT window.
		Average over A-lines.

		This function assumes the size of the region is divisible by 4.
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
			fourierData =numpy.fft.fft(dataWindow, axis = 0)
			fftList.append( fourierData.copy() )
	
		fftSample = numpy.zeros( fourierData.shape)
		for f in range(3):
			fftSample += abs(fftList[f])

		fftSample = fftSample.mean(axis = 1)
		return fftSample	
			
	def FindCenterFrequency(self, centerFreq = None):
		'''Fit the transmit pulse to a Gaussian. Show a B-mode image of the reference phantom
		and pick the focus by visual inspection.'''

		import numpy
		from scipy import optimize

		self.refRf.SetRoiFixedSize(self.windowYmm, self.windowXmm)

		#get the depth of the ROI
		focalDepth = (self.refRf.roiY[0] + self.refRf.roiY[1])//2
		
		#work out the Roi size
		tempRoi =self.refRf.data[focalDepth - self.halfY:focalDepth + self.halfY + 1, self.winCenterX[0] - self.halfX:self.winCenterX[0] + self.halfX + 1]
		points = tempRoi.shape[0]
		points -= points%4
		points /= 2
		
		deltaF = self.refRf.fs/(points*10**6)
		freq = numpy.arange(0, points/2*deltaF, deltaF)
		
		#ask for an initial guess on the center frequency
		from matplotlib import pyplot
		if not centerFreq:
			spec = self.CalculateSpectrumWelch(tempRoi)
			pyplot.plot(freq, spec[0:len(spec)/2])
			pyplot.show()
			centerFreq = input("Enter a center frequency (MHz) " )

		#create functions for LSQ gaussian fit
		fitfunc = lambda p,f: numpy.exp( - (f - p[0])**2 / (2*p[1]**2) )
		errfunc = lambda p,f, y: (fitfunc(p,f) - y )**2
		p0 = numpy.array( [float(centerFreq), .5] )


		#loop through depth calculating power spectrum performing Gaussian fit
		self.centerFreq = 0.
		self.bw = 0.
		count = 0
		for x in self.winCenterX:
			spec = self.CalculateSpectrumWelch(self.refRf.data[focalDepth - self.halfY:focalDepth + self.halfY + 1, x - self.halfX:x + self.halfX + 1] )
			spec /= spec.max()	
			#fit a Gaussian to the power spectrum
			[p1, success] = optimize.leastsq(errfunc, p0, args = (freq,spec[0:len(spec)/2]) )
			self.centerFreq += p1[0]
			self.bw += p1[1]
			count += 1

		self.bw /= count
		self.centerFreq /= count

		#show fitted spectrum
		self.gaussianFilter = fitfunc([self.centerFreq, self.bw],freq)
		pyplot.plot(freq, self.gaussianFilter )
		pyplot.show()
	
		#work out frequency indexes I'll use to work with spectrums
		self.spectrumLowCutoff = int( (self.centerFreq - self.bw)/self.spectrumFreqStep )
		self.spectrumHighCutoff	= int( (self.centerFreq + self.bw)/self.spectrumFreqStep )
	
		
	def CalculateAttenuationImage(self, convertToRgb = True):
		'''Loop through the image and calculate the spectral shift at each depth.
		Perform the operation 1 A-line at a time to avoid repeating calculations.
		Input:
		convertToRgb:  A switch to make the output an RGB image that I can plot directly, but
		I'll lose the attenuation slope values.
		'''

		import numpy
		self.ReadFrame()
		self.refRf.ReadFrame()
		self.FindCenterFrequency()	
		numY = len(self.winCenterY)
		numX = len(self.winCenterX)
		self.attenImage = numpy.zeros( (numY, numX) )	
		startY = self.winCenterY[self.lsqFitPoints//2] 
		startX = self.winCenterX[0]
		stepY = self.winCenterY[1] - self.winCenterY[0]
		stepX = self.winCenterX[1] - self.winCenterX[0]
			
		for x in range(numX):
			if not x:
				from time import time
				t1 = time()
			tempRegionSample = self.data[:, self.winCenterX[x] - self.halfX:self.winCenterX[x] + self.halfX + 1]
			tempRegionRef = self.refRf.data[:, self.winCenterX[x] - self.halfX:self.winCenterX[x] + self.halfX + 1]
			self.attenImage[:, x] = self.CalculateAttenuationAlongBeamLine(tempRegionSample, tempRegionRef)
			if not x:
				t2 = time()
				print "Time for a beamline was" + str(t2 - t1) + " seconds"
				print "Time for all lines is: " + str( (t2-t1)*numX/60 ) + "minutes"
					
		#convert slope value to attenuation value
		self.attenImage *= -8.686/(4*self.bw)
		self.attenImage += self.betaRef
			
		if convertToRgb:
			self.attenImage = self.CreateParametricImage(self.attenImage,[startY, startX], [stepY, stepX] )



	def CalculateAttenuationAlongBeamLine(self, sampleRegion, refRegion):
		'''Input:
		sampleRegion:  RF data from the sample
		refRegion:  RF data from the reference phantom

		'''
		import numpy
		from matplotlib import pyplot
		from scipy.signal import hamming
		import cv
		
		#find out number of points in list
		points = 2*self.halfY
		points -= points%4
		points /= 2
		
		#compute 3 fourier transforms and average them	
		windowFunc = hamming(points).reshape(points,1)

		spectrumList = []

		for y in self.winCenterYpriorLsq:
			#get first spectrum for x-correlation
			maxDataWindow = sampleRegion[y - points:y + points, :]
			fftListSample = []
			for f in range(3):
				dataWindow = maxDataWindow[(points/2)*f:(points/2)*f + points, :]*windowFunc	
				fourierData = numpy.fft.fft(dataWindow, axis = 0)
				fftListSample.append( fourierData.copy() )
		
			fftSample = numpy.zeros( fftListSample[0].shape)	
			for f in range(3):
				fftSample += abs(fftListSample[f])

			fftSample = fftSample.mean(axis = 1)
				
			#######REFERENCE REGION#######
			maxDataWindow = refRegion[y - points:y + points, :]

			#compute 3 fourier transforms and average them	
			fftListRef = []
			for f in range(3):
				dataWindow = maxDataWindow[(points/2)*f:(points/2)*f + points, :]*windowFunc	
				fourierData =numpy.fft.fft(dataWindow, axis = 0)
				fftListRef.append(fourierData.copy() )


			fftRef = numpy.zeros( fftListRef[0].shape)	
			for f in range(3):
				fftRef += abs(fftListRef[f])

			fftRef = fftRef.mean(axis = 1)
		
			###Divide sample spectrum by reference spectrum to perform deconvolution
			#gfr = gaussian filtered ratio
			gfr = fftSample[0:len(fftSample)//2]/fftRef[0:len(fftSample)//2]*self.gaussianFilter	
			spectrumList.append(gfr.copy())

		#work out number of freq points equal to .5 MHz
		deltaF = ((self.fs/2)/10**6)/len(gfr)
		pointsToCut = int( 1/deltaF)

		shift = [0]*(len(self.winCenterY) - 1)

		slope = numpy.zeros( len(self.winCenterY) - self.lsqFitPoints)
		#compute spectral cross correlation over 4 windows and perform
		#least squares fitting
		resultCv = numpy.zeros( (2*pointsToCut + 1, 1 ) )  
		resultCv = cv.fromarray(numpy.float32(resultCv) )

		#Compute spectral shift at adjacent points
		for y in range(len(self.winCenterY) - 1):
			gfr0 = spectrumList[y]
			gfr0[gfr0 == numpy.nan] = 0
			gfr0 = gfr0.reshape( (len(gfr0), 1) )
			image = cv.fromarray(numpy.float32( gfr0 ) )
			gfr1 = spectrumList[y + 1]
			gfr1[gfr1 == numpy.nan] = 0 
			gfr1 = gfr1[pointsToCut: -pointsToCut]
			gfr1 = gfr1.reshape( (len(gfr1), 1) )
			template = cv.fromarray( numpy.float32( gfr1) )	
			cv.MatchTemplate(template, image, resultCv, cv.CV_TM_CCORR_NORMED)
			resultNp = numpy.asarray(resultCv)
			maxShift = resultNp.argmax()
			maxShiftFreq = (maxShift - pointsToCut)*deltaF
			shift[y] = maxShiftFreq 
	
		for y in range(len(self.winCenterY) - self.lsqFitPoints): 	
			#perform linear least squares fit to line
			#want to solve equation y = mx + b for m
			#[x1   1      [m     = [y1
			# x2   1       b ]	y2
			# x3   1]               y3]
			#
			#strain image will be smaller than displacement image
			A = numpy.ones( (self.lsqFitPoints,2) )
			A[:,0] = numpy.arange(0, self.lsqFitPoints)*deltaF
			b = numpy.array(shift[y: y + self.lsqFitPoints])
			out = numpy.linalg.lstsq(A, b)
			slope[y] = out[0][0]

		return slope
