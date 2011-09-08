from rfData import rfClass
from attenuationBsc import attenuation
from faranScattering import faranBsc

class scattererSizeClass(attenuation):
			
	def ComputeBscImage(self, convertToRgb = True):
		'''First compute the attenuation image, then use the attenuation image and the
		reference phantom spectrum to get the spectrum at each point.  Minimize difference
		between this power spectrum with a spectrum calculated from a Gaussian autocorrelation function.'''
		
		import numpy

		self.CalculateAttenuationImage(convertToRgb = False)
		self.InitializeTheoreticalBackscatter()
		numY = len(self.winCenterY)
		numX = len(self.winCenterX)
		self.scatSizeImage = numpy.zeros( (numY, numX) )	
		startY = self.winCenterY[0] 
		startX = self.winCenterX[0]
		stepY = self.winCenterY[1] - self.winCenterY[0]
		stepX = self.winCenterX[1] - self.winCenterX[0]
		
		for yParam,yRf in enumerate(self.winCenterY):
			for xParam, xRf in enumerate(self.winCenterX):
				attenuationDifference = (self.deltaY*yRf)/10*self.attenImage[0:1+yParam,xParam].mean()	
				tempRegionSample = self.data[yRf - self.halfY:yRf + self.halfY + 1, xRf - self.halfX:xRf + self.halfX + 1]
				tempRegionRef = self.refRf.data[yRf - self.halfY:yRf + self.halfY + 1, xRf - self.halfX:xRf + self.halfX + 1]
				self.scatSizeImage[yParam, xParam] = self.CalculateScatSizeOnePoint(tempRegionSample, tempRegionRef, attenuationDifference)
			
	
		#convert scatterer size image to RGB parametric image	

		if convertToRgb:
			self.scatSizeImage = self.CreateParametricImage(self.scatSizeImage,[startY, startX], [stepY, stepX] )

	def InitializeTheoreticalBackscatter(self):
		'''Within the transducer bandwidth calculate the backscatter coefficients over a set of scatter size.
		Each array of backscatter coefficients is placed in a list whose different entries correpsond to different
		scatterer sizes.'''
		#for scatter sizes between 10 and 100 micrometers
		import numpy

		#I get a lot of divide by zeros when computing the theoretical backscatter
		#coefficients

		oldWarningSettings = numpy.seterr(all = 'ignore')
		instance = faranBsc()
		self.bscFaranSizes = numpy.arange(10,100,1)
		self.bscCurveFaranList =[]
		for d in self.bscFaranSizes:
			tempBsc = instance.calculateBSC(self.spectrumInTxdcerBandwidth, d)
			self.bscCurveFaranList.append( tempBsc.copy())

		numpy.seterr(**oldWarningSettings) 

	def CalculateScatSizeOnePoint(self, sampleRegion, refRegion, attenuationDifference):
		'''This function calculates power spectrum ratios and multiplies by an attenuation 
		difference with a known reference phantom.  The units on the attenuation difference are in 
		dB/(cm MHz)'''
		import numpy
		from matplotlib import pyplot
		from scipy.signal import hamming
		import cv
		
		attenuationDifference /= 8.686 #convert to Nepers
		#find out number of points in list
		points = 2*self.halfY
		points -= points%4
		points /= 2
		
		#compute 3 fourier transforms and average them	
		windowFunc = hamming(points).reshape(points,1)

		#get first spectrum for x-correlation
		maxDataWindow = sampleRegion[0:2*points]
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
		maxDataWindow = refRegion[0:2*points]

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
		
		#work out number of freq points equal to .5 MHz
		###Divide sample spectrum by reference spectrum to perform deconvolution
		rpmSpectrum = fftSample[0:len(fftSample)//2]/fftRef[0:len(fftSample)//2]
		
		deltaF = ((self.fs)/10**6)/(2*len(rpmSpectrum) )
		freq = numpy.arange(0, len(rpmSpectrum)*deltaF, deltaF)

		#work out attenuation difference by assu
		rpmSpectrum*=numpy.exp(-4*freq*attenuationDifference)
		diffs = self.ComputeBscCoefficients(rpmSpectrum[self.spectrumLowCutoff:self.spectrumHighCutoff] )		
		
		return self.bscFaranSizes[diffs.argmin()]	
		
	def ComputeBscCoefficients(self, BSCs):
		'''Compute a theoretical backscatter curve over the transducer bandwidth.  Find the 
		logarithm of the difference between the two curves. 
		 
		
		Input:  The frequencies at which the spectrum to be compared with theory is calculated,
			probably should not be outside of transducer bandwidth
			
		Output:  Error and corresponding scatterer sizes'''
		import numpy
		mmse = numpy.zeros( len(self.bscCurveFaranList) )
		for count,BSCt in enumerate(self.bscCurveFaranList):

			#Added small number to avoid taking logarithm of zero
			psi = 10*numpy.log(BSCs) - 10*numpy.log(BSCt )
			psiHat = psi.mean()
			mmse[count] = ((psi - psiHat)**2).mean()

		return mmse
