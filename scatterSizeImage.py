from rfData import rfClass
from attenuationBsc import attenuation
from faranScattering import faranBsc

class scattererSizeClass(attenuation):
			
	def ComputeBscImage(self):
		'''First compute the attenuation image, then use the attenuation image and the
		reference phantom spectrum to get the spectrum at each point.  Minimize difference
		between this power spectrum with a spectrum calculated from a Gaussian autocorrelation function.'''

		self.CalculateAttenuationImage(convertToRgb = False)
		self.InitializeTheoreticalBackscatter()
		numY = len(self.winCenterY)
		numX = len(self.winCenterX)
		self.scatSizeImage = numpy.zeros( (numY, numX) )	
		startY = self.winCenterY[self.lsqFitPoints//2] 
		startX = self.winCenterX[0]
		stepY = self.winCenterY[1] - self.winCenterY[0]
		stepX = self.winCenterX[1] - self.winCenterX[0]
		
		for y in self.winCenterY:
			attenuationDifference = self.deltaY*10**2*y*self.attenImage[0:y,:].mean()	
			for x in self.winCenterX:
				tempRegionSample = self.data[y - self.halfY:y + self.halfY + 1, x - self.halfX:x + self.halfX + 1]
				tempRegionRef = self.refRf.data[y - self.halfY:y + self.halfY + 1, x - self.halfX:x + self.halfX + 1]
				self.scatSizeImage[y, x] = self.CalculateScatSizeOnePoint(tempRegionSample, tempRegionRef)
			
		

	def InitializeTheoreticalBackscatter(self):
		'''Within the transducer bandwidth calculate the backscatter coefficients over a set of scatter size.
		Each array of backscatter coefficients is placed in a list whose different entries correpsond to different
		scatterer sizes.'''
		#for scatter sizes between 10 and 100 micrometers
		import numpy
		instance = faranBsc()
		self.bscFrequencies =numpy.arange(self.spectrumLowCutoff, self.spectrumHighCutoff)*self.spectrumFreqStep
		self.bscFaranSizes = numpy.arange(10,100,1)
		self.bscCurveFaranList =[]
		for d in self.bscFaranSizes:
			tempFreq, tempBsc = instance.calculateBSC(self, 1, .1, 12)
			self.bscCurveFaranList.append( tempBsc)


	def CalculateScatSizeOnePoint(self, sampleRegion, refRegion, attenuationDifference):
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
		diffs, size = self.ComputeBscCoefficients(rpmSpectrum[self.spectrumLowCutoff:self.spectrumHighCutoff] )		
		
		return size[diffs.argmin()]	
		
	def ComputeBscCoefficients(self, BSCs):
		'''Compute a theoretical backscatter curve over the transducer bandwidth.  Find the 
		logarithm of the difference between the two curves. 
		 
		
		Input:  The frequencies at which the spectrum to be compared with theory is calculated,
			probably should not be outside of transducer bandwidth
			
		Output:  Error and corresponding scatterer sizes'''
		import numpy
		for count,s in enumerate(self.bscFaranSizes):

			BSCt = bscCurveFaranList[count]
			psi = 10*numpy.log(BSCs) - 10*numpy.log(BSCt)
			psiHat = psi.mean()
			mmse[count] = ((psi - psiHat)**2).mean()

		return mmse
