from rfData import rfClass
from attenuationBsc import attenuation
from faranScattering import faranBsc

class scattererSizeClass(attenuation):
			
	def ComputeBscImage(self, convertToRgb = True, vmin = None, vmax = None):
		'''First compute the attenuation image, then use the attenuation image and the
		reference phantom spectrum to get the spectrum at each point.  Minimize difference
		between this power spectrum with a spectrum calculated from a Gaussian autocorrelation function.'''
		
		import numpy

		self.CalculateAttenuationImage(convertToRgb = False)
		self.InitializeTheoreticalBackscatter()
		numY = len(self.attenCenterY)
		numX = len(self.blockCenterX)
		self.scatSizeImage = numpy.zeros( (numY, numX) )	
		startY = self.attenCenterY[0] 
		startX = self.blockCenterX[0]
		stepY = self.attenCenterY[1] - self.attenCenterY[0]
		stepX = self.blockCenterX[1] - self.blockCenterX[0]
		
		for yParam,yRf in enumerate(self.attenCenterY):
			for xParam, xRf in enumerate(self.blockCenterX):
				betaDiff = self.attenuationImage[0:yParam + 1,xParam].mean()	
				depthCm = (self.deltaY*yRf)/10
				tempRegionSample = self.data[yRf - self.halfY:yRf + self.halfY + 1, xRf - self.halfX:xRf + self.halfX + 1]
				tempRegionRef = self.refRf.data[yRf - self.halfY:yRf + self.halfY + 1, xRf - self.halfX:xRf + self.halfX + 1]
				self.scatSizeImage[yParam, xParam] = self.CalculateScatSizeOnePoint(tempRegionSample, tempRegionRef, betaDiff, depthCm)
			
	
		#convert scatterer size image to RGB parametric image	
		if vmin:
			self.scatSizeImage[self.scatSizeImage < vmin] = vmin
		if vmax:
			self.scatSizeImage[self.scatSizeImage > vmax] = vmax
		if convertToRgb:
			self.scatSizeImage = self.CreateParametricImage(self.scatSizeImage,[startY, startX], [stepY, stepX] )

	def InitializeTheoreticalBackscatter(self):
		'''Within the transducer bandwidth calculate the backscatter coefficients over a set of scatter size.
		Each array of backscatter coefficients is placed in a list whose different entries correpsond to different
		scatterer sizes.'''
		#for scatter sizes between 10 and 100 micrometers
		import numpy
		faranInstance =  faranBsc()

		#I get a lot of divide by zeros when computing the theoretical backscatter
		#coefficients

		instance = faranBsc()
		self.bscFaranSizes = numpy.arange(10,200,1)
		self.bscCurveFaranList =[]
		for d in self.bscFaranSizes:
			tempBsc = faranInstance.calculateBSC(self.spectrumFreq, d)
			self.bscCurveFaranList.append( tempBsc.copy())


	def CalculateScatSizeOnePoint(self, sampleRegion, refRegion, betaDiff, depthCm):
		'''This function calculates power spectrum ratios and multiplies by an attenuation 
		difference with a known reference phantom.  The units on the attenuation difference are in 
		dB/(cm MHz)'''
		import numpy
		
		betaDiff /= 8.686
		
		fftSample = self.CalculateSpectrumBlock(sampleRegion)
		fftRef = self.CalculateSpectrumBlock(refRegion)				
		###Divide sample spectrum by reference spectrum to perform deconvolution
		rpmSpectrum = fftSample/fftRef*numpy.exp(-4*betaDiff*depthCm*self.spectrumFreq)
		from matplotlib import pyplot
		import pdb
		pdb.set_trace()
		diffs = self.ComputeBscCoefficients(rpmSpectrum )		
		
		return self.bscFaranSizes[diffs.argmin()]	
		
	def ComputeBscCoefficients(self, BSCs):
		'''Compute a theoretical backscatter curve over the transducer bandwidth.  Find the 
		logarithm of the difference between the two curves. 
		 
		
		Input:  The frequencies at which the spectrum to be compared with theory is calculated,
			probably should not be outside of transducer bandwidth
			
		Output:  Error and corresponding scatterer sizes'''
		import numpy
		
		oldWarningSettings = numpy.seterr(all = 'ignore')
		
		mmse = numpy.zeros( len(self.bscCurveFaranList) )
		for count,BSCt in enumerate(self.bscCurveFaranList):

			#Added small number to avoid taking logarithm of zero
			psi = 10*numpy.log(BSCs) - 10*numpy.log(BSCt)
			psiHat = psi.mean()
			mmse[count] = ((psi - psiHat)**2).mean()

		mmse[numpy.isnan(mmse)] = 1.e15
		numpy.seterr(**oldWarningSettings)
		return mmse
