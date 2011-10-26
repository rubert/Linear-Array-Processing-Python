from rfData import rfClass
from attenuation import attenuation
from faranScattering import faranBsc

class scatterSizeClass(attenuation):

    def __init__(self, sampleName, refName, dataType, bscCoefficients, bscFreqArray, numRefFrames = 0, refAttenuation = .5, freqLow = 2., freqHigh = 8., attenuationKernelSizeYmm = 15, blockYmm = 8, blockXmm = 10, overlapY = .85, overlapX = .85, bscFitRadius = 1.0, centerFreqSimulation = 5.0, sigmaSimulation = 1.0 ):

        '''The additional information this code needs on top of attenuation is a set of theoretical
        backscatter coefficients to match against, and the backscatter coefficients of the
        reference phantom.'''


        super(scatterSizeClass, self).__init__(sampleName, refName, dataType,  numRefFrames, refAttenuation, freqLow, freqHigh, attenuationKernelSizeYmm, blockYmm, blockXmm, overlapY , overlapX , bscFitRadius, centerFreqSimulation, sigmaSimulation )
        #Get backscatter coefficients of the reference phantom in the same analysis range as the frequency Kernel
        import numpy
        from scipy.interpolate import interp1d
        bscNew = interp1d(bscFreqArray, bscCoefficients)
        self.bscReference = bscNew(self.spectrumFreq)

        self.bscFaranSizes = numpy.arange(1,150)


    def ComputeScattererSizeImage(self, vmin = None, vmax = None, itkFileName = None):
        '''First compute the attenuation image, then use the attenuation image and the
        reference phantom spectrum to get the spectrum at each point.  Minimize difference
        between this power spectrum with a spectrum calculated from a Gaussian autocorrelation function.
        Assume the attenuation outside the known region of the attenuation image is equal to the mean
        of the attenuation image.'''

        import numpy

        self.CalculateAttenuationImage()

        self.InitializeTheoreticalBackscatter()
        numY = len(self.blockCenterY)
        numX = len(self.blockCenterX)
        self.scatSizeImage = numpy.zeros( (numY, numX) )
        startY = self.blockCenterY[0]
        startX = self.blockCenterX[0]
        stepY = self.blockCenterY[1] - self.blockCenterY[0]
        stepX = self.blockCenterX[1] - self.blockCenterX[0]

        #assume sample and reference have same attenuation in out of bounds region.
        for xParam, xRf in enumerate(self.blockCenterX):
            for yParam,yRf in enumerate(self.blockCenterY):

                if yParam < self.halfLsq:
                    betadiff = 0
                else:
                    tmp = self.attenuationImage[0:yParam - self.halfLsq + 1,xParam].mean()
                    if tmp < 0.:
                        tmp = 0
                    betadiff = tmp - self.betaRef
                betadiff /= 8.686
                depthCm = (self.deltaY*yRf)/10
                rpmSpectrum = self.sampleSpectrum[:,yParam,xParam]/self.refSpectrum[:,yParam]*numpy.exp(4*betadiff*depthCm*self.spectrumFreq)*self.bscReference
                middleInd = self.refSpectrum[:,yParam].argmax()
                lowInd = middleInd - self.radiusInPoints
                if lowInd < 0:
                    lowInd = 0
                highInd = middleInd + self.radiusInPoints
                if highInd > self.refSpectrum.shape[0]:
                    highInd = self.refSpectrum.shape[0]
                diffs = self.CompareBscCoefficients(rpmSpectrum , lowInd, highInd)
                self.scatSizeImage[yParam, xParam] = self.bscFaranSizes[diffs.argmin()]

        print "Mean Scatterer size is: " + str(self.scatSizeImage.mean())
        #convert scatterer size image to RGB parametric image
        if vmin:
            self.scatSizeImage[self.scatSizeImage < vmin] = vmin
        if vmax:
            self.scatSizeImage[self.scatSizeImage > vmax] = vmax

        self.scatSizeImageRGB = self.CreateParametricImage(self.scatSizeImage,[startY, startX], [stepY, stepX] )

        #Write image to itk format
        if itkFileName:
            if 'mhd' not in itkFileName:
                itkFilename += '.mhd'

            import itk
            itkIm = itk.Image.F2.New()
            itkIm.SetRegions(self.scatSizeImage.shape)
            itkIm.Allocate()
            for countY in range(numY):
                for countX in range(numX):
                    itkIm.SetPixel( [countY, countX], self.scatSizeImage[countY, countX])

            itkIm.SetSpacing( [self.deltaY*stepY, self.deltaX*stepY] )
            itkIm.SetOrigin( [startY*self.deltaY, startX*self.deltaX] )

            writer = itk.ImageFileWriter[itk.Image.F2]
            writer.SetInput(itkIm)
            writer.SetFileName(itkFileName)
            writer.Update()

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
        self.bscCurveFaranList =[]
        for d in self.bscFaranSizes:
            tempBsc = faranInstance.calculateBSC(self.spectrumFreq, d)
            self.bscCurveFaranList.append( tempBsc.copy())


    def CompareBscCoefficients(self, BSCs, startIndex, stopIndex):
        '''Compute a theoretical backscatter curve over the transducer bandwidth.  Find the
        logarithm of the difference between the two curves.


        Input:  The frequencies at which the spectrum to be compared with theory is calculated,
                probably should not be outside of transducer bandwidth

        Output:  Error and corresponding scatterer sizes'''
        import numpy

        if startIndex < 0:
            startIndex = 0
        if stopIndex > len(BSCs):
            stopIndex = len(BSCs)

        #Work out start index and stop index
        mmse = numpy.zeros( len(self.bscCurveFaranList) )
        for count,BSCt in enumerate(self.bscCurveFaranList):

            #Added small number to avoid taking logarithm of zero
            psi = 10*numpy.log(BSCs[startIndex:stopIndex]) - 10*numpy.log(BSCt[startIndex:stopIndex])
            psiHat = psi.mean()
            mmse[count] = ((psi - psiHat)**2).mean()

        mmse[numpy.isnan(mmse)] = 1.e15

        return mmse
