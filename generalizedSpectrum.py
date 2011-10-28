from rfData import rfClass

class collapsedAverageImage(rfClass):

    def __init__(self, fname, ftype, freqLowMHz = 0.0, freqHighMHz = 15.0, windowYmm = 8, windowXmm = 8, overlapY = .75, overlapX = .75):
        '''At each block this function will compute a power spectrum.  It will then fit that power spectrum to a Gaussian, and
        set the bounds on the chirpz transform to run from mu -sigma to mu + sigma.
        Input parameters:
        freqLowMHz:  A low cut-off for the chirpz transform of the power spectrum calculation.
        freqHighMHz:  The high cut-off for the chirpz transform of the power spectrum calcuation.
        windowYmm:  Block size in axial direction in mm
        windowXmm:  Block size in lateral direction in mm
        overlapY:  Axial overlap between blocks as a percentage.
        overlapX:  Lateral overlap between blocks as a percentage.
        '''


        import numpy
        #Set up RF data object
        super(collapsedAverageImage, self).__init__(fname, ftype)

        #figure out how many 6 mm by 4 mm windows fit into image
        self.gsWindowYmm = windowYmm
        self.gsWindowXmm = windowXmm
        self.windowX =int( windowYmm/self.deltaX)
        self.windowY =int( windowXmm/self.deltaY)

        #make the windows odd numbers
        if not self.windowY%2:
            self.windowY +=1

        if not self.windowX%2:
            self.windowX +=1

        self.halfY = (self.windowY -1)/2
        self.halfX = (self.windowX -1)/2

        #overlap the windows axially and laterally
        stepY = int(  (1-overlapY)*self.windowY )
        startY = self.halfY
        stopY = self.points - self.halfY - 1
        self.winCenterY = range(startY, stopY, stepY)
        self.numY = len(self.winCenterY)

        stepX = int( (1-overlapX)*self.windowX )
        startX = self.halfX
        stopX = self.lines - self.halfX - 1
        self.winCenterX = range(startX, stopX, stepX)
        self.numX = len(self.winCenterX)

        self.caStepX = stepX
        self.caStepY = stepY

        #Work out parameters for chirpz transform
        self.freqHigh = freqHighMHz
        self.freqLow = freqLowMHz

        #figure out block size in points
        self.bartlettY = 2*self.halfY
        self.bartlettY -= self.bartlettY%4
        self.bartlettY /= 2
        freqStep = (freqHighMHz - freqLowMHz)/self.bartlettY
        freqStepLowRes = self.fs/(10.**6*self.bartlettY)
        self.spectrumFreqLowRes = numpy.arange(self.bartlettY)*freqStepLowRes

        self.spectrumFreq = numpy.arange(0,self.bartlettY)*freqStep  + freqLowMHz
        fracUnitCircle = (freqHighMHz - freqLowMHz)/(self.fs/10**6)
        self.cztW = numpy.exp(1j* (-2*numpy.pi*fracUnitCircle)/self.bartlettY )
        self.cztA = numpy.exp(1j* (2*numpy.pi*freqLowMHz/(self.fs/10**6) ) )


    def ComputeCollapsedAverageImage(self, vMax = None, itkFileName = None):

        import numpy
        self.ReadFrame()
        self.caImage = numpy.zeros( (self.numY, self.numX) )
        self.centerFrequency= numpy.zeros( (self.numY, self.numX) )
        self.sigmaImage= numpy.zeros( (self.numY, self.numX) )

        startY = self.halfY
        startX = self.halfX

        #work out time to compute a point, then spit out time to calculate image
        from time import time
        y = x = 0
        t1 = time()
        tempRegion = self.data[self.winCenterY[y] - self.halfY:self.winCenterY[y] + self.halfY + 1, self.winCenterX[x] - self.halfX:self.winCenterX[x] + self.halfX + 1]
        self.CalculateGeneralizedSpectrum(region = tempRegion)
        t2 = time()
        print "Elapsed time was: " + str(t2-t1) + "seconds"
        print "Estimate time to compute an entire image is: "  + str( (t2-t1)*self.numY*self.numX/3600. ) + " hours"
        counter = 0
        #Compute whole image
        for y in range(self.numY):
            for x in range(self.numX):
                tempRegion = self.data[self.winCenterY[y] - self.halfY:self.winCenterY[y] + self.halfY + 1, self.winCenterX[x] - self.halfX:self.winCenterX[x] + self.halfX + 1]
                areaUnderCurve, sigma, mu = self.CalculateGeneralizedSpectrum(region = tempRegion)
                self.caImage[y,x] =areaUnderCurve
                self.centerFrequency[y,x] = sigma
                self.sigmaImage[y,x] = mu

                counter += 1
                print str(counter) + " completed of " + str(self.numY*self.numX)
        
        if vMax:
            self.caImage[self.caImage > vMax] = vMax

        self.caImageRGB = self.CreateParametricImage(self.caImage,[startY, startX], [self.caStepY, self.caStepX] )

        #Write image to itk format
        if itkFileName:

            import itk
            itkIm = itk.Image.F2.New()
            itkIm.SetRegions(self.caImage.shape)
            itkIm.Allocate()
            for countY in range(self.numY):
                for countX in range(self.numX):
                    itkIm.SetPixel( [countY, countX], self.caImage[countY, countX])

            itkIm.SetSpacing( [self.deltaY*self.caStepY, self.deltaX*self.caStepX] )
            itkIm.SetOrigin( [startY*self.deltaY, startX*self.deltaX] )
            writer = itk.ImageFileWriter.IF2.New()
            writer.SetInput(itkIm)
            writer.SetFileName(itkFileName + '.mhd')
            writer.Update()


            for countY in range(self.numY):
                for countX in range(self.numX):
                    itkIm.SetPixel( [countY, countX], self.centerFrequency[countY, countX])

            writer.SetInput(itkIm)
            writer.SetFileName(itkFileName + 'centerFreq.mhd')
            writer.Update()

            for countY in range(self.numY):
                for countX in range(self.numX):
                    itkIm.SetPixel( [countY, countX], self.sigmaImage[countY, countX])

            writer.SetInput(itkIm)
            writer.SetFileName(itkFileName + 'sigma.mhd')
            writer.Update()

    def ComputeCollapsedAverageInRegion(self, fname):
        '''Create a windowY by windowX region, then compute the collapsed average in that
        region.  Save a collapsed average image.'''


        if 'png' not in fname:
            fname += '.png'
        self.SetRoiFixedSize(self.gsWindowXmm, self.gsWindowYmm)
        self.ReadFrame()
        dataBlock = self.data[self.roiY[0]:self.roiY[1], self.roiX[0]:self.roiX[1]]
        self.CalculateGeneralizedSpectrum(dataBlock)

        from matplotlib import pyplot
        pyplot.plot(self.CAaxis, self.CA)
        pyplot.ylabel('Magnitude')
        pyplot.xlabel('Frequency difference (MHz)')
        pyplot.savefig(fname)


    def CalculateGeneralizedSpectrum(self, region):
        '''Use 3 50% overlapping windows axially to compute Fourier transforms.
        Calculate Power spectrum first.  Normalize it to have a max value of 1.
        Fit this to a Gaussian.  f(x) = exp[ -(x - mu)**2 / 2 (sigma**2)]
        A Gaussian falls to -6 dB (1/4 of max value) at a little more than one sigma '''

        from scipy.signal import hann,convolve
        from scipy import optimize, interpolate
        import numpy
        from chirpz import chirpz

        points = region.shape[0]
        points -= points%4
        points /= 2

        maxDataWindow = region[0:2*points, :]
        windowFunc = hann(points+2)[1:-1].reshape(points,1)
        powerSpectrum = numpy.zeros( points )

        #first compute the power spectrum, normalize it to maximum value
        for r in range(3):
            dataWindow = maxDataWindow[points/2*r:points/2*r + points, :]*windowFunc
            fourierData = numpy.zeros( dataWindow.shape )
            for l in range(maxDataWindow.shape[1]):
                fourierData[:,l] = abs(chirpz(dataWindow[:,l], self.cztA, self.cztW, points))
            powerSpectrum += fourierData.mean(axis = 1)

        powerSpectrum /= powerSpectrum.max()

        #now fit the spectrum to a gaussian
        errfunc = lambda param,x,y: y - numpy.exp(- (x - param[0])**2/(2*param[1]**2) )
        param0 = (5., 1.0)
        args = (self.spectrumFreq,powerSpectrum)
        param, message = optimize.leastsq(errfunc, param0, args)
        mu = param[0]
        sigma = param[1]

        #set low and high frequency cutoffs based on output mu, sigma
        lowCut = mu - sigma
        highCut = mu + sigma

        freqStep = (highCut - lowCut)/self.bartlettY
        spectrumFreq = numpy.arange(0,self.bartlettY)*freqStep  + lowCut
        fracUnitCircle = (highCut - lowCut)/(self.fs/10**6)
        cztW = numpy.exp(1j* (-2*numpy.pi*fracUnitCircle)/self.bartlettY )
        cztA = numpy.exp(1j* (2*numpy.pi*lowCut/(self.fs/10**6) ) )


        #compute 3 fourier transforms and average them
        #Cutting off the zero-value end points of the hann window
        #so it matches Matlab's definition of the function
        normGS = numpy.zeros( (points, points), numpy.complex128)

        ##############
        ###BLOCK 1####
        ##############
        #Work out the delay between the first point and the maximum intensity point
        dataWindow = maxDataWindow[0:points, :]*windowFunc
        tempPoints = dataWindow.shape[0]
        tempLines = dataWindow.shape[1]

        fourierData = numpy.zeros( dataWindow.shape ) + 1j*numpy.zeros( dataWindow.shape)
        for l in range(tempLines):
            fourierData[:,l] = chirpz(dataWindow[:,l], cztA, cztW, points)


        #get point delay in seconds
        maxPointDelay = dataWindow.argmax(axis = 0 )/self.fs
        #Work out frequency spacing for DFT points
        delta = numpy.outer(spectrumFreq*10**6,maxPointDelay)
        phase = numpy.exp(1j*2*numpy.pi*delta)
        gsList1 = []

        for l in range(tempLines):
            outerProd = numpy.outer(fourierData[:,l]*phase[:,l], fourierData[:,l].conjugate()*phase[:,l].conjugate() )
            #normalize by the diagonal elements
            PSD = abs(numpy.diag(outerProd))
            PSDproduct = numpy.sqrt( numpy.outer( PSD, PSD ) )
            normGS = outerProd/PSDproduct
            gsList1.append(normGS.copy() )

        ############
        ###BLOCK 2##
        ############
        #step ahead by 50% the window size and add in more contributions
        dataWindow = maxDataWindow[points/2: points/2 + points, :]*windowFunc
        tempPoints = dataWindow.shape[0]
        tempLines = dataWindow.shape[1]
        fourierData = numpy.zeros( dataWindow.shape ) + 1j*numpy.zeros( dataWindow.shape)
        for l in range(tempLines):
            fourierData[:,l] = chirpz(dataWindow[:,l], cztA, cztW, points)

        maxPointDelay = dataWindow.argmax(axis = 0 )/self.fs


        #Work out frequency spacing for DFT points
        deltaF = spectrumFreq[1] - spectrumFreq[0]
        delta = numpy.outer(spectrumFreq*10**6,maxPointDelay)
        phase = numpy.exp(1j*2*numpy.pi*delta)

        gsList2 = []
        for l in range(tempLines):
            outerProd = numpy.outer(fourierData[:,l]*phase[:,l], fourierData[:,l].conjugate()*phase[:,l].conjugate() )
            PSD = abs(numpy.diag(outerProd))
            PSDproduct = numpy.sqrt( numpy.outer( PSD, PSD ) )
            #normalize by the diagonal elements
            normGS = outerProd/PSDproduct
            gsList2.append(normGS.copy())


        ############
        ###BLOCK 3##
        ############
        #step ahead by 50% the window size and add in more contributions
        dataWindow = maxDataWindow[points:2*points, :]*windowFunc
        tempPoints = dataWindow.shape[0]
        tempLines = dataWindow.shape[1]
        fourierData = numpy.zeros( dataWindow.shape ) + 1j*numpy.zeros( dataWindow.shape)
        for l in range(tempLines):
            fourierData[:,l] = chirpz(dataWindow[:,l], cztA, cztW, points)

        maxPointDelay = dataWindow.argmax(axis = 0 )/self.fs
        #Work out frequency spacing for DFT points
        deltaF = spectrumFreq[1] - spectrumFreq[0]
        delta = numpy.outer(spectrumFreq*10**6,maxPointDelay)
        phase = numpy.exp(1j*2*numpy.pi*delta)

        gsList3 = []

        for l in range(tempLines):
            outerProd = numpy.outer(fourierData[:,l]*phase[:,l], fourierData[:,l].conjugate()*phase[:,l].conjugate() )
            PSD = abs(numpy.diag(outerProd))
            PSDproduct = numpy.sqrt( numpy.outer( PSD, PSD ) )
            #normalize by the diagonal elements
            normGS = outerProd/PSDproduct
            gsList3.append(normGS.copy())

        GS = numpy.zeros( (points, points) ) + 1j*numpy.zeros( (points, points) )
        for l in range(tempLines):
            GS += gsList1[l] + gsList2[l] + gsList3[l]



        ########################################
        ########COMPUTE ENERGY IN OFF DIAGONAL##
        #######DIVIDE BY AREA OCCUPIED##########
        ########################################
        numFreq = len(spectrumFreq)
        counts = numpy.zeros(numFreq)
        CA = 0.
        for f1 in range(numFreq):
            for f2 in range(numFreq):
                if abs(f1 - f2) > 0:
                    CA += abs(GS[f1,f2])

        CA = CA/( (numFreq**2 - numFreq)*deltaF)
        return CA, mu, sigma
