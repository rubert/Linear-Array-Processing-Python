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
        self.windowY =int( windowYmm/self.deltaY)
        self.windowX =int( windowXmm/self.deltaX)

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
        self.psdPoints = 2*self.halfY
        self.psdPoints -= self.psdPoints%4
        self.gsPoints = self.psdPoints
        self.psdPoints /= 2
        freqStep = (freqHighMHz - freqLowMHz)/self.psdPoints

        self.psdFreq = numpy.arange(0,self.psdPoints)*freqStep  + freqLowMHz
        fracUnitCircle = (freqHighMHz - freqLowMHz)/(self.fs/10**6)
        self.cztW = numpy.exp(1j* (-2*numpy.pi*fracUnitCircle)/self.psdPoints )
        self.cztA = numpy.exp(1j* (2*numpy.pi*freqLowMHz/(self.fs/10**6) ) )


    def ComputeCollapsedAverageImage(self, vMax = None, itkFileName = None):

        import numpy
        self.ReadFrame()
        self.CaAvgImage =numpy.zeros( (self.numY, self.numX,7) )
        self.CaSlopeImage = numpy.zeros( (self.numY, self.numX,7) )
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
                caBinned,caSlope, sigma, mu = self.CalculateGeneralizedSpectrum(region = tempRegion)
                self.CaAvgImage[y,x, :] = caBinned
                self.CaSlopeImage[y,x] =caSlope
                self.centerFrequency[y,x] = sigma
                self.sigmaImage[y,x] = mu

                counter += 1
                print str(counter) + " completed of " + str(self.numY*self.numX)
        
        if vMax:
            self.caImage[self.caImage > vMax] = vMax
        
        #Write image to itk format
        if itkFileName:

            import itk
            p = itk.Point.F2()
            itkIm = itk.Image.F2.New()
            itkIm.SetRegions(self.sigmaImage.shape)
            itkIm.Allocate()
            itkIm.SetSpacing( [self.deltaY*self.caStepY, self.deltaX*self.caStepX] )
            itkIm.SetOrigin( [startY*self.deltaY, startX*self.deltaX] )
            
            writer = itk.ImageFileWriter.IF2.New()
           
            for freqBin in range(7):
                for countY in range(self.numY):
                    for countX in range(self.numX):
                        itkIm.SetPixel( [countY, countX], self.CaAvgImage[countY, countX, freqBin])

                writer.SetInput(itkIm)
                writer.SetFileName(itkFileName + 'caBin' + str(freqBin) + '.mhd')
                writer.Update()
            
            for freqBin in range(7):
                for countY in range(self.numY):
                    for countX in range(self.numX):
                        itkIm.SetPixel( [countY, countX], self.CaSlopeImage[countY, countX, freqBin])

                writer.SetInput(itkIm)
                writer.SetFileName(itkFileName + 'CaSlope' + str(freqBin) + '.mhd')
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

        import numpy

        self.SetRoiFixedSize(self.gsWindowXmm, self.gsWindowYmm)
        self.ReadFrame()
        dataBlock = self.data[self.roiY[0]:self.roiY[1], self.roiX[0]:self.roiX[1]]
        out1, out2, mu,sigma = self.CalculateGeneralizedSpectrum(dataBlock)

        from matplotlib import pyplot
        pyplot.plot(self.CAaxis, self.CA, 'k')
        pyplot.ylabel('Magnitude')
        pyplot.xlabel('Frequency difference (MHz)')
        pyplot.savefig(fname + 'collapsedAverage.pdf')
        pyplot.close()
   
        from matplotlib import pyplot
        import pdb
        pdb.set_trace()
        fig = pyplot.figure()
        axes = fig.add_subplot(111)
        spacingAxis = 1540./(2*self.CAaxis*10**6)*10**3
        axes.plot(self.CAaxis, self.CA, 'k')
        step = len(self.CAaxis)//5
        print "step is " + str(step)
        xticks = self.CAaxis[0::step]
        axes.set_xticks(xticks)
        xticklabelsFloats = spacingAxis[0::step]
        xtickLabels = [ '{:.2f}'.format(j) for j in xticklabelsFloats ]
        axes.set_xticklabels(xtickLabels )
        axes.set_ylabel('Magnitude')
        axes.set_xlabel('Scatterer Spacing (mm)')
        pyplot.savefig(fname + 'collapsedAveragePlotBySpacing.png')
        pyplot.close()
       
        #color different regions from by each sigma
        '''fig = pyplot.figure()
        axes = fig.add_subplot(111)
        spacingAxis = 1540./(2*self.CAaxis*10**6)*10**3
        axes.plot(self.CAaxis, self.CA, 'k')
        step = len(self.CAaxis)//5
        print "step is " + str(step)
        xticks = self.CAaxis[0::step]
        axes.set_xticks(xticks)
        xticklabelsFloats = spacingAxis[0::step]
        xtickLabels = [ '{:.2f}'.format(j) for j in xticklabelsFloats ]
        axes.set_xticklabels(xtickLabels )
        axes.set_ylabel('Magnitude')
        axes.set_xlabel('Scatterer Spacing (mm)')
        for i in range(3):
            axes.axvspan(sigma*(2*i), sigma*(2*i+1), alpha = .5) '''

        pyplot.savefig(fname + 'collapsedAveragePlotBySpacingColored.png')
        pyplot.close()
        
        pyplot.imshow( numpy.flipud(abs(self.GS)), extent= [self.gsFreq.min(), self.gsFreq.max(), self.gsFreq.min(),
        self.gsFreq.max()] )        
        pyplot.xlabel('Frequency (MHz)')
        pyplot.ylabel('Frequency (MHz)')
        pyplot.colorbar()
        pyplot.savefig(fname + 'generalizedSpectrum.png')
        pyplot.close()

    def CalculateGeneralizedSpectrum(self, region):
        '''Use 3 50% overlapping windows axially to compute PSD.
        Calculate Power spectrum first.  Normalize it to have a max value of 1.
        Fit this to a Gaussian.  f(x) = exp[ -(x - mu)**2 / 2 (sigma**2)]
        A Gaussian falls to -6 dB (1/4 of max value) at a little more than one sigma.
        
        Use full window length to compute axial signal.'''

        from scipy.signal import hann,convolve
        from scipy import optimize, interpolate
        import numpy
        from chirpz import chirpz

        maxDataWindow = region[0:self.gsPoints, :]
        
        windowFuncPsd = hann(self.psdPoints+2)[1:-1].reshape(self.psdPoints,1)
        powerSpectrum = numpy.zeros( self.psdPoints )

        #first compute the power spectrum, normalize it to maximum value
        for r in range(3):
            dataWindow = maxDataWindow[self.psdPoints/2*r:self.psdPoints/2*r + self.psdPoints, :]*windowFuncPsd
            fourierData = numpy.zeros( dataWindow.shape )
            for l in range(maxDataWindow.shape[1]):
                fourierData[:,l] = abs(chirpz(dataWindow[:,l], self.cztA, self.cztW, self.psdPoints))
            powerSpectrum += fourierData.mean(axis = 1)

        powerSpectrum /= powerSpectrum.max()

        #now fit the spectrum to a gaussian
        errfunc = lambda param,x,y: y - numpy.exp(- (x - param[0])**2/(2*param[1]**2) )
        param0 = (5., 1.0)
        args = (self.psdFreq,powerSpectrum)
        param, message = optimize.leastsq(errfunc, param0, args)
        mu = param[0]
        sigma = param[1]
        if sigma < .75:
            sigma = .75

        #set low and high frequency cutoffs based on output mu, sigma
        #3 Sigmas is 1% intensity on PSD is -20 dB
        lowCut = mu - 3*sigma
        if lowCut < 0:
            lowCut = 0
        highCut = mu + 3*sigma
        if highCut > self.fs/10**6:
            highCut = self.fs/10**6

        freqStep = (highCut - lowCut)/self.psdPoints
        spectrumFreq = numpy.arange(0,self.psdPoints)*freqStep  + lowCut
        fracUnitCircle = (highCut - lowCut)/(self.fs/10**6)
        cztW = numpy.exp(1j* (-2*numpy.pi*fracUnitCircle)/self.psdPoints )
        cztA = numpy.exp(1j* (2*numpy.pi*lowCut/(self.fs/10**6) ) )

        self.gsFreq = spectrumFreq

        ###########################
        ###Time averaged GS########
        ###########################
        GS = numpy.zeros( (self.psdPoints, self.psdPoints,3 ) , numpy.complex128)
        
        for r in range(3):
            dataWindow = maxDataWindow[self.psdPoints/2*r:self.psdPoints/2*r + self.psdPoints, :]
            maxPointDelay = dataWindow.argmax(axis = 0 )/self.fs
            dataWindow*=windowFuncPsd
            tempPoints = dataWindow.shape[0]
            tempLines = dataWindow.shape[1]
            fourierData = numpy.zeros( dataWindow.shape, numpy.complex128 )
            
            for l in range(tempLines):
                fourierData[:,l] = chirpz(dataWindow[:,l], cztA, cztW, self.psdPoints)

            #get point delay in seconds
            delta = numpy.outer(spectrumFreq*10**6,maxPointDelay)
            phase = numpy.exp(1j*2*numpy.pi*delta)
            
            for l in range(tempLines):
                outerProd = numpy.outer(fourierData[:,l]*phase[:,l], fourierData[:,l].conjugate()*phase[:,l].conjugate() )
                GS[:,:,r] += outerProd/abs(outerProd)

        numSegs = tempLines*3
        
        GS = GS.sum(axis = 2)/numSegs
        ##############################################
        ########COMPUTE COLLAPSED AVERAGE#############
        ##############################################
        #Along diagonal lines the scatterer spacing/frequency difference is the same
        #exclude all the entries that have a larger scatter spacing/lower frequency difference
        #than 1.5 mm:  .51 MHz 
        #Exclude all sizes less than .5 mm, which is 1.54 MHz
        numFreq = len(spectrumFreq)
        self.CA = numpy.zeros(numFreq)
        counts = numpy.zeros(numFreq)
       
        for f1 in range(numFreq):
            for f2 in range(numFreq):
                d = abs(f1-f2)
                self.CA[d] += abs(GS[f1,f2])
                counts[d] += 1

        
        self.CA /= counts
        self.CAaxis = numpy.arange(numFreq)*freqStep
        self.GS = GS

        #######################################
        ########COMPUTE THE AVERAGE OF THE#####
        ########COLLAPSED AVERAGE##############
        ###########AND THE SLOPE###############
        ###CA avg. bins:
        ###0.0-1.0 sigma
        ###1.0-2.0 sigma
        ###2.0-3.0 sigma
        ###3.0-4.0 sigma
        ###4.0-5.0 sigma
        ###5.0-6.0 sigma
        ###0.0-6.0 sigma
        #######################################

        ###Work out leakage location from PSD#######
        self.psdLimit = (1540./(2*self.gsWindowYmm*10**-3/2 ))/10**6
        secondDeriv = self.CA[0:-2] - 2*self.CA[1:-1] + self.CA[2:]
        psdInd = int(self.psdLimit/freqStep) 
        psdInd = secondDeriv[0:psdInd+10].argmax() + 1

        print "The PSD leakage stops at: " + str(secondDeriv[0:psdInd+10].argmax()  +1 )
        
        CAbinned = numpy.zeros(7)
        CAcounts = numpy.zeros(7)
        CAslope = numpy.zeros(7)
        
        for f, val in enumerate(self.CA):
            if f*freqStep > self.psdLimit and f*freqStep < 1.0*sigma:
                CAbinned[0] += val
                CAcounts[0] += 1
            elif f*freqStep > 1.0*sigma and f*freqStep < 2.0*sigma:
                CAbinned[1] += val
                CAcounts[1] += 1
            elif f*freqStep > 2.0*sigma and f*freqStep < 3.0*sigma:
                CAbinned[2] += val
                CAcounts[2] += 1
            elif f*freqStep > 3.0*sigma and f*freqStep < 4.0*sigma:
                CAbinned[3] += val
                CAcounts[3] += 1
            elif f*freqStep > 4.0*sigma and f*freqStep < 5.0*sigma:
                CAbinned[4] += val
                CAcounts[4] += 1
            elif f*freqStep > 5.0*sigma and f*freqStep < 6.0*sigma:
                CAbinned[5] += val
                CAcounts[5] += 1
            
            if f*freqStep > self.psdLimit:
                CAbinned[6] += val
                CAcounts[6] += 1


        for ind,count in enumerate(CAcounts):
            if val > 0:
                CAbinned[ind] /= count 

        
        ##compute slope of Collapsed Average
        for ind in range(7):
            if ind == 0 or ind == 6:
                lowCut = int(self.psdLimit/freqStep)
            else:
                lowCut = int(sigma*ind/freqStep)
            
            highCut = int(sigma*(ind+1)/freqStep)

            CAslope[ind] = self.ComputeSlope( self.CA[lowCut:highCut], self.CAaxis[lowCut:highCut])
        
        return CAbinned, CAslope, mu, sigma


    def ComputeSlope(self, yData, xData):
        from numpy import zeros, ones, array
        from numpy.linalg import lstsq
        #want to solve equation y = mx + b for m
        #[x1   1      [m     = [y1
        # x2   1       b ]      y2
        # x3   1]               y3]
        #
        #strain image will be smaller than displacement image
        A = ones( (len(xData),2) )
        A[:,0] = xData
        

        b = yData
        out = lstsq(A, b)
        xVec = out[0]
        return xVec[0]
