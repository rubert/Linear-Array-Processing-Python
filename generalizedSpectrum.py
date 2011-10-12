from rfData import rfClass

class collapsedAverageImage(rfClass):

	def __init__(self, fname, ftype, freqLowMHz = 2.0, freqHighMHz = 12.0, windowYmm = 8, windowXmm = 8, overlapY = .75, overlapX = .75):

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
		
		self.caImage = numpy.zeros( (self.numY, self.numX) )
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
		self.spectrumFreqLowRes	= numpy.arange(self.bartlettY)*freqStepLowRes

		self.spectrumFreq = numpy.arange(0,self.bartlettY)*freqStep  + freqLowMHz	
		fracUnitCircle = (freqHighMHz - freqLowMHz)/(self.fs/10**6)
		self.cztW = numpy.exp(1j* (-2*numpy.pi*fracUnitCircle)/self.bartlettY ) 
		self.cztA = numpy.exp(1j* (2*numpy.pi*freqLowMHz/(self.fs/10**6) ) )
	
	
	def ComputeCollapsedAverageImage(self, vMax = None, itkFileName = None):
		'''Using 4 mm by 4 mm windows, create a collapsed average image.'''

		import numpy
		self.ReadFrame()

		
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
			
		#Compute whole image	
		for y in range(self.numY):
			for x in range(self.numX):
				tempRegion = self.data[self.winCenterY[y] - self.halfY:self.winCenterY[y] + self.halfY + 1, self.winCenterX[x] - self.halfX:self.winCenterX[x] + self.halfX + 1]
				self.CalculateGeneralizedSpectrum(region = tempRegion)
				self.caImage[y,x] = self.areaUnderCA
		
		if vMax:
			self.caImage[self.caImage > vMax] = vMax	
		self.caImageRGB = self.CreateParametricImage(self.caImage,[startY, startX], [self.caStepY, self.caStepX] )	

		#Write image to itk format
		if itkFileName:
			if 'mhd' not in itkFileName:
				itkFilename += '.mhd'
		
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
			writer.SetFileName(itkFileName)
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
		'''Accept a block of data and compute the collapsed average in the region.
		Average over 3 sub-blocks to compute GS'''
		
		from scipy.signal import hann,convolve
		import numpy
		from chirpz import chirpz
		points = region.shape[0]
		points -= points%4
		points /= 2
		
		maxDataWindow = region[0:2*points, :]

		#compute 3 fourier transforms and average them	
		#Cutting off the zero-value end points of the hann window
		#so it matches Matlab's definition of the function
		windowFunc = hann(points+2)[1:-1].reshape(points,1)
	
		##############
		###BLOCK 1####
		##############
		#Work out the delay between the first point and the maximum intensity point
		dataWindow = maxDataWindow[0:points, :]*windowFunc	
		tempPoints = dataWindow.shape[0]
		tempLines = dataWindow.shape[1]
		
		fourierData = numpy.zeros( dataWindow.shape ) + 1j*numpy.zeros( dataWindow.shape)
		for l in range(tempLines):	
			fourierData[:,l] = chirpz(dataWindow[:,l], self.cztA, self.cztW, points)
		
	
		#get point delay in seconds	
		maxPointDelay = dataWindow.argmax(axis = 0 )/self.fs
		#Work out frequency spacing for DFT points
		delta = numpy.outer(self.spectrumFreq*10**6,maxPointDelay)
		phase = numpy.exp(1j*2*numpy.pi*delta)	
		gsList1 = []
		
		for l in range(tempLines):
			outerProd = numpy.outer(fourierData[:,l]*phase[:,l], fourierData[:,l].conjugate()*phase[:,l].conjugate() )
			outerProd /= abs(outerProd)
			outerProd[outerProd == numpy.nan] = 0. + 1j*0.
			gsList1.append(outerProd.copy() )

		############
		###BLOCK 2##
		############
		#step ahead by 50% the window size and add in more contributions
		dataWindow = maxDataWindow[points/2: points/2 + points, :]*windowFunc
		tempPoints = dataWindow.shape[0]
		tempLines = dataWindow.shape[1]
		fourierData = numpy.zeros( dataWindow.shape ) + 1j*numpy.zeros( dataWindow.shape)
		for l in range(tempLines):	
			fourierData[:,l] = chirpz(dataWindow[:,l], self.cztA, self.cztW, points)
	
		maxPointDelay = dataWindow.argmax(axis = 0 )/self.fs
		

		#Work out frequency spacing for DFT points
		deltaF = self.spectrumFreq[1] - self.spectrumFreq[0]
		delta = numpy.outer(self.spectrumFreq*10**6,maxPointDelay)
		phase = numpy.exp(1j*2*numpy.pi*delta)	
	
		gsList2 = []	
		for l in range(tempLines):
			outerProd = numpy.outer(fourierData[:,l]*phase[:,l], fourierData[:,l].conjugate()*phase[:,l].conjugate() )
			outerProd /= abs(outerProd)
			outerProd[outerProd == numpy.nan] = 0. + 1j*0.
			gsList2.append(outerProd)


		############
		###BLOCK 3##
		############
		#step ahead by 50% the window size and add in more contributions
		dataWindow = maxDataWindow[points:2*points, :]*windowFunc		
		tempPoints = dataWindow.shape[0]
		tempLines = dataWindow.shape[1]
		fourierData = numpy.zeros( dataWindow.shape ) + 1j*numpy.zeros( dataWindow.shape)
		for l in range(tempLines):	
			fourierData[:,l] = chirpz(dataWindow[:,l], self.cztA, self.cztW, points)

		maxPointDelay = dataWindow.argmax(axis = 0 )/self.fs
		#Work out frequency spacing for DFT points
		deltaF = self.spectrumFreq[1] - self.spectrumFreq[0]
		delta = numpy.outer(self.spectrumFreq*10**6,maxPointDelay)
		phase = numpy.exp(1j*2*numpy.pi*delta)	
	
		gsList3 = []	
		for l in range(tempLines):
			outerProd = numpy.outer(fourierData[:,l]*phase[:,l], fourierData[:,l].conjugate()*phase[:,l].conjugate() )
			outerProd /= abs(outerProd)
			outerProd[outerProd == numpy.nan] = 0. + 1j*0.
			gsList3.append(outerProd)

		GS = numpy.zeros( (points, points) ) + 1j*numpy.zeros( (points, points) )
		for l in range(tempLines):
			GS += gsList1[l] + gsList2[l] + gsList3[l]

		
		######################################################
		########scale GS by its standard deviation############
		######################################################
		SL = numpy.zeros( (points, points) ) + 1j*numpy.zeros( (points, points) ) 
	
		for l in range(tempLines):
			SL += (GS - gsList1[l])**2 + (GS - gsList2[l])**2 + (GS - gsList3[l])**2  
		
		SL = numpy.sqrt( (1/(tempPoints*tempLines*3 - 1) )*SL )

		for f1 in range(len(self.spectrumFreq)):
			for f2 in range(len(self.spectrumFreq)):	
				if abs(SL[f1,f2]) > 1.E-8:
					GS[f1,f2] /= abs(SL[f1,f2])	
	


		########################################
		########COMPUTE COLLAPSED AVERAGE#######
		########################################
		numFreq = len(self.spectrumFreq)	
		counts = numpy.zeros(numFreq)
		CA = numpy.zeros(numFreq) + 1j*numpy.zeros(numFreq)
		for f1 in range(numFreq):
			for f2 in range(numFreq):
				d = abs(f2-f1)		
				if d < numFreq:
					CA[d] += GS[f1,f2]
					counts[d] += 1

		self.CAaxis = numpy.arange(numFreq)*(self.spectrumFreq[1] - self.spectrumFreq[0])
		self.CA = abs(CA)/counts

		#compute area under the collapsed average curve
		self.areaUnderCA = 0.

		for point in range(len(self.CA) - 1):
			self.areaUnderCA += (self.CA[point] + self.CA[point + 1])/2
