from rfData import rfClass

class collapsedAverageImage(rfClass):

	def ComputeCollapsedAverageImage(self, windowYmm = 8, windowXmm = 8, overlapY = .75, overlapX = .75):
		'''Using 4 mm by 4 mm windows, create a collapsed average image.'''

		import numpy
		self.ReadFrame()
		#figure out how many 6 mm by 4 mm windows fit into image		
		windowX =int( windowYmm/self.deltaX)
		windowY =int( windowXmm/self.deltaY)
		
		#make the windows odd numbers
		if not windowY%2:
			windowY +=1

		if not windowX%2:
			windowX +=1

		halfY = (windowY -1)/2
		halfX = (windowX -1)/2
			
		#overlap the windows axially by 50%
		stepY = int(  (1-overlapY)*windowY )
		startY = halfY
		stopY = self.points - halfY - 1
		winCenterY = range(startY, stopY, stepY)
		numY = len(winCenterY)
	
		stepX = int( (1-overlapX)*windowX )	
		startX = halfX
		stopX = self.lines - halfX - 1
		winCenterX = range(startX, stopX, stepX)
		numX = len(winCenterX)
		
		self.caImage = numpy.zeros( (numY, numX) )
		self.caStepX = stepX
		self.caStepY = stepY
			
		#work out time to compute a point, then spit out time to calculate image
		from time import time
		y = x = 0
		t1 = time()
		tempRegion = self.data[winCenterY[y] - halfY:winCenterY[y] + halfY + 1, winCenterX[x] - halfX:winCenterX[x] + halfX + 1]
		self.CalculateGeneralizedSpectrum(region = tempRegion)
		t2 = time()
		print "Elapsed time was: " + str(t2-t1) + "seconds"
		print "Estimate time to compute an entire image is: "  + str( (t2-t1)*numY*numX/3600. ) + " hours"
				
		for y in range(numY):
			for x in range(numX):
				tempRegion = self.data[winCenterY[y] - halfY:winCenterY[y] + halfY + 1, winCenterX[x] - halfX:winCenterX[x] + halfX + 1]
				self.CalculateGeneralizedSpectrum(region = tempRegion)
				self.caImage[y,x] = self.areaUnderCA
		
			
		self.caImage = self.CreateParametricImage(self.caImage,[startY, startX], [self.caStepY, self.caStepX] )	
	
	
	def CalculateGeneralizedSpectrum(self, region):
		'''Accept a block of data and compute the collapsed average in the region.'''
		
		import numpy	
		from scipy.signal import hamming

		points = region.shape[0]
		points -= points%4
		maxDataWindow = region[0:points, :]
		points /= 2
		
		#Work out the delay between the first point and the maximum intensity point
		windowFunc = hamming(points).reshape(points,1)
		dataWindow = maxDataWindow[0:points, :]*windowFunc	
		maxPointDelay = dataWindow.argmax(axis = 0 )/self.fs
		fourierData =numpy.fft.fft(dataWindow, axis = 0)
		tempPoints = dataWindow.shape[0]
		tempLines = dataWindow.shape[1]
		
		#Work out frequency spacing for DFT points
		deltaF = self.fs/dataWindow.shape[0]
		freq =numpy.arange(0, self.fs, deltaF)
		delta = numpy.outer(freq,maxPointDelay)
		phase = numpy.exp(1j*2*numpy.pi*delta)	
		gsList1 = []
		
		for l in range(tempLines):
			outerProd = numpy.outer(fourierData[:,l]*phase[:,l], fourierData[:,l].conjugate()*phase[:,l].conjugate() )
			outerProd[abs(outerProd) < 1E-8] =  0. + 1j*0.
			outerProd /= abs(outerProd)
			outerProd[outerProd == numpy.nan] = 0. + 1j*0.
			gsList1.append(outerProd.copy() )

		#step ahead by 50% the window size and add in more contributions
		dataWindow = maxDataWindow[points/2: points/2 + points, :]*windowFunc
		
		maxPointDelay = dataWindow.argmax(axis = 0 )/self.fs
		fourierData =numpy.fft.fft(dataWindow, axis = 0)
		tempPoints = dataWindow.shape[0]
		tempLines = dataWindow.shape[1]

		#Work out frequency spacing for DFT points
		deltaF = self.fs/dataWindow.shape[0]
		freq = numpy.arange(0, self.fs, deltaF)
		delta = numpy.outer(freq,maxPointDelay)
		phase = numpy.exp(1j*2*numpy.pi*delta)	
	
		gsList2 = []	
		for l in range(tempLines):
			outerProd = numpy.outer(fourierData[:,l]*phase[:,l], fourierData[:,l].conjugate()*phase[:,l].conjugate() )
			outerProd[abs(outerProd) < 1E-8] =  0. + 1j*0.
			outerProd /= abs(outerProd)
			outerProd[outerProd == numpy.nan] = 0. + 1j*0.
			gsList2.append(outerProd)


		#step ahead by 50% the window size and add in more contributions
		dataWindow = maxDataWindow[points:2*points, :]*windowFunc		
		
		maxPointDelay = dataWindow.argmax(axis = 0 )/self.fs
		fourierData =numpy.fft.fft(dataWindow, axis = 0)
		tempPoints = dataWindow.shape[0]
		tempLines = dataWindow.shape[1]

		#Work out frequency spacing for DFT points
		deltaF = self.fs/dataWindow.shape[0]
		freq = numpy.arange(0, self.fs, deltaF)
		delta = numpy.outer(freq,maxPointDelay)
		phase = numpy.exp(1j*2*numpy.pi*delta)	
	
		gsList3 = []	
		for l in range(tempLines):
			outerProd = numpy.outer(fourierData[:,l]*phase[:,l], fourierData[:,l].conjugate()*phase[:,l].conjugate() )
			outerProd[abs(outerProd) < 1E-8] =  0. + 1j*0.
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

		for f1 in range(len(freq)):
			for f2 in range(len(freq)):	
				if abs(SL[f1,f2]) > 1.E-8:
					GS[f1,f2] /= abs(SL[f1,f2])	
	
		#only show min:max MHz
		minFMHz = 5.
		maxFMHz = 13.
		fIndexMax = int(maxFMHz*1.E6/deltaF)
		fIndexMin = int(minFMHz*1.E6/deltaF)
				
		gsDb = abs(GS)
		gsDb -= gsDb.max()


		########################################
		########COMPUTE COLLAPSED AVERAGE#######
		########################################
		numFreq = len(range(fIndexMin, fIndexMax))	
		counts = numpy.zeros(numFreq)
		CA = numpy.zeros(numFreq) + 1j*numpy.zeros(numFreq)
		for f1 in range(fIndexMin,fIndexMax):
			for f2 in range(fIndexMin, fIndexMax):
				d = abs(f2-f1)		
				if d < numFreq:
					CA[d] += GS[f1,f2]
					counts[d] += 1


		self.CA = abs(CA)/counts
		self.freqDiff = numpy.arange(0, numFreq)*deltaF/1.E6

		#normalize to value at origin
		self.CA = self.CA/self.CA[0]

		#compute area under the collapsed average curve
		self.areaUnderCA = 0.

		for point in range(len(self.CA) - 1):
			self.areaUnderCA += (self.CA[point] + self.CA[point + 1])/2*deltaF/1.E6

