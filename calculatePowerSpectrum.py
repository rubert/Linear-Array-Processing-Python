from rfData import rfClass

class powerSpectrum(rfClass):

	def calculatePowerSpectrum(self, show = True):
		'''Give a power spectrum for the RF data using Welch's method'''


		from numpy import fft,zeros, arange
		from scipy import signal
		from matplotlib import pyplot
		
		self.readFrame()
		#going to use half the window size for Welch method
		windowY = max(rfClass.roiY) - min(rfClass.roiY)
		welchSub = windowY/2
		welchStep = welchSub/2

		tempLow = min(rfClass.roiY)
		tempAverage = zeros( welchSub )
		win = signal.hamming(welchSub)
		#loop through welch sub-windows
		for sw in range(4): #0,1,2,3
			temp = self.data[tempLow:tempLow + welchSub, min(rfClass.roiX):max(rfClass.roiX)]*win.reshape((welchSub,1))
			temp = temp - temp.mean(axis = 0)
			fourierData =fft.fft(temp, axis = 0)
			tempAverage += abs(fourierData).mean(axis = 1)
			tempLow = tempLow + welchStep
		
		#Work out sampling for DFT
		T = welchSub/self.fs
		deltaF = 1/T
		freqRange = deltaF*welchSub
		self.spectrumFrequencies = arange(0, freqRange, deltaF)/1E6
		self.spectrum = tempAverage
		if show:
			#show power spectrum
			fig = pyplot.figure() 
			ax = fig.add_subplot(1,1,1)
			ax.plot(self.spectrumFrequencies[0:welchSub/2], self.spectrum[0:welchSub/2])
			ax.set_title('Spectrum magnitude versus Frequency (MHz)')
			pyplot.show()


	def writeSpectrumToFile(self, fname, title = ' '):
		'''Write the power spectrum and the analysis region to file'''

		from numpy import fft,zeros, arange
		from scipy import signal
		from matplotlib import pyplot
		
		self.readFrame()
		#going to use half the window size for Welch method
		windowY = max(rfClass.roiY) - min(rfClass.roiY)
		welchSub = windowY/2
		welchStep = welchSub/2

		tempLow = min(rfClass.roiY)
		tempAverage = zeros( welchSub )
		win = signal.hamming(welchSub)
		#loop through welch sub-windows
		for sw in range(4): #0,1,2,3
			fourierData =fft.fft(self.data[tempLow:tempLow + welchSub, min(rfClass.roiX):max(rfClass.roiX)]*win.reshape((welchSub,1)), axis = 0)
			tempAverage += abs(fourierData).mean(axis = 1)
			tempLow = tempLow + welchStep
		
		#Work out sampling for DFT
		T = welchSub/self.fs
		deltaF = 1/T
		freqRange = deltaF*welchSub
		#show power spectrum
		from matplotlib.font_manager import FontProperties

		fig = pyplot.figure() 
		fig.canvas.set_window_title(title)
		ax = fig.add_subplot(1,2,1)
		f = arange(0, freqRange, deltaF)/1E6
		ax.plot(f[0:welchSub/2], tempAverage[0:welchSub/2])
		ax.set_title('Spectrum vs Frequency (MHz)')
		bMode = self.createRoiArray()	
		#import matplotlib and create plot
		import matplotlib.pyplot as plt
		import matplotlib.cm as cm
		palette = cm.gray
		palette.set_bad('r')

		ax2 = fig.add_subplot(1,2,2)
		ax2.set_title('Analysis Region')
		ax2.imshow(bMode, cmap = palette, extent = [0, self.fovX, self.fovY, 0])
		#Add title above all else	
		t = fig.text(0.5,
		0.95, title,
		horizontalalignment='center',
		fontproperties=FontProperties(size=16))
		pyplot.savefig(fname)	
		pyplot.show()	
	

	def calculateGeneralizedSpectrum(self, show = True):
		'''First pick a point from the image using ginput.  Then, compute the
		generalized spectrum around this region.  Use the generalized spectrum to compute
		the collapsed average.'''
		from numpy import zeros, arange, fft, exp, pi, outer, log10, sqrt
		from matplotlib import pyplot

		self.yLow = min(rfClass.roiY)
		self.yHigh = max(rfClass.roiY)
		self.xLow = min(rfClass.roiX)
		self.xHigh = max(rfClass.roiX)

		self.readFrame()
		#cut window size in half
		points = (self.yHigh - self.yLow)/2
		print "Points used is equal to: " +  str(points)
		#Work out the delay between the first point and the maximum intensity point
		dataWindow = self.data[self.yLow:self.yLow + points, self.xLow:self.xHigh]
		maxPointDelay = dataWindow.argmax(axis = 0 )/self.fs
		fourierData =fft.fft(dataWindow, axis = 0)
		tempPoints = dataWindow.shape[0]
		tempLines = dataWindow.shape[1]
		print "The number of lines used is equal to: " + str(tempLines)
		
		#Work out frequency spacing for DFT points
		deltaF = self.fs/dataWindow.shape[0]
		freq = arange(0, self.fs, deltaF)
		delta = outer(freq,maxPointDelay)
		phase = exp(1j*2*pi*delta)	
		print "Computing normalized GS"
		GS = zeros( (tempPoints, tempPoints) ) + 1j*zeros( (tempPoints, tempPoints) ) 
		for f1 in range(len(freq)):
			for f2 in range(len(freq)):
				for l in range(tempLines):
					outerProd = fourierData[f1,l]*phase[f1,l]*fourierData[f2,l].conjugate()*phase[f2,l].conjugate()
					if abs(outerProd) < 1E-8:
						outerProd = 0. + 1j*0.
					else:
						outerProd /= abs(outerProd)
					GS[f1,f2] += outerProd

		print "GS estimate computed over first time interval"
		#step ahead by 50% the window size and add in more contributions
		dataWindow = self.data[self.yLow + len(freq)/2: self.yLow + len(freq)/2 + tempPoints, self.xLow:self.xHigh]
		
		maxPointDelay = dataWindow.argmax(axis = 0 )/self.fs
		fourierData =fft.fft(dataWindow, axis = 0)
		tempPoints = dataWindow.shape[0]
		tempLines = dataWindow.shape[1]

		#Work out frequency spacing for DFT points
		deltaF = self.fs/dataWindow.shape[0]
		freq = arange(0, self.fs, deltaF)
		delta = outer(freq,maxPointDelay)
		phase = exp(1j*2*pi*delta)	
		
		for f1 in range(len(freq)):
			for f2 in range(len(freq)):
				for l in range(tempLines):
					outerProd = fourierData[f1,l]*phase[f1,l]*fourierData[f2,l].conjugate()*phase[f2,l].conjugate()
					if abs(outerProd) < 1E-8:
						outerProd = 0. + 1j*0.
					else:
						outerProd /= abs(outerProd)
					GS[f1,f2] += outerProd


		#step ahead by 50% the window size and add in more contributions
		dataWindow = self.data[self.yLow + len(freq): self.yLow + len(freq) + tempPoints, self.xLow:self.xHigh]
		
		maxPointDelay = dataWindow.argmax(axis = 0 )/self.fs
		fourierData =fft.fft(dataWindow, axis = 0)
		tempPoints = dataWindow.shape[0]
		tempLines = dataWindow.shape[1]

		#Work out frequency spacing for DFT points
		deltaF = self.fs/dataWindow.shape[0]
		freq = arange(0, self.fs, deltaF)
		delta = outer(freq,maxPointDelay)
		phase = exp(1j*2*pi*delta)	
		
		for f1 in range(len(freq)):
			for f2 in range(len(freq)):
				for l in range(tempLines):
					outerProd = fourierData[f1,l]*phase[f1,l]*fourierData[f2,l].conjugate()*phase[f2,l].conjugate()
					if abs(outerProd) < 1E-8:
						outerProd = 0. + 1j*0.
					else:
						outerProd /= abs(outerProd)
					GS[f1,f2] += outerProd


		######################################################
		########scale GS by its standard deviation############
		######################################################
		print "Computing standard deviation of GS estimates"
		#Work out the delay between the first point and the maximum intensity point
		dataWindow = self.data[self.yLow:self.yLow + points, self.xLow:self.xHigh]
		maxPointDelay = dataWindow.argmax(axis = 0 )/self.fs
		fourierData = fft.fft(dataWindow, axis = 0)
		tempPoints = dataWindow.shape[0]
		tempLines = dataWindow.shape[1]

		#Work out frequency spacing for DFT points
		deltaF = self.fs/dataWindow.shape[0]
		freq = arange(0, self.fs, deltaF)
		delta = outer(freq,maxPointDelay)
		phase = exp(1j*2*pi*delta)
		
		SL = zeros( (tempPoints, tempPoints) ) + 1j*zeros( (tempPoints, tempPoints) ) 
		for f1 in range(len(freq)):
			for f2 in range(len(freq)):
				for l in range(tempLines):
					outerProd = fourierData[f1,l]*phase[f1,l]*fourierData[f2,l].conjugate()*phase[f2,l].conjugate()
					if abs(outerProd) < 1E-8:
						outerProd = 0. + 1j*0.
					else:
						outerProd /= abs(outerProd)

					SL[f1,f2] += (GS[f1,f2] - outerProd)**2


		#step ahead by 50% the window size and add in more contributions
		dataWindow = self.data[self.yLow + len(freq)/2: self.yLow + len(freq)/2 + tempPoints, self.xLow:self.xHigh]
		
		maxPointDelay = dataWindow.argmax(axis = 0 )/self.fs
		fourierData =fft.fft(dataWindow, axis = 0)
		tempPoints = dataWindow.shape[0]
		tempLines = dataWindow.shape[1]

		#Work out frequency spacing for DFT points
		deltaF = self.fs/dataWindow.shape[0]
		freq = arange(0, self.fs, deltaF)
		delta = outer(freq,maxPointDelay)
		phase = exp(1j*2*pi*delta)
		
		for f1 in range(len(freq)):
			for f2 in range(len(freq)):
				for l in range(tempLines):
					outerProd = fourierData[f1,l]*phase[f1,l]*fourierData[f2,l].conjugate()*phase[f2,l].conjugate()
					if abs(outerProd) < 1E-8:
						outerProd = 0. + 1j*0.
					else:
						outerProd /= abs(outerProd)

					SL[f1,f2] += (GS[f1,f2] - outerProd)**2


		#step ahead by 50% the window size and add in more contributions
		dataWindow = self.data[self.yLow + len(freq): self.yLow + len(freq) + tempPoints, self.xLow:self.xHigh]
		
		maxPointDelay = dataWindow.argmax(axis = 0 )/self.fs
		fourierData =fft.fft(dataWindow, axis = 0)
		tempPoints = dataWindow.shape[0]
		tempLines = dataWindow.shape[1]

		#Work out frequency spacing for DFT points
		deltaF = self.fs/dataWindow.shape[0]
		freq = arange(0, self.fs, deltaF)
		delta = outer(freq,maxPointDelay)
		phase = exp(1j*2*pi*delta)
		
		for f1 in range(len(freq)):
			for f2 in range(len(freq)):
				for l in range(tempLines):
					outerProd = fourierData[f1,l]*phase[f1,l]*fourierData[f2,l].conjugate()*phase[f2,l].conjugate()
					if abs(outerProd) < 1E-8:
						outerProd = 0. + 1j*0.
					else:
						outerProd /= abs(outerProd)

					SL[f1,f2] += (GS[f1,f2] - outerProd)**2


		
		SL = sqrt( (1/(tempPoints*tempLines*3 - 1) )*SL )

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
		if show:
			fig2 = pyplot.figure() 
			ax2 = fig2.add_subplot(1,1,1)
			cax = ax2.imshow(gsDb[fIndexMin:fIndexMax, fIndexMin:fIndexMax], extent=[minFMHz, maxFMHz,  minFMHz, maxFMHz], interpolation = 'nearest', origin = 'lower')
			ax2.set_title('Generalized spectrum ')
			cbar = fig2.colorbar(cax)
			pyplot.show()


		########################################
		########COMPUTE COLLAPSED AVERAGE#######
		########################################
		print "Computing collapsed average"
		numFreq = len(range(fIndexMin, fIndexMax))	
		counts = zeros(numFreq)
		CA = zeros(numFreq) + 1j*zeros(numFreq)
		for f1 in range(fIndexMin,fIndexMax):
			for f2 in range(fIndexMin, fIndexMax):
				d = abs(f2-f1)		
				if d < numFreq:
					CA[d] = GS[f1,f2]
					counts[d] += 1


		self.CA = abs(CA)/counts
		self.freqDiff = arange(0, numFreq)*deltaF/1.E6
		if show:
			fig3 = pyplot.figure()
			ax3 = fig3.add_subplot(1,1,1)
			ax3.plot(self.freqDiff,self.CA)
			pyplot.show()
