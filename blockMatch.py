#first define the class unique priority queue
from Queue import PriorityQueue
import heapq
from rfData import rfClass
import cv
import numpy as np

class UniquePriorityQueue(PriorityQueue):
    def _init(self, maxsize):
        PriorityQueue._init(self, maxsize)
        self.values = set()

    def _put(self, item, heappush=heapq.heappush):
        if item[1] not in self.values:
            self.values.add(item[1])
            PriorityQueue._put(self, item, heappush)
#        else:  I should check to see if the quality is better, if so, remove old item
#            print 'dupe',item[1]

    def _get(self, heappop=heapq.heappop):
        item = PriorityQueue._get(self, heappop)
        self.values.remove(item[1])
        return item


class blockMatchClass(rfClass):


	def __init__(self, fname, dataType,  windowYmm = 1.0, windowXmm = 10., rangeYmm = 1.0, rangeXmm = .6, overlap = .65, strainKernelmm = 6.0):
	
		super(blockMatchClass, self).__init__(fname, dataType)	

		self.windowYmm = windowYmm	
		self.windowXmm = windowXmm
		self.rangeYmm =  rangeYmm
		self.rangeXmm = rangeXmm
	
		#axial block size	
		self.windowY = int(self.windowYmm/self.deltaY)
		if not self.windowY%2:
			self.windowY += 1
		self.halfY = self.windowY//2
	
		#lateral block size
		self.windowX = int(self.windowXmm/self.deltaX)
		if not self.windowX%2:
			self.windowX +=1
		self.halfX = self.windowX//2
		

		#search range
		self.rangeY = int(self.rangeYmm/self.deltaY)
		self.rangeX = int(self.rangeXmm/self.deltaX)
		self.smallRangeY = 3
		self.smallRangeX = 2

		#overlap and step
		self.overlap = overlap
		self.stepY =int( (1 - overlap)*self.windowY )
		####work out boundaries
		startY = 0 + self.rangeY + self.smallRangeY + self.halfY
		stepY =int( (1 - self.overlap)*self.windowY )
		stopY = self.points - self.halfY - self.smallRangeY - self.rangeY #last pixel that can fit a window

		startX = self.halfX + self.rangeX + self.smallRangeX
		stopX = self.lines - self.halfX - self.smallRangeX - self.rangeX

		self.startY = startY
		self.stopY =stopY
		self.startX = startX
		self.stopX = stopX
	###create arrays containing window centers in rf data coordinates
		self.windowCenterY = range(startY,stopY, stepY)
		self.numY = len(self.windowCenterY)

		self.windowCenterX = range(startX,stopX)
		self.numX = len(self.windowCenterX)

	
	#work out strainwindow in rf pixels, divide by step to convert to displacement pixels
	#make odd number
		self.strainWindow = int(strainKernelmm/(self.deltaY*stepY) )
		if not self.strainWindow%2:
			self.strainWindow += 1
		self.halfLsq = self.strainWindow//2	
		self.strainwindowmm = self.strainWindow*self.deltaY*self.stepY
			
	######allocate numpy arrays to store quality, dpy, dpx
		self.dpY = np.zeros( (self.numY, self.numX) )
		self.dpX = np.zeros( (self.numY, self.numX) )
		self.quality = np.zeros( (self.numY, self.numX) )
		self.processed = np.zeros( (self.numY, self.numX) )
		self.regionImage = np.zeros( (self.numY, self.numX) )
			
			
	####seed is tuple containing (quality, ptidx, region)
		self.seedList = UniquePriorityQueue()

	#work out number of seed points, and index into dp arrays
		seedsY = int(.05*self.numY)
		if seedsY == 0:
			seedsY = 1
		seedsX = int(.05*self.numX)
		if seedsX == 0:
			seedsX = 1
		strideY = self.numY/seedsY
		strideX = self.numX/seedsX
		self.seedsY = range(0,self.numY, strideY)  
		self.seedsX = range(0,self.numX, strideX)
	#regionarray to determine how many points grown from a particular seed
		self.numRegions = len(self.seedsY)*len(self.seedsX)
		self.regionArray = np.ones( self.numRegions )
	
	#for determining bad seeds
		self.threshold = round( (self.numY/10)*(self.numX/10) )


	def CreateStrainImage(self, preFrame = 0, postFrame = 1, vMax = .02):
		'''With the given parameters, pre, and post RF data create a strain image.'''

		#get pre and post frame
		self.ReadFrame(preFrame)
		self.pre = self.data.copy()
		self.ReadFrame(postFrame)
		self.post = self.data.copy()
		#Perform tracking
		self.TrackSeeds()
		self.TrackNonSeeds()
		#self.DropOutCorrection()
		self.DisplacementToStrain()

		#take the strain image and create a scan-converted strain image.
		startY =self.windowCenterY[self.halfLsq]
		startX =self.windowCenterX[0]
		stepY =self.stepY
		stepX =1

		self.strain[self.strain> vMax] = vMax
		self.strainRGB = self.CreateParametricImage(self.strain,[startY, startX], [stepY, stepX], colormap = 'gray' )

	
	def Sub2ind(self,shape, row, col):
		"""Takes shape as a tuple with two entries"""

		return shape[1]*row + col

	def SubSampleFit(self,f):
		"""This function takes a 3 by 1 array assumes a quadratic function and finds the maximum as if the function were continuous.
	The assumption is that f(x) = ax^2 + bx + c"""

		c = f[1]
		a =  (f[0] + f[2] )/2. - c
		b = -(f[0] - f[2] )/2. 
		delta = -b/(2*a)
		if abs(delta) > 1:
			delta = 0

		return delta

	def AddToSeedList(self,maxCC, y,x, intDpY,intDpX, region ):
		"""Add up to four neighbors to a point to the seed list.  Perform bounds checking to make sure the neighbors aren't
		located out of the image.  Notice that the Normalized cross-correlation goes in as the quality metric, and that the
		sign will be reversed inside the function as necessary. """
				
		#add point above
		if y > 0 and not self.processed[y-1,x]:
			self.seedList.put( (-maxCC, self.Sub2ind(self.dpY.shape, y-1, x), intDpY,intDpX, region)) 


		#add point below
		if y < self.numY - 1 and not self.processed[y+1,x]:
			self.seedList.put((-maxCC, self.Sub2ind(self.dpY.shape, y+1, x), intDpY,intDpX,region) ) 

		#add point to left
		if x > 0 and not self.processed[y,x-1]:
			
			self.seedList.put( (-maxCC, self.Sub2ind(self.dpY.shape, y, x-1),intDpY,intDpX,region) )  

		#add point to right
		if x < self.numX - 1 and not self.processed[y, x+1]:

			self.seedList.put( (-maxCC, self.Sub2ind(self.dpY.shape, y, x+1),intDpY,intDpX,region) )



	def TrackSeeds(self):
		"""Perform block-matching only on the seed points """
		#pdb.set_trace()
	#allocate array to hold results of cross correlation
		resultNp = np.float32( np.zeros( (2*self.rangeY + 1, 2*self.rangeX + 1) ) )
		resultCv = cv.fromarray(resultNp)


		region = 0

	#perform tracking on seed points
		for y in self.seedsY:
			for x in self.seedsX:

				self.processed[y,x] = True
				self.regionImage[y,x] = region

				#GRAB SLICE OF DATA FOR CC CONVERT TO CV FRIENDLY FORMAT
				startBlockY = self.windowCenterY[y] - self.halfY
				stopBlockY = self.windowCenterY[y] + self.halfY + 1
				startBlockX = self.windowCenterX[x] - self.halfX
				stopBlockX = self.windowCenterX[x] + self.halfX + 1
				template = cv.fromarray( np.float32(self.pre[startBlockY:stopBlockY, startBlockX:stopBlockX ] )  )
				
				startBlockY = self.windowCenterY[y] - self.halfY - self.rangeY
				stopBlockY = self.windowCenterY[y] + self.halfY + self.rangeY + 1
				startBlockX = self.windowCenterX[x] - self.halfX - self.rangeX
				stopBlockX = self.windowCenterX[x] + self.halfX + self.rangeX + 1
				image = cv.fromarray( np.float32( self.post[startBlockY:stopBlockY, startBlockX:stopBlockX]  )   )
				
				cv.MatchTemplate(template, image, resultCv, cv.CV_TM_CCORR_NORMED )
				resultNp = np.asarray(resultCv)


				#FIND MAXIMUM, PERFORM SUB-SAMPLE FITTING
				maxInd = resultNp.argmax()
				maxCC = resultNp.max()
				maxY, maxX = np.unravel_index(maxInd, resultNp.shape)
				self.quality[y,x] = maxCC
				
				#perform sub-sample fit
				if maxY > 0 and maxY < 2*self.rangeY - 1:
					deltaY = self.SubSampleFit(resultNp[maxY-1:maxY+2, maxX])	
				else:
					deltaY = 0.0

				
				if maxX > 0 and maxX < 2*self.rangeX - 1:
					deltaX = self.SubSampleFit(resultNp[maxY, maxX-1:maxX + 2] )
				else:
					deltaX = 0.0


				self.dpY[y,x] = maxY - self.rangeY + deltaY
				self.dpX[y,x] = maxX - self.rangeX + deltaX

				intDpY = int(round(self.dpY[y,x]))
				intDpX = int(round(self.dpX[y,x]))

				#PUT ITEM IN SEED LIST
				self.AddToSeedList(maxCC, y,x, intDpY,intDpX, region)

				region +=1


	def TrackNonSeeds(self):
		"""Perform block-matching on all the non-seed points, using a very small search range. """

	#re-allocate array to hold results of cross correlation
		resultNp = np.float32( np.zeros( (2*self.smallRangeY + 1, 2*self.smallRangeX + 1) ) )
		resultCv = cv.fromarray(resultNp)



	#perform tracking on other points
		while self.seedList.qsize() > 0:
		
			(tempQuality, pointInd, iniDpY, iniDpX, region) = self.seedList.get()
			self.regionArray[region] += 1
			(y,x) = np.unravel_index(pointInd, self.dpY.shape)
			self.processed[y,x] = True
			self.regionImage[y,x] = region
			
			#GRAB SLICE OF DATA FOR CC CONVERT TO CV FRIENDLY FORMAT
			startBlockY = self.windowCenterY[y] - self.halfY
			stopBlockY = self.windowCenterY[y] + self.halfY + 1
			startBlockX = self.windowCenterX[x] - self.halfX
			stopBlockX = self.windowCenterX[x] + self.halfX + 1
			template = cv.fromarray( np.float32(self.pre[startBlockY:stopBlockY, startBlockX:stopBlockX ] )  )
			
			startBlockY = self.windowCenterY[y] - self.halfY - self.smallRangeY + iniDpY
			stopBlockY = self.windowCenterY[y] + self.halfY + self.smallRangeY + iniDpY + 1
			startBlockX = self.windowCenterX[x] - self.halfX - self.smallRangeX + iniDpX
			stopBlockX = self.windowCenterX[x] + self.halfX + self.smallRangeX + iniDpX + 1
			image = cv.fromarray( np.float32( self.post[startBlockY:stopBlockY, startBlockX:stopBlockX]  )   )
			
			cv.MatchTemplate(template, image, resultCv, cv.CV_TM_CCORR_NORMED )
			resultNp = np.asarray(resultCv)


			#FIND MAXIMUM, PERFORM SUB-SAMPLE FITTING
			maxInd = resultNp.argmax()
			maxCC = resultNp.max()
			maxY, maxX = np.unravel_index(maxInd, resultNp.shape)
			self.quality[y,x] = maxCC
			
			#perform sub-sample fit
			#fit to f(x) = ax^2 + bx + c in both directions
			if maxY > 0 and maxY < 2*self.smallRangeY - 1:
				deltaY = self.SubSampleFit(resultNp[maxY-1:maxY+2, maxX])	
			else:
				deltaY = 0.0

			
			if maxX > 0 and maxX < 2*self.smallRangeX - 1:
				deltaX = self.SubSampleFit(resultNp[maxY, maxX-1:maxX + 2] )
			else:
				deltaX = 0.0

			#Add in y displacement, keep within allowed range
			self.dpY[y,x] = maxY - self.smallRangeY + deltaY + iniDpY
			if self.dpY[y,x] > self.rangeY:
				self.dpY[y,x] = self.rangeY
			if self.dpY[y,x] < -self.rangeY:
				self.dpY[y,x] = -self.rangeY
		
			#add in x displacement, keep in range	
			self.dpX[y,x] = maxX - self.smallRangeX + deltaX + iniDpX
			if self.dpX[y,x] > self.rangeX:
				self.dpX[y,x] = self.rangeX
			if self.dpX[y,x] < -self.rangeX:
				self.dpX[y,x] = -self.rangeX

			intDpY = int(round(self.dpY[y,x]))
			intDpX = int(round(self.dpX[y,x]))

			#PUT ITEM IN SEED LIST
			self.AddToSeedList(maxCC, y,x, intDpY,intDpX, region) 



	def DropOutCorrection(self):
		"""Re-run the algorithm, but replace the displacements in all the regions not grown from good seeds"""
		self.processed[:] = 0
	#determine bad regions
		goodSeeds = np.zeros( self.numRegions ) 
		goodSeeds[self.regionArray > self.threshold ] = 1 	

		region = 0
	#rerun good seed points
		for y in self.seedsY:
			for x in self.seedsX:

				if goodSeeds[region]:
					self.processed[y,x] = 1

					intDpY = int(round(self.dpY[y,x]))
					intDpX = int(round(self.dpX[y,x]))

					#PUT ITEM IN SEED LIST
					self.AddToSeedList(self.quality[y,x], y,x, intDpY,intDpX, region)



				region += 1

	#re-allocate array to hold results of cross correlation
		resultNp = np.float32( np.zeros( (2*self.smallRangeY + 1, 2*self.smallRangeX + 1) ) )
		resultCv = cv.fromarray(resultNp)

	#rerun algorithm, if point was grown from good seed maintain it

		while self.seedList.qsize() > 0:
		
			#pdb.set_trace()	
			(tempQuality, pointInd, iniDpY, iniDpX, region) = self.seedList.get()
			(y,x) = np.unravel_index(pointInd, self.dpY.shape)
			self.processed[y,x] = 1

			#Re-process if not originally from a good seed
			if not goodSeeds[ self.regionImage[y,x] ]:
				self.regionImage[y,x] = region
			
				#GRAB SLICE OF DATA FOR CC CONVERT TO CV FRIENDLY FORMAT
				startBlockY = self.windowCenterY[y] - self.halfY
				stopBlockY = self.windowCenterY[y] + self.halfY + 1
				startBlockX = self.windowCenterX[x] - self.halfX
				stopBlockX = self.windowCenterX[x] + self.halfX + 1
				template = cv.fromarray( np.float32(self.pre[startBlockY:stopBlockY, startBlockX:stopBlockX ] )  )
				
				startBlockY = self.windowCenterY[y] - self.halfY - self.smallRangeY + iniDpY
				stopBlockY = self.windowCenterY[y] + self.halfY + self.smallRangeY + iniDpY + 1
				startBlockX = self.windowCenterX[x] - self.halfX - self.smallRangeX + iniDpX
				stopBlockX = self.windowCenterX[x] + self.halfX + self.smallRangeX + iniDpX + 1
				image = cv.fromarray( np.float32( self.post[startBlockY:stopBlockY, startBlockX:stopBlockX]  )   )
				
				cv.MatchTemplate(template, image, resultCv, cv.CV_TM_CCORR_NORMED )
				resultNp = np.asarray(resultCv)


				#FIND MAXIMUM, PERFORM SUB-SAMPLE FITTING
				maxInd = resultNp.argmax()
				maxCC = resultNp.max()
				maxY, maxX = np.unravel_index(maxInd, resultNp.shape)
				self.quality[y,x] = maxCC
				
				#perform sub-sample fit
				#fit to f(x) = ax^2 + bx + c in both directions
				if maxY > 0 and maxY < 2*self.smallRangeY - 1:
					deltaY = self.SubSampleFit(resultNp[maxY-1:maxY+2, maxX])	
				else:
					deltaY = 0.0

				
				if maxX > 0 and maxX < 2*self.smallRangeX - 1:
					deltaX = self.SubSampleFit(resultNp[maxY, maxX-1:maxX + 2] )
				else:
					deltaX = 0.0


				self.dpY[y,x] = maxY - self.smallRangeY + deltaY + iniDpY
				if self.dpY[y,x] > self.rangeY:
					self.dpY[y,x] = self.rangeY
				if self.dpY[y,x] < -self.rangeY:
					self.dpY[y,x] = -self.rangeY
				
				self.dpX[y,x] = maxX - self.smallRangeX + deltaX + iniDpX
				if self.dpX[y,x] > self.rangeX:
					self.dpX[y,x] = self.rangeX
				if self.dpX[y,x] < -self.rangeX:
					self.dpX[y,x] = -self.rangeX

				intDpY = int(round(self.dpY[y,x]))
				intDpX = int(round(self.dpX[y,x]))

				#PUT ITEM IN SEED LIST
				self.AddToSeedList(maxCC, y,x, intDpY,intDpX, region)


			else:

				intDpY = int(round(self.dpY[y,x]))
				intDpX = int(round(self.dpY[y,x]))

				#PUT ITEM IN SEED LIST
				self.AddToSeedList(self.quality[y, x], y,x, intDpY,intDpX, region)


	def DisplacementToStrain(self):
		from numpy import zeros, ones, array
		from numpy.linalg import lstsq
		#want to solve equation y = mx + b for m
		#[x1   1      [m     = [y1
		# x2   1       b ]	y2
		# x3   1]               y3]
		#
		#strain image will be smaller than displacement image

		self.strain = zeros( (self.numY - self.strainWindow + 1, self.numX) ) 
		A = ones( (self.strainWindow,2) )
		colOne = array( range(0, self.strainWindow) )
		A[:,0] = colOne
		halfWindow = (self.strainWindow-1)/2 
		self.startYstrain = self.startY + self.stepY*halfWindow 
		self.stopYstrain = self.stopY - self.stepY*halfWindow

		for y in range( halfWindow,self.numY - halfWindow):
			for x in range(self.numX):
				b = self.dpY[y - halfWindow: y + halfWindow + 1, x]
				out = lstsq(A, b)
				xVec = out[0]
				self.strain[y - halfWindow,x] = xVec[0]


		
		self.strain = self.strain/self.stepY
		self.strain = abs(self.strain)
