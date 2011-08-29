#first define the class unique priority queue
from Queue import PriorityQueue
import heapq
import cv
import numpy as np
import pdb


class UniquePriorityQueue(PriorityQueue):
    def _init(self, maxsize):
#        print 'init'
        PriorityQueue._init(self, maxsize)
        self.values = set()

    def _put(self, item, heappush=heapq.heappush):
#        print 'put',item
        if item[1] not in self.values:
#            print 'uniq',item[1]
            self.values.add(item[1])
            PriorityQueue._put(self, item, heappush)
#        else:  I should check to see if the quality is better, if so, remove old item
#            print 'dupe',item[1]

    def _get(self, heappop=heapq.heappop):
#        print 'get'
        item = PriorityQueue._get(self, heappop)
#        print 'got',item
        self.values.remove(item[1])
        return item


class blockMatchingParameters():

	def __init__(self):	
		self.windowY = 51
		self.halfY = 45
		self.windowX = 5
		self.halfX = 2
		self.rangeY = 45
		self.smallRangeY = 3
		self.rangeX = 8
		self.smallRangeX = 2
		self.overlap = .6
		self.strainWindow = 13


class blockMatchClass():
	
	def __init__(self, pre, post, params = None):
		self.pre = pre
		self.post = post
		(self.points, self.lines) = pre.shape

		if params == None:
			self.params = blockMatchingParameters()
		else:
			self.params = params

	
	####work out boundaries
		startY = 0 + self.params.rangeY + self.params.smallRangeY + self.params.halfY
		stepY =int( (1 - self.params.overlap)*self.params.windowY )
		stopY = self.points - self.params.halfY - self.params.smallRangeY - self.params.rangeY #last pixel that can fit a window

		startX = self.params.halfX + self.params.rangeX + self.params.smallRangeX
		stopX = self.lines - self.params.halfX - self.params.smallRangeX - self.params.rangeX

		self.stepY = stepY
		self.startY = startY
		self.stopY =stopY
		self.startX = startX
		self.stopX = stopX
	###create arrays containing window Centers in RF data coordinates
		self.windowCenterY = range(startY,stopY, stepY)
		self.numY = len(self.windowCenterY)

		self.windowCenterX = range(startX,stopX)
		self.numX = len(self.windowCenterX)

	######Allocate numpy arrays to store quality, dpY, dpX
		self.dpY = np.zeros( (self.numY, self.numX) )
		self.dpX = np.zeros( (self.numY, self.numX) )
		self.quality = np.zeros( (self.numY, self.numX) )
		self.processed = np.zeros( (self.numY, self.numX) )
		self.regionImage = np.zeros( (self.numY, self.numX) )
			
			
	####seed is tuple containing (quality, ptIdx, region)
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
	#regionArray to determine how many points grown from a particular seed
		self.numRegions = len(self.seedsY)*len(self.seedsX)
		self.regionArray = np.ones( self.numRegions )
	
	#for determining bad seeds
		self.threshold = round( (self.numY/10)*(self.numX/10) )

	def sub2ind(self,shape, row, col):
		"""Takes shape as a tuple with two entries"""

		return shape[1]*row + col

	def subSampleFit(self,f):
		"""This function takes a 3 by 1 array assumes a quadratic function and finds the maximum as if the function were continuous.
	The assumption is that f(x) = ax^2 + bx + c"""

		c = f[1]
		a =  (f[0] + f[2] )/2. - c
		b = -(f[0] - f[2] )/2. 
		delta = -b/(2*a)
		if abs(delta) > 1:
			delta = 0

		return delta

	def addToSeedList(self,maxCC, y,x, intDpY,intDpX, region ):
		"""Add up to four neighbors to a point to the seed list.  Perform bounds checking to make sure the neighbors aren't
		located out of the image.  Notice that the Normalized cross-correlation goes in as the quality metric, and that the
		sign will be reversed inside the function as necessary. """
				
		#add point above
		if y > 0 and not self.processed[y-1,x]:
			self.seedList.put( (-maxCC, self.sub2ind(self.dpY.shape, y-1, x), intDpY,intDpX, region)) 


		#add point below
		if y < self.numY - 1 and not self.processed[y+1,x]:
			self.seedList.put((-maxCC, self.sub2ind(self.dpY.shape, y+1, x), intDpY,intDpX,region) ) 

		#add point to left
		if x > 0 and not self.processed[y,x-1]:
			
			self.seedList.put( (-maxCC, self.sub2ind(self.dpY.shape, y, x-1),intDpY,intDpX,region) )  

		#add point to right
		if x < self.numX - 1 and not self.processed[y, x+1]:

			self.seedList.put( (-maxCC, self.sub2ind(self.dpY.shape, y, x+1),intDpY,intDpX,region) )



	def trackSeeds(self):
		"""Perform block-matching only on the seed points """
		#pdb.set_trace()
	#allocate array to hold results of cross correlation
		resultNp = np.float32( np.zeros( (2*self.params.rangeY + 1, 2*self.params.rangeX + 1) ) )
		resultCv = cv.fromarray(resultNp)


		region = 0

	#perform tracking on seed points
		for y in self.seedsY:
			for x in self.seedsX:

				self.processed[y,x] = True
				self.regionImage[y,x] = region

				#GRAB SLICE OF DATA FOR CC CONVERT TO CV FRIENDLY FORMAT
				startBlockY = self.windowCenterY[y] - self.params.halfY
				stopBlockY = self.windowCenterY[y] + self.params.halfY + 1
				startBlockX = self.windowCenterX[x] - self.params.halfX
				stopBlockX = self.windowCenterX[x] + self.params.halfX + 1
				template = cv.fromarray( np.float32(self.pre[startBlockY:stopBlockY, startBlockX:stopBlockX ] )  )
				
				startBlockY = self.windowCenterY[y] - self.params.halfY - self.params.rangeY
				stopBlockY = self.windowCenterY[y] + self.params.halfY + self.params.rangeY + 1
				startBlockX = self.windowCenterX[x] - self.params.halfX - self.params.rangeX
				stopBlockX = self.windowCenterX[x] + self.params.halfX + self.params.rangeX + 1
				image = cv.fromarray( np.float32( self.post[startBlockY:stopBlockY, startBlockX:stopBlockX]  )   )
				
				cv.MatchTemplate(template, image, resultCv, cv.CV_TM_CCORR_NORMED )
				resultNp = np.asarray(resultCv)


				#FIND MAXIMUM, PERFORM SUB-SAMPLE FITTING
				maxInd = resultNp.argmax()
				maxCC = resultNp.max()
				maxY, maxX = np.unravel_index(maxInd, resultNp.shape)
				self.quality[y,x] = maxCC
				
				#perform sub-sample fit
				if maxY > 0 and maxY < 2*self.params.rangeY - 1:
					deltaY = self.subSampleFit(resultNp[maxY-1:maxY+2, maxX])	
				else:
					deltaY = 0.0

				
				if maxX > 0 and maxX < 2*self.params.rangeX - 1:
					deltaX = self.subSampleFit(resultNp[maxY, maxX-1:maxX + 2] )
				else:
					deltaX = 0.0


				self.dpY[y,x] = maxY - self.params.rangeY + deltaY
				self.dpX[y,x] = maxX - self.params.rangeX + deltaX

				intDpY = int(round(self.dpY[y,x]))
				intDpX = int(round(self.dpX[y,x]))

				#PUT ITEM IN SEED LIST
				self.addToSeedList(maxCC, y,x, intDpY,intDpX, region)

				region +=1


	def trackNonSeeds(self):
		"""Perform block-matching on all the non-seed points, using a very small search range. """

	#re-allocate array to hold results of cross correlation
		resultNp = np.float32( np.zeros( (2*self.params.smallRangeY + 1, 2*self.params.smallRangeX + 1) ) )
		resultCv = cv.fromarray(resultNp)



	#perform tracking on other points
		while self.seedList.qsize() > 0:
		
			#pdb.set_trace()	
			(tempQuality, pointInd, iniDpY, iniDpX, region) = self.seedList.get()
			self.regionArray[region] += 1
			(y,x) = np.unravel_index(pointInd, self.dpY.shape)
			self.processed[y,x] = True
			self.regionImage[y,x] = region
			
			#GRAB SLICE OF DATA FOR CC CONVERT TO CV FRIENDLY FORMAT
			startBlockY = self.windowCenterY[y] - self.params.halfY
			stopBlockY = self.windowCenterY[y] + self.params.halfY + 1
			startBlockX = self.windowCenterX[x] - self.params.halfX
			stopBlockX = self.windowCenterX[x] + self.params.halfX + 1
			template = cv.fromarray( np.float32(self.pre[startBlockY:stopBlockY, startBlockX:stopBlockX ] )  )
			
			startBlockY = self.windowCenterY[y] - self.params.halfY - self.params.smallRangeY + iniDpY
			stopBlockY = self.windowCenterY[y] + self.params.halfY + self.params.smallRangeY + iniDpY + 1
			startBlockX = self.windowCenterX[x] - self.params.halfX - self.params.smallRangeX + iniDpX
			stopBlockX = self.windowCenterX[x] + self.params.halfX + self.params.smallRangeX + iniDpX + 1
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
			if maxY > 0 and maxY < 2*self.params.smallRangeY - 1:
				deltaY = self.subSampleFit(resultNp[maxY-1:maxY+2, maxX])	
			else:
				deltaY = 0.0

			
			if maxX > 0 and maxX < 2*self.params.smallRangeX - 1:
				deltaX = self.subSampleFit(resultNp[maxY, maxX-1:maxX + 2] )
			else:
				deltaX = 0.0


			self.dpY[y,x] = maxY - self.params.smallRangeY + deltaY + iniDpY
			if self.dpY[y,x] > self.params.rangeY:
				self.dpY[y,x] = self.params.rangeY
			if self.dpY[y,x] < -self.params.rangeY:
				self.dpY[y,x] = -self.params.rangeY
			
			self.dpX[y,x] = maxX - self.params.smallRangeX + deltaX + iniDpX
			if self.dpX[y,x] > self.params.rangeX:
				self.dpX[y,x] = self.params.rangeX
			if self.dpX[y,x] < -self.params.rangeX:
				self.dpX[y,x] = -self.params.rangeX

			intDpY = int(round(self.dpY[y,x]))
			intDpX = int(round(self.dpX[y,x]))

			#PUT ITEM IN SEED LIST
			self.addToSeedList(maxCC, y,x, intDpY,intDpX, region) 



	def dropOutCorrection(self):
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
					self.addToSeedList(self.quality[y,x], y,x, intDpY,intDpX, region)



				region += 1

	#re-allocate array to hold results of cross correlation
		resultNp = np.float32( np.zeros( (2*self.params.smallRangeY + 1, 2*self.params.smallRangeX + 1) ) )
		resultCv = cv.fromarray(resultNp)

	#rerun algorithm, if point was grown from good seed maintain it
		pdb.set_trace()

		while self.seedList.qsize() > 0:
		
			#pdb.set_trace()	
			(tempQuality, pointInd, iniDpY, iniDpX, region) = self.seedList.get()
			(y,x) = np.unravel_index(pointInd, self.dpY.shape)
			self.processed[y,x] = 1

			#Re-process if not originally from a good seed
			if not goodSeeds[ self.regionImage[y,x] ]:
				self.regionImage[y,x] = region
			
				#GRAB SLICE OF DATA FOR CC CONVERT TO CV FRIENDLY FORMAT
				startBlockY = self.windowCenterY[y] - self.params.halfY
				stopBlockY = self.windowCenterY[y] + self.params.halfY + 1
				startBlockX = self.windowCenterX[x] - self.params.halfX
				stopBlockX = self.windowCenterX[x] + self.params.halfX + 1
				template = cv.fromarray( np.float32(self.pre[startBlockY:stopBlockY, startBlockX:stopBlockX ] )  )
				
				startBlockY = self.windowCenterY[y] - self.params.halfY - self.params.smallRangeY + iniDpY
				stopBlockY = self.windowCenterY[y] + self.params.halfY + self.params.smallRangeY + iniDpY + 1
				startBlockX = self.windowCenterX[x] - self.params.halfX - self.params.smallRangeX + iniDpX
				stopBlockX = self.windowCenterX[x] + self.params.halfX + self.params.smallRangeX + iniDpX + 1
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
				if maxY > 0 and maxY < 2*self.params.smallRangeY - 1:
					deltaY = self.subSampleFit(resultNp[maxY-1:maxY+2, maxX])	
				else:
					deltaY = 0.0

				
				if maxX > 0 and maxX < 2*self.params.smallRangeX - 1:
					deltaX = self.subSampleFit(resultNp[maxY, maxX-1:maxX + 2] )
				else:
					deltaX = 0.0


				self.dpY[y,x] = maxY - self.params.smallRangeY + deltaY + iniDpY
				if self.dpY[y,x] > self.params.rangeY:
					self.dpY[y,x] = self.params.rangeY
				if self.dpY[y,x] < -self.params.rangeY:
					self.dpY[y,x] = -self.params.rangeY
				
				self.dpX[y,x] = maxX - self.params.smallRangeX + deltaX + iniDpX
				if self.dpX[y,x] > self.params.rangeX:
					self.dpX[y,x] = self.params.rangeX
				if self.dpX[y,x] < -self.params.rangeX:
					self.dpX[y,x] = -self.params.rangeX

				intDpY = int(round(self.dpY[y,x]))
				intDpX = int(round(self.dpX[y,x]))

				#PUT ITEM IN SEED LIST
				self.addToSeedList(maxCC, y,x, intDpY,intDpX, region)


			else:

				intDpY = int(round(self.dpY[y,x]))
				intDpX = int(round(self.dpY[y,x]))

				#PUT ITEM IN SEED LIST
				self.addToSeedList(self.quality[y, x], y,x, intDpY,intDpX, region)


	def displacementToStrain(self):
		from numpy import zeros, ones, array
		from numpy.linalg import lstsq
		#want to solve equation y = mx + b for m
		#[x1   1      [m     = [y1
		# x2   1       b ]	y2
		# x3   1]               y3]
		#
		#strain image will be smaller than displacement image

		self.strain = zeros( (self.numY - self.params.strainWindow + 1, self.numX) ) 
		A = ones( (self.params.strainWindow,2) )
		colOne = array( range(0, self.params.strainWindow) )
		A[:,0] = colOne
		halfWindow = (self.params.strainWindow-1)/2 
		self.startYstrain = self.startY + self.stepY*halfWindow 
		self.stopYstrain = self.stopY - self.stepY*halfWindow

		for y in range( halfWindow,self.numY - halfWindow):
			for x in range(self.numX):
				b = self.dpY[y - halfWindow: y + halfWindow + 1, x]
				out = lstsq(A, b)
				xVec = out[0]
				self.strain[y - halfWindow,x] = xVec[0]


		
		self.strain = self.strain/self.stepY
