#first define the class unique priority queue
from Queue import PriorityQueue
import heapq

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

#to perform block matching use openCV library

class blockMatchingParameters():

	def __init__(self):	
		self.windowY = 91
		self.halfY = 45
		self.windowX = 5
		self.halfX = 2
		self.rangeY = 45
		self.smallRangeY = 3
		self.rangeX = 8
		self.smallRangeX = 2
		self.overlap = .6

def sub2ind(shape, row, col):
	"""Takes shape as a tuple with two entries"""

	return shape[1]*row + col

def performBlockMatch(pre,post):
	"""This function takes RF data as input and calculates lateral and axial displacements.  The input arrays should be doubles. """
	

	import cv
	import numpy as np
	import pdb

	#pdb.set_trace()
	params = blockMatchingParameters()

	(points, lines) = pre.shape

####work out boundaries
	startY = 0 + params.rangeY + params.smallRangeY + params.halfY
	stepY =int( (1 - params.overlap)*params.windowY )
	stopY = points - params.halfY - params.smallRangeY - params.rangeY #last pixel that can fit a window

	startX = params.halfX + params.rangeX + params.smallRangeX
	stopX = lines - params.halfX - params.smallRangeX - params.rangeX

###create arrays containing window Centers in RF data coordinates
	windowCenterY = range(startY,stopY, stepY)
	numY = len(windowCenterY)
	stopY = windowCenterY[-1]

	windowCenterX = range(startX,stopX)
	numX = len(windowCenterX)

######Allocate numpy arrays to store quality, dpY, dpX
	dpY = np.zeros( (numY, numX) )
	dpX = np.zeros( (numY, numX) )
	quality = np.zeros( (numY, numX) )
	processed = np.zeros( (numY, numX) )
	regionSet = np.zeros( (numY, numX) )
		
		
####seed is tuple containing (quality, ptIdx, region)
	seedList = UniquePriorityQueue()

#work out number of seed points
	seedsY = int(.05*numY)
	if seedsY == 0:
		seedsY = 1
	seedsX = int(.05*numX)
	if seedsX == 0:
		seedsX = 1
	strideY = numY/seedsY
	strideX = numX/seedsX
#indexing into dp arrays
	seedsY = range(0,numY, strideY)  
	seedsX = range(0,numX, strideX)


#regionArray to determine how many points grown from a particular seed
	regionArray = np.zeros( (numY, numX) )


#allocate array to hold results of cross correlation
	resultNp = np.float32( np.zeros( (2*params.rangeY + 1, 2*params.rangeX + 1) ) )
	resultCv = cv.fromarray(resultNp)


	region = 0

#perform tracking on seed points
	for y in seedsY:
		for x in seedsX:

			processed[y,x] = True

			#GRAB SLICE OF DATA FOR CC CONVERT TO CV FRIENDLY FORMAT
			startBlockY = windowCenterY[y] - params.halfY
			stopBlockY = windowCenterY[y] + params.halfY + 1
			startBlockX = windowCenterX[x] - params.halfX
			stopBlockX = windowCenterX[x] + params.halfX + 1
			template = cv.fromarray( np.float32(pre[startBlockY:stopBlockY, startBlockX:stopBlockX ] )  )
			
			startBlockY = windowCenterY[y] - params.halfY - params.rangeY
			stopBlockY = windowCenterY[y] + params.halfY + params.rangeY + 1
			startBlockX = windowCenterX[x] - params.halfX - params.rangeX
			stopBlockX = windowCenterX[x] + params.halfX + params.rangeX + 1
			image = cv.fromarray( np.float32( post[startBlockY:stopBlockY, startBlockX:stopBlockX]  )   )
			
			cv.MatchTemplate(template, image, resultCv, cv.CV_TM_CCORR_NORMED )
			resultNp = np.asarray(resultCv)


			#FIND MAXIMUM, PERFORM SUB-SAMPLE FITTING
			maxInd = resultNp.argmax()
			maxCC = resultNp.max()
			maxY, maxX = np.unravel_index(maxInd, resultNp.shape)
			quality[y,x] = maxCC
			
			#perform sub-sample fit
			#fit to f(x) = ax^2 + bx + c in both directions
			if maxY > 0 and maxY < 2*params.rangeY - 1:
				
				a = (resultNp[maxY -1, maxX] + resultNp[maxY + 1, maxX] )/2. - resultNp[maxY, maxX]
				b =- (resultNp[maxY - 1, maxX] - resultNp[maxY + 1, maxX] )/2. 
				deltaY = -b/(2*a)

			else:
				deltaY = 0.0

			if maxX > 0 and maxX < 2*params.rangeX - 1:
				a = (resultNp[maxY, maxX - 1] + resultNp[maxY, maxX+ 1] )/2. - resultNp[maxY, maxX]
				b = -(resultNp[maxY, maxX - 1] - resultNp[maxY, maxX+ 1] )/2. 
				deltaX = -b/(2*a)

			else:
				deltaX = 0.0


			dpY[y,x] = maxY - params.rangeY + deltaY
			dpX[y,x] = maxX - params.rangeX + deltaX

			intDpY = int(round(dpY[y,x]))
			intDpX = int(round(dpY[y,x]))

			#PUT ITEM IN SEED LIST
			
			#add point above
			if y > 0 and not processed[y-1,x]:
				seedList.put( (-maxCC, sub2ind(dpY.shape, y-1, x), intDpY,intDpX, region)) 


			#add point below
			if y < numY - 1 and not processed[y+1,x]:
				seedList.put((-maxCC, sub2ind(dpY.shape, y+1, x), intDpY,intDpX,region) ) 

			#add point to left
			if x > 0 and not processed[y,x-1]:
				
				seedList.put( (-maxCC, sub2ind(dpY.shape, y, x-1),intDpY,intDpX,region) )  

			#add point to right
			if x < numX - 1 and not processed[y, x+1]:

				seedList.put( (-maxCC, sub2ind(dpY.shape, y, x+1),intDpY,intDpX,region) )



			region +=1


#re-allocate array to hold results of cross correlation
	resultNp = np.float32( np.zeros( (2*params.smallRangeY + 1, 2*params.smallRangeX + 1) ) )
	resultCv = cv.fromarray(resultNp)



#perform tracking on other points
	
		
	while seedList.qsize() > 0:
	
		#pdb.set_trace()	
		(tempQuality, pointInd, iniDpY, iniDpX, region) = seedList.get()
		(y,x) = np.unravel_index(pointInd, dpY.shape)
		processed[y,x] = True
		
		#GRAB SLICE OF DATA FOR CC CONVERT TO CV FRIENDLY FORMAT
		startBlockY = windowCenterY[y] - params.halfY
		stopBlockY = windowCenterY[y] + params.halfY + 1
		startBlockX = windowCenterX[x] - params.halfX
		stopBlockX = windowCenterX[x] + params.halfX + 1
		template = cv.fromarray( np.float32(pre[startBlockY:stopBlockY, startBlockX:stopBlockX ] )  )
		
		startBlockY = windowCenterY[y] - params.halfY - params.smallRangeY + iniDpY
		stopBlockY = windowCenterY[y] + params.halfY + params.smallRangeY + iniDpY + 1
		startBlockX = windowCenterX[x] - params.halfX - params.smallRangeX + iniDpX
		stopBlockX = windowCenterX[x] + params.halfX + params.smallRangeX + iniDpX + 1
		image = cv.fromarray( np.float32( post[startBlockY:stopBlockY, startBlockX:stopBlockX]  )   )
		
		cv.MatchTemplate(template, image, resultCv, cv.CV_TM_CCORR_NORMED )
		resultNp = np.asarray(resultCv)


		#FIND MAXIMUM, PERFORM SUB-SAMPLE FITTING
		maxInd = resultNp.argmax()
		maxCC = resultNp.max()
		maxY, maxX = np.unravel_index(maxInd, resultNp.shape)
		quality[y,x] = maxCC
		
		#perform sub-sample fit
		#fit to f(x) = ax^2 + bx + c in both directions
		if maxY > 0 and maxY < 2*params.smallRangeY - 1:
			
			a = (resultNp[maxY -1, maxX] + resultNp[maxY + 1, maxX] )/2. - resultNp[maxY, maxX]
			b =- (resultNp[maxY - 1, maxX] - resultNp[maxY + 1, maxX] )/2. 
			deltaY = -b/(2*a)

		else:
			deltaY = 0.0

		if maxX > 0 and maxX < 2*params.smallRangeX - 1:
			a = (resultNp[maxY, maxX - 1] + resultNp[maxY, maxX+ 1] )/2. - resultNp[maxY, maxX]
			b = -(resultNp[maxY, maxX - 1] - resultNp[maxY, maxX+ 1] )/2. 
			deltaX = -b/(2*a)

		else:
			deltaX = 0.0


		dpY[y,x] = maxY - params.smallRangeY + deltaY + iniDpY
		if dpY[y,x] > params.rangeY:
			dpY[y,x] = params.rangeY
		if dpY[y,x] < -params.rangeY:
			dpY[y,x] = -params.rangeY
		
		dpX[y,x] = maxX - params.smallRangeX + deltaX + iniDpX
		if dpX[y,x] > params.rangeX:
			dpX[y,x] = params.rangeX
		if dpX[y,x] < -params.rangeX:
			dpX[y,x] = -params.rangeX

		intDpY = int(round(dpY[y,x]))
		intDpX = int(round(dpX[y,x]))

		#PUT ITEM IN SEED LIST
		
		#add point above
		if y > 0 and not processed[y-1,x]:
			seedList.put( (-maxCC, sub2ind(dpY.shape, y-1, x), intDpY,intDpX, region)) 


		#add point below
		if y < numY - 1 and not processed[y+1,x]:
			seedList.put((-maxCC, sub2ind(dpY.shape, y+1, x), intDpY,intDpX,region) ) 

		#add point to left
		if x > 0 and not processed[y,x-1]:
			
			seedList.put( (-maxCC, sub2ind(dpY.shape, y, x-1),intDpY,intDpX,region) )  

		#add point to right
		if x < numX - 1 and not processed[y, x+1]:

			seedList.put( (-maxCC, sub2ind(dpY.shape, y, x+1),intDpY,intDpX,region) )



	return (dpY, dpX, quality)
