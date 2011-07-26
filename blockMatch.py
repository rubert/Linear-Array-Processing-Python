#to perform block matching use openCV library

import cv
import numpy as np
import heapq 

def performBlockMatch(pre,post,params):
"""This function takes RF data as input and calculates lateral and axial displacements. """

(points, lines) = pre.shape

####work out boundaries
startY = 0 + params.rangeY + params.halfY
stepY = (1 - params.overlap)*params.windowY
stopY = points - params.halfY #last pixel that can fit a window

startX = params.halfX
stopX = lines - params.halfX

###create arrays containing window Centers in RF data coordinates
windowCenterY = range(startY:stepY:stopY)
numY = len(windowCenterY)
stopY = windowCenterY[-1]

windowCenterX = range(startX:stopX)
numX = len(windowCenterX)

######Allocate numpy arrays to store quality, dpY, dpX
dpY = np.zeros( (numY, numX) )
dpX = np.zeros( (numY, numX) )
quality = np.zeros( (numY, numX) )
processed = np.zeros( (numY, numX) )
regionSet = np.zeros( (numY, numX) )
	
	

seedList = []

#work out number of seed points
seedsY = int(.05*numY)
seedsX = int(.05*numX)
strideY = numY/seedsY*stepY
strideX = numX/seedsX
seedsY = range(startY:strideY:stopY)
seedsX = range(startX:strideX:stopX)


#regionArray to determine how many points grown from a particular seed
regionArray = np.zeros( (numY, numX) )


#allocate array to hold results of cross correlation
resultNp = np.float32( zeros(2*rangeY + 1, 2*rangeX + 1) )
resultCv = cv.fromarray(resultNp)


region = 0

#perform tracking on seed points
for y in range(seedsY):
	for x in range(seedsX):



		#GRAB SLICE OF DATA FOR CC CONVERT TO CV FRIENDLY FORMAT
		startBlockY = windowCenterY[y] - params.halfY
		stopBlockY = windowCenterY[y] + params.halfY
		startBlockX = windowCenterX[x] - params.halfX
		stopBlockX = windowCenterX[x] + params.halfX
		template = cv.fromarray( np.float32(pre[startBlockY:stopBlockY, startBlockX:stopBlockX ] )  )
		
		startBlockY = windowCenterY[y] - params.halfY - params.rangeY
		stopBlockY = windowCenterY[y] + params.halfY + params.rangeY
		startBlockX = windowCenterX[x] - params.halfX - params.rangeX
		stopBlockX = windowCenterX[x] + params.halfX + params.rangeX
		image = cv.fromarray( np.float32( post[startBlockY:stopBlockY, startBlockX:stopBlockX]  )   )
		
		cv.MatchTemplate(template, image, resultCv, cv.CV_TM_CCORR_NORMED )
		resultNp = np.asarray(resultCv)


		#FIND MAXIMUM, PERFORM SUB-SAMPLE FITTING
		maxInd = resultNp.argmax()
		maxCC = resultNp.max()
		maxY, maxX = np.unravel_index(maxCC, resultNp.shape)
		
		#perform sub-sample fit
		if maxY > 0 and maxY < 2*rangeY:
			

		if maxX > 0 and maxX < 2*rangeX:




		#PUT ITEM IN SEED LIST
		heapq.heappush(seedList, (maxCC, maxY, maxX, region) )


		#add point above
		if indY > 0:
			heapq.heappush(seedList, (maxCC, maxY-1, maxX, region) ) 


		#add point below
		if indY < numY:
			heapq.heappush(seedList, (maxCC, maxY + 1, maxX, region) )

		#add point to left
		if indX > 0
			heapq.heappush(seedList, (maxCC, maxY , maxX - 1,region) )

		#add point to right
		if indX < numX:
			heapq.heappush(seedList, (maxCC, maxY, maxX + 1,region) )


		region +=1



#perform tracking on other points
for y in range(seedsY):
	for x in range(seedsX):


		#take point off Seed set, remove seed if it has already been processed
		(temp, indY, indX ) = heapq.heappop(seedList)
		while processed[indY, indX]:
			(temp, indY, indX) = heapq.heappop(seedList)	



		#GRAB SLICE OF DATA FOR CC CONVERT TO CV FRIENDLY FORMAT
		template = cv.fromarray( np.float32(pre[] )  )
		image = cv.fromarray( np.float32( post[] )   )
		cv.MatchTemplate(template, image, resultCv, cv.CV_TM_CCORR_NORMED )
		resultNp = np.asarray(resultCv)


		#FIND MAXIMUM, PERFORM SUB-SAMPLE FITTING
			

		#PUT ITEM IN SEED LIST
		heapq.heappush(seedList, (maxCC, maxY, maxX) )


		#add point above
		if indY > 0 and processed[maxY-1, maxX] = 0:
			heapq.heappush(seedList, (maxCC, maxY-1, maxX) ) 


		#add point below
		if indY < numY and processed[maxY+1, maxX] = 0:
			heapq.heappush(seedList, (maxCC, maxY + 1, maxX) )

		#add point to left
		if indX > 0 and processed[maxY, maxX-1] = 0:
			heapq.heappush(seedList, (maxCC, maxY , maxX - 1) )

		#add point to right
		if indX < numX and processed[maxY, maxX+1] = 0:
			heapq.heappush(seedList, (maxCC, maxY, maxX + 1)





#perform 2nd-pass drop out correction
