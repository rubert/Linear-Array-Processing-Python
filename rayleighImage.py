from rfData import rfClass
import numpy as np
from scipy.signal import hilbert

class rayleighClass(rfClass):

    def __init__(self, fname, ftype,windowYmm = 8, windowXmm = 8, overlapY = .75, overlapX = .75):


        #Set up RF data object
        super(rayleighClass, self).__init__(fname, ftype)

        #figure out how many windows fit into image
        self.windowYmm = windowYmm
        self.windowXmm = windowXmm
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

        self.stepX = stepX
        self.stepY = stepY

    def ComputeRayleighImage(self, vMax = None, itkFileName = None):

        self.ReadFrame()
        self.rayleighImage =np.zeros( (self.numY, self.numX) )


        #Compute whole image
        counter = 0
        for y, centerY in enumerate(self.winCenterY):
            for x, centerX in enumerate(self.winCenterX):
                tempRegion = self.data[centerY - self.halfY:centerY + self.halfY + 1, centerX - self.halfX:centerX + self.halfX + 1]
                
                envelope = abs(hilbert(tempRegion, axis = 0))
                SNR = envelope.mean()/envelope.std()
                self.rayleighImage[y,x] = abs(SNR - 1.91)
                
                counter += 1
                print str(counter) + " completed of " + str(self.numY*self.numX)


        startY = self.winCenterY[0]
        startX = self.winCenterY[1]
        self.rayleighRGB = self.CreateParametricImage(self.rayleighImage,[startY, startX], [self.stepY, self.stepX], colormap =
                'jet' )
        
