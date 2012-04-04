import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import hilbert

class rfClass(object):

    def __init__(self, filename, dataType, centerFreqSimulation = 5.0, sigmaSimulation = 1.0, fsSimulation = 40.0):
        '''Input:
           filename:  The name of the file or directory the data is from
           dataType:  A lowercase string indicating the type of the data, choices are:
                        ei, elasticity imaging on the Seimens S2000
                        rfd, research mode on the Seimens S2000
                        rf, research mode with clinical software on Ultrasonix SonixTouch
                        sim, file format for simulated data with linear array transducers, one frame only
                        multiSim, file format for simulated data with linear array transducers, multiple frames
                         '''
        

        if not( dataType == 'ei' or dataType == 'rfd' or dataType == 'rf' or dataType =='sim' or dataType == 'multiSim'
        or dataType == 'multiFocus', 'wobbler2D'):
            print 'Error.  Datatype must be ei,rfd,rf, sim, multiSim, or multiFocus '
            return

        self.fname = filename
        self.dataType = dataType

        self.soundSpeed = 1540. #m/s
        self.data = None

        if dataType == 'ei':
            import os
            startdir = os.getcwd()
            os.chdir(self.fname)
            
            self.imageType == 'la'
            metaData = open('metaData.txt', 'r')
            metaList = [0,0,0,0]
            for i in range(4):
                metaList[i] = int(metaData.readline() )

            metaData.close()
            nFrames = len( [name for name in os.listdir('.') if os.path.isfile(name)])
            nFrames -= 1 #because of metaData.txt
            points = len( range( metaList[3]/2, metaList[2] + metaList[3]/2 ) )
            aLines = metaList[1]

            self.nFrames = nFrames
            self.soundSpeed = 1540.  #m/s
            self.fs = 40.*10**6  #Hz
            self.deltaX = 40./aLines  #assume a 40 mm FOV
            self.deltaY = self.soundSpeed/(2*self.fs)*10**3
            self.fovX = self.deltaX*(aLines-1)
            self.fovY = self.deltaY*(points -1 )
            os.chdir(startdir)

            self.points = points
            self.lines = aLines

        if dataType == 'rfd':
            import readrfd
            hInfo = readrfd.headerInfo(filename)
            self.fs = hInfo[1]
            if hInfo[4] == 'phased sector':
                self.imageType = 'ps'
                self.deltaX = hInfo[5]
                self.aLineSpacing = hInfo[6]
            else:
                self.imageType = 'la'
                self.deltaX = (1./hInfo[2])*10.
            #read a frame to find out number of points and A lines
            tempFrame = readrfd.readrfd(filename, 1)
            self.nFrames = readrfd.rfdframes(filename)
            self.deltaY = self.soundSpeed/(2*self.fs)*10**3
            self.fovX = self.deltaX*(tempFrame.shape[1] - 1)
            self.fovY = self.deltaY*(tempFrame.shape[0] -1 )
            self.points = tempFrame.shape[0]
            self.lines = tempFrame.shape[1]


        if dataType == 'rf':
            #The header is 19 32-bit integers
            self.imageType = 'la'
            dataFile = open(self.fname, 'r')
            header = np.fromfile(dataFile, np.int32, 19)
            self.header = header
            self.fs = header[15]
            self.deltaY = self.soundSpeed/(2*self.fs)*10**3
            self.points = header[3]
            self.lines = header[2]
            self.deltaX = 40./self.lines #parameter in header is broken, so assume a 40 mm FOV
            self.fovX = self.lines*self.deltaX
            self.fovY = self.points*self.deltaY
            self.nFrames = header[1]
        
        if dataType == 'wobbler2D':
            #The header is 19 32-bit integers

            self.imageType = 'wobbler'
            dataFile = open(self.fname, 'r')
            header = np.fromfile(dataFile, np.int32, 19)
            self.header = header
            self.fs = header[15]
            self.deltaY = self.soundSpeed/(2*self.fs)*10**3
            self.points = header[3]
            self.lines = header[2]
            self.deltaX = 40./self.lines #parameter in header is broken, so assume a 40 mm FOV
            self.fovX = self.lines*self.deltaX
            self.fovY = self.points*self.deltaY
            self.nFrames = header[1]

            self.deltaTheta = np.pi/180*.410
            self.R1 = 39.5

        if dataType == 'sim':

            f = open(filename, 'rb')

            self.imageType == 'la'
            
            #in Hz
            self.freqstep =float( np.fromfile(f, np.double,1) )
            self.points = int( np.fromfile(f, np.int32,1) )
            self.lines = int( np.fromfile(f, np.int32,1) )
            tempReal = np.fromfile(f,np.double, self.points*self.lines )
            tempImag = np.fromfile(f, np.double, self.points*self.lines )
            f.close()

            self.simFs = fsSimulation
            self.sigma = sigmaSimulation
            self.centerFreq = centerFreqSimulation
            self.deltaX = .2  #beamspacing in mm
            self.fovX = self.lines*self.deltaX

        if dataType == 'multiSim':

            f = open(filename, 'rb')

            self.imageType == 'la'
            #in Hz
            self.freqstep =float( np.fromfile(f, np.double,1) )
            self.points = int( np.fromfile(f, np.int32,1) )
            self.lines = int( np.fromfile(f, np.int32,1) )
            self.nFrames = int(np.fromfile(f, np.int32,1) )
            tempReal = np.fromfile(f,np.double, self.points*self.lines )
            tempImag = np.fromfile(f, np.double, self.points*self.lines )
            f.close()

            self.simFs = fsSimulation
            self.sigma = sigmaSimulation
            self.centerFreq = centerFreqSimulation
            self.deltaX = .2  #beamspacing in mm
            self.fovX = self.lines*self.deltaX

        if dataType == 'multiFocus':

            self.imageType == 'la'
            self.fs =  40.0E6
            #self.deltaX = 0.09
            self.deltaX = 0.15
            #read a frame to find out number of points and A lines
            import os
            self.nFrames = len( os.listdir(filename) )
            tempFrame = np.load(filename + '/001.npy')
            self.deltaY = self.soundSpeed/(2*self.fs)*10**3
            self.fovX = self.deltaX*(tempFrame.shape[1] - 1)
            self.fovY = self.deltaY*(tempFrame.shape[0] -1 )
            self.points = tempFrame.shape[0]
            self.lines = tempFrame.shape[1]
       
        if self.imageType == 'ps':
            rPos = self.deltaY*np.arange(self.points)
            thetaPos = self.deltaX*np.arange(-self.lines//2,-self.lines//2 + self.lines)

            THETA, R = np.meshgrid(thetaPos, rPos)

            ##for plotting
            self.X = R*np.sin(THETA) + self.aLineSpacing*np.arange(-self.lines//2,-self.lines//2 + self.lines)
            self.Y = R*np.cos(THETA)
       
        if self.imageType == 'wobbler':
            thetaPos = -self.lines/2*self.deltaTheta + np.arange(self.lines)*self.deltaTheta
            rPos = self.R1 + np.arange(self.points)*self.deltaY
            THETA, R = np.meshgrid(thetaPos, rPos)
            
            self.X = R*np.sin(THETA)
            self.Y = R*np.cos(THETA)
        
        
        self.roiX = [0, self.lines-1]
        self.roiY = [0, self.points - 1]
            
    def __iter__(self):
        self.iterIndex = 0
        return self

    def next(self):
        if self.index == self.nFrames:
            raise StopIteration
        self.ReadFrame(self.iterIndex)
        self.iterIndex += 1
        return self.data.copy()

    def ReadFrame(self, frameNo = 0):
        '''Read a single frame from the input file, the method varies depending on the file type that
        was read in.  This can handle linear array data from the Seimens S2000 recorded in either mode
        or the Ultrasonix scanner recorded using the clinical software.
        For a simulated data file, the 2nd and 3rd parameters matter, not for real data.

        Input::
        frameNo.  The RF data frame from a file obtained off a real ultrasound scanner.
        CenterFreq: (Hz)  The center frequency of the transmitted pulse for simulated ultrasound data.
        Sigma. (Hz)  The standard deviation of the Gaussian shaped transmit pulse for simulated data.
        samplingFrequency. (Hz)  The simulated sampling frequency for a simulated data set.  Higher
        sampling frequencies are achieved through zero-padding.'''

        import os,numpy as np

        if self.dataType == 'ei':
            startdir = os.getcwd()
            os.chdir(self.fname)

            metaData = open('metaData.txt', 'r')
            metaList = [0,0,0,0]
            for i in range(4):
                metaList[i] = int(metaData.readline() )

            metaData.close()
            points = len( range( metaList[3]/2, metaList[2] + metaList[3]/2 ) )
            aLines = metaList[1]
            self.data = np.zeros( (points, aLines) )
            readThisMany = (metaList[2]+ metaList[3]/2)*metaList[1]

            tempIn = open(str(frameNo), 'r')
            frame = fromfile(tempIn,  np.int16, readThisMany)
            frame = frame.reshape( ( metaList[2]+ metaList[3]/2, metaList[1] )  , order = 'F' )
            self.data[:,:] = frame[metaList[3]/2:metaList[3]/2 + metaList[2], :].copy('C')

            os.chdir(startdir)


        if self.dataType == 'rfd':
            import readrfd
            self.data = readrfd.readrfd(self.fname, frameNo + 1)

        if self.dataType == 'rf' or self.dataType == 'wobbler2D':
            #Read in the header, then read in up to one less frame than necessary
            dataFile = open(self.fname, 'r')
            header = np.fromfile(dataFile, np.int32, 19)
            if frameNo == 0:
                temp = np.fromfile( dataFile, np.int16, self.points*self.lines)
                self.data = temp.reshape( (self.points, self.lines),order='F')
            else:
                dataFile.seek( 2*self.points*self.lines*(frameNo-1), 1)
                temp = np.fromfile( dataFile, np.int16, self.points*self.lines)
                self.data = temp.reshape( (self.points, self.lines), order='F')


        if self.dataType == 'sim':
            ####To get the actual sampling frequency we will zero-pad up to the
            ###desired sampling frequency

            f = open(self.fname, 'rb')

            #in Hz
            self.freqstep =float( np.fromfile(f, np.double,1) )
            tmpPoints = int( np.fromfile(f, np.int32,1) )
            self.lines = int( np.fromfile(f, np.int32,1) )
            tempReal = np.fromfile(f,np.double, tmpPoints*self.lines )
            tempImag = np.fromfile(f, np.double, tmpPoints*self.lines )
            f.close()

            pointsToDesiredFs = int(self.simFs*10**6/self.freqstep)
            self.fs = pointsToDesiredFs*self.freqstep
            self.freqData = (tempReal - 1j*tempImag).reshape( (tmpPoints, self.lines), order = 'F' )

            import math
            f = np.arange(0,tmpPoints*self.freqstep, self.freqstep)
            self.pulseSpectrum = np.exp(-(f - self.centerFreq*10**6)*(f - self.centerFreq*10**6)/(2*(self.sigma*10**6)**2) )
            temp = self.freqData*self.pulseSpectrum.reshape( (tmpPoints, 1) )

            zeroPadded = np.zeros( (pointsToDesiredFs, self.lines) ) + 1j*np.zeros( (pointsToDesiredFs, self.lines) )
            zeroPadded[0:tmpPoints, :] = temp
            self.data = np.fft.ifft(zeroPadded,axis=0 ).real
            self.points = pointsToDesiredFs

            self.fs = self.freqstep*self.points
            self.deltaY = 1540. / (2*self.fs)*10**3
            self.fovY = self.deltaY*self.points


        if self.dataType == 'multiSim':
            ####To get the actual sampling frequency we will zero-pad up to the
            ###desired sampling frequency

            f = open(self.fname, 'rb')

            #in Hz
            self.freqstep =float( np.fromfile(f, np.double,1) )
            tmpPoints = int( np.fromfile(f, np.int32,1) )
            self.lines = int( np.fromfile(f, np.int32,1) )
            self.nFrames = int(np.fromfile(f, np.int32, 1) )

            #now jump ahead by as many frames as necessary
            #assume 64 bit, 8 byte doubles
            f.seek(8*2*tmpPoints*self.lines*frameNo, 1)
            #read in desired frame
            tempReal = np.fromfile(f,np.double, tmpPoints*self.lines )
            tempImag = np.fromfile(f, np.double, tmpPoints*self.lines )
            f.close()

            pointsToDesiredFs = int(self.simFs*10**6/self.freqstep)
            self.fs = pointsToDesiredFs*self.freqstep
            self.freqData = (tempReal - 1j*tempImag).reshape( (tmpPoints, self.lines), order = 'F' )

            import math
            f = np.arange(0,tmpPoints*self.freqstep, self.freqstep)
            self.pulseSpectrum = np.exp(-(f - self.centerFreq*10**6)*(f - self.centerFreq*10**6)/(2*(self.sigma*10**6)**2) )
            temp = self.freqData*self.pulseSpectrum.reshape( (tmpPoints, 1) )
            zeroPadded = np.zeros( (pointsToDesiredFs, self.lines) ) + 1j*np.zeros( (pointsToDesiredFs, self.lines) )
            zeroPadded[0:tmpPoints, :] = temp
            self.data = np.fft.ifft(zeroPadded,axis=0 ).real
            self.points = pointsToDesiredFs

            self.fs = self.freqstep*self.points
            self.deltaY = 1540. / (2*self.fs)*10**3
            self.fovY = self.deltaY*self.points

        if self.dataType == 'multiFocus':
            self.data = np.load(self.fname + '/' + str(frameNo + 1).zfill(3) + '.npy')

    def SaveBModeImage(self, fname, image = 0, itkFileName = None, noPng = False):

        """Create a B-mode image, display it, and save it as a Png. The data can also be saved in Itk format."""

        import matplotlib.cm as cm

        self.ReadFrame(image)
        bMode = np.log10(abs(hilbert(self.data, axis = 0)))
        bMode = bMode - bMode.max()
        bMode[ bMode < -3] = -3
        #import matplotlib and create plot
        if not noPng: 
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            
            if self.imageType == 'la':
                ax.imshow(bMode, cmap = cm.gray, vmin = -3, vmax = 0, extent = [0, self.fovX,  self.fovY, 0])
                ax.set_xlabel('Width (mm) ' )
                ax.set_ylabel('Depth (mm) ' )
                ax.set_xticks( np.arange(0, self.fovX, 10) )
                ax.set_yticks( np.arange(0, self.fovY, 10) )
                plt.savefig(fname + '.png')
                plt.close()

            if self.imageType == 'ps' or self.imageType == 'wobbler':

                ax.pcolormesh(self.X,self.Y, bMode, vmin = -3, vmax = 0, cmap = cm.gray)
                ax.set_axis_bgcolor("k")
                plt.ylim(plt.ylim()[::-1])
                plt.savefig(fname + '.png')
                plt.close()
            
        if itkFileName:

            import SimpleITK as sitk 
            
            itkIm = sitk.GetImageFromArray(bMode)
            itkIm.SetSpacing( [self.deltaY, self.deltaX] )
            itkIm.SetOrigin( [0., 0.] )
            
            writer = sitk.ImageFileWriter()
            if itkFileName[-4:] != '.mhd':
                itkFileName += '.mhd'
            writer.SetFileName(itkFileName)
            
            writer.Execute(itkIm)
    
    def DisplayBmodeImage(self, frameNo = 0):
        """Create a B-mode image and immediately display it.  This is useful for interactive viewing of particular
        frames from a data set."""

        self.ReadFrame(frameNo)
        temp = self.data

        #import signal processing modules and generate Numpy array
        bMode = np.log10(abs(hilbert(temp, axis = 0)))
        bMode = bMode - bMode.max()

        #import matplotlib and create plot
        import matplotlib.cm as cm
        fig = plt.figure()
        fig.canvas.set_window_title("B-mode image " )
        ax = fig.add_subplot(1,1,1)
        if self.imageType == 'la':
            ax.imshow(bMode, cmap = cm.gray, vmin = -3, vmax = 0, extent = [0, self.fovX,  self.fovY, 0])

        if self.imageType == 'ps' or self.imageType == 'wobbler':

            ax.pcolormesh(self.X,self.Y, bMode, vmin = -3, vmax = 0, cmap = cm.gray)
            ax.set_axis_bgcolor("k")
            plt.ylim(plt.ylim()[::-1])
        
        plt.show()
    
    def CreateParametricImage(self, paramImage, origin, spacing, inPixels = True, frameNo = 0, colormap = 'jet', vmin = None, vmax = None):
        '''Input:
           paramImage: The image to place within the B-mode image
           origin:  The B-mode pixel at which the upper left hand corner of the parametric image is located.
                        [row, column]
           spacing:  The number of B-mode pixels separating paramImage's pixels
                        [rows, columns]
           colormap:  The colormap of the parametric image'''

        
=======
   
    def ParametricImageResolutionToBmodeResolution(self, paramImage, origin, spacing, inPixels = True, vmin= None, vmax
    = None):
        import numpy as np 
>>>>>>> ff47221521632f6d5effdbc4eb8e0235bcf848df
        from scipy import interpolate

        if not inPixels:
            origin = np.array(origin); spacing = np.array(spacing)
            
            origin[0] /= self.deltaY
            spacing[0] /= self.deltaY
            origin[1] /= self.deltaX
            spacing[1] /= self.deltaX

            origin = origin.round().astype('int')
            spacing = spacing.round().astype('int')

        if vmin:
            paramImage[paramImage < vmin] = vmin

        if vmax:
            paramImage[paramImage > vmax] = vmax
        
        self.paramValMax = paramImage.max()
        self.paramValMin = paramImage.min()
        #work out size of region in B-mode image
        bModeSizeY = spacing[0]*(paramImage.shape[0] - 1) + 1
        bModeSizeX = spacing[1]*(paramImage.shape[1] - 1) + 1

        paramImageUpInY = np.zeros( (bModeSizeY, paramImage.shape[1] ) )
        paramImageUp = np.zeros( (bModeSizeY, bModeSizeX) )

        #Upsample strain image to same spacing as B-mode image.
        paramRfYIndexes = np.arange(origin[0], origin[0] + spacing[0]*paramImage.shape[0], spacing[0] )
        paramRfYIndexesNew = np.arange(origin[0], origin[0] + spacing[0]*(paramImage.shape[0]-1) + 1 )
        paramRfXIndexes = np.arange(origin[1], origin[1] + spacing[1]*paramImage.shape[1], spacing[1] )
        paramRfXIndexesNew = np.arange(origin[1], origin[1] + spacing[1]*(paramImage.shape[1]-1) + 1 )

        #use old image to interpolate
        for x in range(paramImage.shape[1]):
            interp = interpolate.interp1d(paramRfYIndexes, paramImage[:,x] )
            paramImageUpInY[:, x] = interp(paramRfYIndexesNew )

        for y in range(paramImageUp.shape[0]):
            interp = interpolate.interp1d(paramRfXIndexes, paramImageUpInY[y,:] )
            paramImageUp[y, :] = interp( paramRfXIndexesNew )

        return paramImageUp

    def CreateParametricImage(self, paramImage, origin, spacing, inPixels = True, frameNo = 0, colormap = 'jet', vmin = None, vmax = None):
        '''Input:
           paramImage: The image to place within the B-mode image
           origin:  The B-mode pixel at which the upper left hand corner of the parametric image is located.
                        [row, column]
           spacing:  The number of B-mode pixels separating paramImage's pixels
                        [rows, columns]
           colormap:  The colormap of the parametric image'''

         
        self.ParametricImageResolutionToBmodeResolution(paramImage, origin, spacing, inPixels)
 
        #Convert array containing param values to RGBALpha array
        from matplotlib import cm

        palette = cm.ScalarMappable()
        palette.set_cmap(colormap)
        tempIm = palette.to_rgba(paramImageUp)

        #Create B-mode image
        self.ReadFrame(frameNo)
        from scipy.signal import hilbert
        from numpy import log10
        bMode = log10(abs(hilbert(self.data, axis = 0)))
        bMode = bMode - bMode.max()
        bMode[bMode < -3]  = -3
        palette = cm.ScalarMappable()
        palette.set_cmap('gray')
        bMode = palette.to_rgba(bMode)

        bMode[origin[0]:origin[0] + tempIm.shape[0],origin[1]:origin[1] + tempIm.shape[1], :] = tempIm

        return bMode


    def SetRoiFixedSize(self, windowX = 4, windowY = 4, parametricImage = None, origin = None, spacing = None):
        #The Roi size here is given in mm.  For a phased array probe windowX is the
        #number of A-lines in an ROI.


        self.ReadFrame()
        yExtent = self.deltaY*self.points
        xExtent = self.deltaX*self.lines
        bmode = np.log10(abs(hilbert(self.data, axis=0) ) )
        bmode = bmode - bmode.max()
        
        #Compute the center of the bounding box and its size in pixels
        if self.imageType == 'la':        
            if parametricImage == None:
                plt.imshow(bmode, extent = [0,xExtent,yExtent,0], cmap = 'gray', vmin = -3, vmax = 0 )
            else:
                bMode = self.CreateParametricImage(paramImage, origin, spacing)
                plt.imshow(bmode, extent = [0,xExtent,yExtent,0], cmap = 'gray', vmin = -3, vmax = 0 )
 
            x= plt.ginput()
            plt.close()
            boxX = int( x[0][0]/self.deltaX )
            boxY = int( x[0][1]/self.deltaY )
            windowY = int(windowY/self.deltaY)
            windowX = int(windowX/self.deltaX)
            
            xLow = boxX - windowX/2
            if xLow < 0:
                xLow = 0
            xHigh = xLow + windowX
            if xHigh > self.lines:
                xHigh = self.lines
            xLow = xHigh - windowX

            yLow = boxY - windowY/2
            if yLow < 0:
                yLow = 0
            yHigh = yLow + windowY
            if yHigh > self.points+1:
                yHigh = self.points+1
            yLow = yHigh - windowY
       
        #Compute center of bounding box and its size, phased array probe
        if self.imageType == 'ps':
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            if parametricImage == None:
                ax.pcolormesh(self.X, self.Y, bmode, cmap = 'gray', vmin = -3, vmax = 0)
            else:
                ax.pcolormesh(self.X, self.Y, bmode, cmap = 'gray', vmin = -3, vmax = 0, hold = True)
                paramUp = self.ParametricImageResolutionToBmodeResolution(paramImage, origin, spacing)
                xx = self.X[origin[0]:origin[0] + paramUp.shape[0], origin[1]:origin[1] + paramUp.shape[1]]
                yy = self.Y[origin[0]:origin[0] + paramUp.shape[0], origin[1]:origin[1] + paramUp.shape[1]]
                ax.pcolormesh(xx,yy,paramUp)

            ax.set_axis_bgcolor("k")
            plt.ylim(plt.ylim()[::-1])
            x = plt.ginput()
            plt.close() 
           
            #pain in the ass
            #Find X-coordinates of beamlines at depth Y
            xPos = x[0][0]
            depth = x[0][1]

            angles = self.deltaX*np.arange(-self.lines//2,-self.lines//2 + self.lines)
            startX = self.aLineSpacing*np.arange(-self.lines//2,-self.lines//2 + self.lines)
            beamLineXPos = startX + np.sin(angles)*depth
            
            closestLine = abs(beamLineXPos - xPos).argmin()
            closestPoint = abs(np.cos(angles[closestLine])*np.arange(self.points)*self.deltaY  - depth).argmin()

            extentRPixels = int(windowY/self.deltaY)
            
            xLow = closestLine - windowX//2
            if xLow < 0:
                xLow = 0
            xHigh = closestLine + windowX//2
            if xHigh > self.lines:
                xHigh = self.lines
            xLow = xHigh - windowX

            yLow = closestPoint - extentRPixels//2
            if yLow < 0:
                yLow = 0
            yHigh = yLow + extentRPixels
            
            if yHigh > self.points:
                yHigh = self.points
            yLow = yHigh - extentRPixels
        
        if self.imageType == 'wobbler':
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            ax.pcolormesh(self.X, self.Y, bmode, cmap = 'gray', vmin = -3, vmax = 0)
            ax.set_axis_bgcolor("k")
            plt.ylim(plt.ylim()[::-1])
            x = plt.ginput()
            plt.close() 
           
            #pain in the ass
            #Find X-coordinates of beamlines at depth Y
            xPos = x[0][0]
            depth = x[0][1]

            rPos = np.sqrt(xPos**2 + depth**2)
            thetaPos = np.arctan(xPos/depth)

            closestLine = int(thetaPos/self.deltaTheta + (self.lines/2) )
            closestPoint = int( (rPos - self.R1)/self.deltaY )

            extentRPixels = int(windowY/self.deltaY)
            
            xLow = closestLine - windowX//2
            if xLow < 0:
                xLow = 0
            xHigh = closestLine + windowX//2
            if xHigh > self.lines:
                xHigh = self.lines
            xLow = xHigh - windowX

            yLow = closestPoint - extentRPixels//2
            if yLow < 0:
                yLow = 0
            yHigh = yLow + extentRPixels
            
            if yHigh > self.points:
                yHigh = self.points
            yLow = yHigh - extentRPixels

        self.roiX = [xLow, xHigh]
        self.roiY = [yLow, yHigh]

        #Draw Analysis Window
        import matplotlib.cm as cm
        palette = cm.gray
        palette.set_bad('r')
        from numpy import ma, zeros
        bmode = ma.array(bmode)
        bmode.mask = zeros(bmode.shape)


        bmode.mask[yLow:yHigh, xLow] = True
        bmode.mask[yLow:yHigh, xHigh] = True
        #since the pixels are so much thinner in the axial direction we use ten of them to
        #draw the ROI
        lineThickY = 10
        if yLow < lineThickY:
            yLow = lineThickY
        bmode.mask[yLow-lineThickY:yLow, xLow:xHigh] = True

        if yHigh > bmode.mask.shape[0] - lineThickY:
            yHigh = bmode.mask.shape[0] - lineThickY
        bmode.mask[yHigh-lineThickY:yHigh, xLow:xHigh] = True

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        
        if self.imageType == 'la':
            plt.imshow(bmode, cmap = palette, extent = [0, xExtent, yExtent, 0 ],vmin = -3, vmax = 0)
            plt.xticks( np.arange(0,self.fovX, 10) )
            plt.yticks( np.arange(0,self.fovY, 10) )
            plt.xlabel('Width (mm)')
            plt.ylabel('Depth (mm)')
            plt.show()
        
        if self.imageType == 'ps' or self.imageType == 'wobbler':
            ax.pcolormesh(self.X, self.Y,bmode, cmap = palette, vmin = -3, vmax = 0)
            ax.set_axis_bgcolor("k")
            #plt.xticks( np.arange(0,self.fovX, 10) )
            #plt.yticks( np.arange(0,self.fovY, 10) )
            plt.xlabel('Width (mm)')
            plt.ylabel('Depth (mm)')
            plt.ylim(plt.ylim()[::-1])
            plt.show()
   

    def SetRoiBoxSelect(self):
        """
        Do a mouseclick somewhere, move the mouse to some destination, release
        the button. This class gives click- and release-events and also draws
        a line or a box from the click-point to the actual mouseposition
        (within the same axes) until the button is released. Within the
        method 'self.ignore()' it is checked wether the button from eventpress
        and eventrelease are the same.

        """
        #Check that a Data set has been loaded
        if self.fname == None:
            return

        if self.imageType != 'la':
            print 'Function only valid for linear array probes. '
            return

        from matplotlib.widgets import RectangleSelector
        current_ax = plt.subplot(111) # make a new plotingrangej

        def on_select(eclick, erelease):
            self.roiX = [0,0]
            self.roiY = [0,0]

            self.roiX[0] = int(erelease.xdata/self.deltaX)
            self.roiX[1] = int(eclick.xdata/self.deltaX)
            self.roiX.sort()

            self.roiY[0] = int(eclick.ydata/self.deltaY)
            self.roiY[1] = int(erelease.ydata/self.deltaY)
            self.roiY.sort()

        # drawtype is 'box' or 'line' or 'none'
        rectprops = dict(facecolor='red', edgecolor = 'red',
         alpha=0.5, fill=False)
        
        rs = RectangleSelector(current_ax, on_select,
                                       drawtype='box', useblit=True,
                                       button=[1,3], # don't use middle button
                                       minspanx=0, minspany=0,
                                       spancoords='data',
                                       rectprops = rectprops)

        #could be image sequence or just a 2-D image
        import types
        if type(self.data) == types.NoneType:
            self.ReadFrame(0)
            temp = self.data


        from scipy.signal import hilbert
        from numpy import log10
        bMode = log10(abs(hilbert(temp, axis = 0)))
        bMode = bMode - bMode.max()
        bMode[bMode < -3] = -3

        #import matplotlib and create plot
        import matplotlib.cm as cm

        plt.imshow(bMode, cmap = cm.gray, extent = [0, self.fovX, self.fovY, 0])
        plt.show()
        
    def CreateRoiArray(self):
        #Check that a Data set has been loaded
        #could be image sequence or just a 2-D image
        if self.data == None:
            self.ReadFrame(0)
        temp = self.data

        #Create masked array with B-mode data, mask entries for box
        from scipy.signal import hilbert
        from numpy import log10
        bMode = log10(abs(hilbert(temp, axis = 0)))
        bMode = bMode - bMode.max()
        bMode[bMode < -3] = -3

        from numpy import ma, zeros
        bMode = ma.array(bMode)
        bMode.mask = zeros(bMode.shape)
        #replace some B-mode data with sentinels
        xLow = min(self.roiX); xHigh = max(self.roiX)
        yLow = min(self.roiY); yHigh = max(self.roiY)

        bMode.mask[yLow:yHigh, xLow] = True
        bMode.mask[yLow:yHigh, xHigh] = True
        #since the pixels are so much thinner in the axial direction we use ten of them to
        #draw the ROI
        lineThickY = 10
        if yLow < lineThickY:
            yLow = lineThickY
        
        bMode.mask[yLow-lineThickY:yLow, xLow:xHigh] = True

        if yHigh > bMode.mask.shape[0] - lineThickY:
            yHigh = bMode.mask.shape[0] - lineThickY
        
        bMode.mask[yHigh-lineThickY:yHigh, xLow:xHigh] = True

        return bMode

    def ShowRoiImage(self):
        bMode = self.CreateRoiArray()
        #import matplotlib and create plot
        import matplotlib.cm as cm
        palette = cm.gray
        palette.set_bad('r')
        
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        if self.imageType == 'ps' or self.imageType == 'wobbler':
            ax.pcolormesh(self.X, self.Y,bMode, cmap = palette, vmin = -3, vmax = 0)
            ax.set_axis_bgcolor("k")
            #plt.xticks( np.arange(0,self.fovX, 10) )
            #plt.yticks( np.arange(0,self.fovY, 10) )
            plt.xlabel('Width (mm)')
            plt.ylabel('Depth (mm)')
            plt.ylim(plt.ylim()[::-1])
        else:
            im = plt.imshow(bMode, cmap = palette, extent = [0, self.fovX, self.fovY, 0])
            plt.show()



    def SaveRoiImage(self, fname):
        bMode = self.CreateRoiArray()
        #import matplotlib and create plot
        import matplotlib.cm as cm
        palette = cm.gray
        palette.set_bad('r')

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
       
        if self.imageType == 'ps' or self.imageType == 'wobbler':
            ax.pcolormesh(self.X, self.Y,bMode, cmap = palette, vmin = -3, vmax = 0)
            ax.set_axis_bgcolor("k")
            #plt.xticks( np.arange(0,self.fovX, 10) )
            #plt.yticks( np.arange(0,self.fovY, 10) )
            plt.xlabel('Width (mm)')
            plt.ylabel('Depth (mm)')
            plt.ylim(plt.ylim()[::-1])
            if self.imageType == 'wobbler':
                from matplotlib.ticker import MaxNLocator
                plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
        else:
            im = plt.imshow(bMode, cmap = palette, extent = [0, self.fovX, self.fovY, 0])
                    
        plt.savefig(fname)
        plt.close()
