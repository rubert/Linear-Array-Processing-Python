from blockMatch import blockMatchParams
from blockMatch import blockMatchClass
from rfData import rfClass

class blockMatchingData(rfClass):

	def __init__(self, fname, ftype):

		super(attenuation, self).__init__(fname, ftype)
	
		#set parameters for block matching
		self.params = blockMatchingParameters()
		self.strainLimit = .02
		self.strain = None
		self.minimumROIBoundaryY = 0 
		self.minimumROIBoundaryX = 0

	
	def CreateStrainImage(self, preFrameNo, postFrameNo):
		'''Takes two frame numbers as input and produces displacement quality and strain images for plotting'''
		self.preFrame = preFrameNo
		from numpy import arange, meshgrid, zeros
		from scipy import interpolate
		from blockMatch import blockMatchClass
	
		self.ReadFrame(preFrameNo)
		pre = self.data.copy()	
		self.ReadFrame(postFrameNo)
		post = self.data.copy()	
		blockMatchInstance = blockMatchClass(pre, post, self.params)
		
		blockMatchInstance.trackSeeds()
		blockMatchInstance.trackNonSeeds()
		blockMatchInstance.displacementToStrain()

		self.startY = blockMatchInstance.startY
		self.stopY = blockMatchInstance.stopY
		self.startX = blockMatchInstance.startX
		self.stopX = blockMatchInstance.stopX

		self.startYStrain = blockMatchInstance.startYstrain 
		self.stopYStrain = blockMatchInstance.stopYstrain
		self.stepY = blockMatchInstance.stepY

		strainRfIndexes = arange(self.startYStrain, self.stopYStrain, self.stepY)
		strainRfIndexesNew = arange(self.startYStrain, self.stopYStrain - self.stepY)
		
		self.strain = self.createParametricImage(blockMatchInstance.strain, origin, spacing, colormap = 'jet'):
		
		if not rfClass.roiX == None:		
			self.strain[0:min(rfClass.roiY), :, :] = bMode[0:min(rfClass.roiY), :, :]
			self.strain[max(rfClass.roiY):-1, :, :] = bMode[max(rfClass.roiY):-1, :, :]
			self.strain[:, 0:min(rfClass.roiX), :] = bMode[:, 0:min(rfClass.roiX), :]
			self.strain[:, max(rfClass.roiX):-1, :] = bMode[:, max(rfClass.roiX):-1, :]
					
		self.fovStrainX = self.deltaX*self.strain.shape[1]
		self.fovStrainY = self.deltaY*self.strain.shape[0] 


	def PlotBmodeStrain(self):
		'''Plot b-mode and strain images together, saves them to output directory.'''

		self.ReadFrame(self.preFrame)	
		temp = self.data
		#import signal processing modules and generate Numpy array
		from scipy.signal import hilbert
		from numpy import log10
		bMode = log10(abs(hilbert(temp, axis = 0)))
		bMode = bMode - bMode.max()

		#import matplotlib and create plot
		import matplotlib.pyplot as plt

		fig = plt.figure()
		ax = fig.add_subplot(121)
		import matplotlib.cm as cm
		plt.imshow(bMode, cmap = cm.gray, vmin = -3, vmax = 0, extent = [0, self.fovX,  self.fovY, 0])
		ax2 = fig.add_subplot(122)	
	
		
		plt.imshow(self.strain, extent = [0, self.fovStrainX, self.fovStrainY, 0] )
		plt.show()

	def WriteBmodeStrain(self, fname):
		self.ReadFrame(self.preFrame)
		temp = self.data

		#import signal processing modules and generate Numpy array
		from scipy.signal import hilbert
		from numpy import log10
		bMode = log10(abs(hilbert(temp, axis = 0)))
		bMode = bMode - bMode.max()

		#import matplotlib and create plot
		import matplotlib.pyplot as plt
		fig = plt.figure()
		ax = fig.add_subplot(121)
		import matplotlib.cm as cm
		plt.imshow(bMode, cmap = cm.gray, vmin = -3, vmax = 0, extent = [0, self.fovX,  self.fovY, 0])
		ax2 = fig.add_subplot(122)	
		
		plt.imshow(self.strain, extent = [0, self.fovStrainX, self.fovStrainY, 0] )
		plt.savefig(fname)
		plt.close()

	def CreateWriteStrainItk(self, preFrameNo, fname):
		'''This function calculates strain for a fiven frame number and writes the file to fname + .mhd'''
		import itk, numpy
	        converter = itk.PyBuffer[itk.Image.F2]	
		writerType = itk.ImageFileWriter[itk.Image.F2]
		writer = writerType.New()
		from numpy import arange, meshgrid, zeros
		from scipy import interpolate
		from blockMatch import blockMatchClass
	
		self.ReadFrame(preFrameNo)
		pre = self.data.copy()	
		postFrameNo = preFrameNo + 1
		self.ReadFrame(postFrameNo)
		post = self.data.copy()	
		blockMatchInstance = blockMatchClass(pre, post, self.params)
		
		blockMatchInstance.trackSeeds()
		blockMatchInstance.trackNonSeeds()
		blockMatchInstance.displacementToStrain()

		self.startY = blockMatchInstance.startY
		self.stopY = blockMatchInstance.stopY
		self.startX = blockMatchInstance.startX
		self.stopX = blockMatchInstance.stopX

		self.startYStrain = blockMatchInstance.startYstrain 
		self.stopYStrain = blockMatchInstance.stopYstrain
		self.stepY = blockMatchInstance.stepY

		strainRfIndexes = arange(self.startYStrain, self.stopYStrain, self.stepY)
		strainRfIndexesNew = arange(self.startYStrain, self.stopYStrain - self.stepY)
		
		self.strain = zeros( (self.data.shape[0], self.data.shape[1] ) )

		'''Upsample strain image to same spacing as B-mode image.'''
		for x in range(self.startX, self.stopX):	
			interp = interpolate.interp1d(strainRfIndexes, blockMatchInstance.strain[:,x - self.startX] )
			self.strain[strainRfIndexesNew, x] = interp(strainRfIndexesNew )

	
		if not rfClass.roiX == None:		
			strain = numpy.float32( abs(self.strain[ min(rfClass.roiY):max(rfClass.roiY), min(rfClass.roiX):max(rfClass.roiX) ] ) )	
		else:	
			strain = numpy.float32(abs(self.strain)	)		
	
		strain[strain > self.strainLimit] = self.strainLimit	
		im = converter.GetImageFromArray(strain)	
		im.SetSpacing( (self.deltaX, self.deltaY) )
		writer.SetFileName(fname)
		writer.SetInput( im )
		writer.Update()
