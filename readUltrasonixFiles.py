def readUltrasonixData(fname):
	import numpy as np

	f = open(fname, 'rb')
	header = {'lines': 0 , 'samples': 0, 'fpa': 0, 'numAngles':0, 'degPerAngle':0 }
	#first read header
	header['lines'] = np.fromfile(f, np.int32, 1)  
	header['samples'] = np.fromfile(f, np.int32, 1)  
	header['fpa'] = np.fromfile(f, np.int32, 1)  
	header['numAngles'] = np.fromfile(f, np.int32, 1)  
	header['degPerAngle'] = np.fromfile(f, np.double, 1 )

	frames = np.fromfile(f, np.int16, header['lines']*header['samples']*header['fpa']*header['numAngles'])

	frames = frames.reshape( (header['samples'], header['lines'], header['fpa']*header['numAngles']), order = 'F' )

	return (frames, header)
