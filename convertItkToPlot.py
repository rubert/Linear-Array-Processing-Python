def CreateParametricImage(self, bModeName, parametricName, colormap = 'jet', show = False, outFileName = None):
        '''Input:
        bModeName:
        parametricName:
        colormap:
        show:
        outFileName

        If show is true, display the parametric image immediately.
        If outFileName is true, save a PNG image of the output.
         '''

        #read in B-mode image/parametric image
        import SimpleItk as sitk

        reader = sitk.ImageFileReader.New()
        reader.SetFileName(bModeName)
        bModeItk = reader.Execute()
        bMode = sitk.GetArrayFromImage(bModeItk)
        bModeSpacing = bModeItk.GetSpacing()
        bModeOrigin = bModeItk.GetOrigin()

        reader.SetFileName(parametricName)
        paramImItk = reader.Execute()
        paramImage = sitk.GetArrayFromImage(paramImItk)
        paramSpacing = paramImItk.GetSpacing()
        paramOrigin = paramImItk.GetOrigin()
        

        spacing = [ round(paramSpacing[0]/bModeSpacing[0]), round(paramSpacing[1]/bModeSpacing[0]) ]
        origin = [ round(paramOrigin[0]/bModeSpacing[0]), round(paramOrigin[1]/bModeSpacing[1]) ]

        self.paramValMax = paramImage.max()
        self.paramValMin = paramImage.min()

        point,lines = bMode.shape

        #work out size of region in B-mode image
        bModeSizeY = spacing[0]*(paramImage.shape[0] - 1) + 1
        bModeSizeX = spacing[1]*(paramImage.shape[1] - 1) + 1
        
        from numpy import arange,zeros
        from scipy import interpolate

        paramImageUpInY = zeros( (bModeSizeY, paramImage.shape[1] ) )
        paramImageUp = zeros( (bModeSizeY, bModeSizeX) )

        #Upsample strain image to same spacing as B-mode image.
        paramRfYIndexes = arange(origin[0], origin[0] + spacing[0]*paramImage.shape[0], spacing[0] )
        paramRfYIndexesNew = arange(origin[0], origin[0] + spacing[0]*(paramImage.shape[0]-1) + 1 )
        paramRfXIndexes = arange(origin[1], origin[1] + spacing[1]*paramImage.shape[1], spacing[1] )
        paramRfXIndexesNew = arange(origin[1], origin[1] + spacing[1]*(paramImage.shape[1]-1) + 1 )

        #use old image to interpolate
        for x in range(paramImage.shape[1]):
            interp = interpolate.interp1d(paramRfYIndexes, paramImage[:,x] )
            paramImageUpInY[:, x] = interp(paramRfYIndexesNew )


        for y in range(paramImageUp.shape[0]):
            interp = interpolate.interp1d(paramRfXIndexes, paramImageUpInY[y,:] )
            paramImageUp[y, :] = interp( paramRfXIndexesNew )



        '''Convert array containing param values to RGBALpha array'''
        from matplotlib import cm
        from matplotlib import pyplot

        palette = cm.ScalarMappable()
        palette.set_cmap(colormap)
        tempIm = palette.to_rgba(paramImageUp)

        palette = cm.ScalarMappable()
        palette.set_cmap('gray')
        bMode = palette.to_rgba(bMode)
        bMode[origin[0]:origin[0] + tempIm.shape[0],origin[1]:origin[1] + tempIm.shape[1], :] = tempIm

        if show:
            pyplot.imshow( bMode, extent = [0, bModeSpacing[1]*lines, bModeSpacing[0]*points] )
            pyplot.show()

if name == '__main__':

    print "RUnning as a script"
