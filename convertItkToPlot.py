def CreateParametricImage(bModeName, parametricName, colormap = 'jet', show = False, outFileName = None,
maxParametricValue = None):
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
        import SimpleITK as sitk

        reader = sitk.ImageFileReader()
        reader.SetFileName(bModeName)
        bModeItk = reader.Execute()
        bMode = sitk.GetArrayFromImage(bModeItk).T
        bModeSpacing = bModeItk.GetSpacing()
        bModeOrigin = bModeItk.GetOrigin()

        reader.SetFileName(parametricName)
        paramImItk = reader.Execute()
        paramImage = sitk.GetArrayFromImage(paramImItk).T
        paramSpacing = paramImItk.GetSpacing()
        paramOrigin = paramImItk.GetOrigin()
       
        if maxParametricValue:
            paramImage[paramImage > maxParametricValue] = maxParametricValue

        print "Mean value in parametric image is: " + str(paramImage.mean() )

        spacing = [ int( round(paramSpacing[0]/bModeSpacing[0]) ), int( round(paramSpacing[1]/bModeSpacing[1]) ) ]
        origin = [int( round(paramOrigin[0]/bModeSpacing[0]) ),int( round(paramOrigin[1]/bModeSpacing[1]) ) ]

        points,lines = bMode.shape

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
            interp = interpolate.interp1d(paramRfYIndexes, paramImage[:,x], kind = 'nearest' )
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
            fig =pyplot.figure()
            ax = fig.add_subplot(111)
            cax = ax.imshow( bMode, extent = [0, bModeSpacing[1]*lines, bModeSpacing[0]*points, 0] , vmin =
            paramImage.min(), vmax = paramImage.max() )
            cbar = fig.colorbar(cax)
            pyplot.show()

if __name__ == '__main__':

    import sys

    if len(sys.argv) <4:
        CreateParametricImage(sys.argv[1], sys.argv[2], show = True)
    else:
        CreateParametricImage(sys.argv[1], sys.argv[2], show = True, maxParametricValue = float(sys.argv[3]) )
