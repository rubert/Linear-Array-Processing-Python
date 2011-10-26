class faranBsc(object):


    def writeBackScatterFile(self,fname,d, deltaFreq = .01,maxFreq = 20, lowCutoff = .1, sosm = 1540.,soss = 5570.,sosshear = 3374.7,rhom = 1020.,rhos = 2540.,maxang = 180,maxn = 25):
        '''This is used to write a backscatter file to be read by the ultrasound simulation program.  The parameters are identical to
        CalculateBSC(), except for the filename.
        Input:
        fname = the name of the output file to be written.
        d = the diameter of the scatterer in microns
        deltaFreq: (MHz)  The frequency spacing of the calculation.
        maxFreq:  (MHz)  The maximum frequency for which the BSC is calculated.
        lowCutoff: (MHz)  The low frequency for which the BSC is calculated.  If this threshold is too low
        calculations may blow up and give NaN.  The program assumes all BSC below lowCutoff are equal to zero.'''


        import numpy
        freq = numpy.arange(0, maxFreq, deltaFreq)
        #make all backscatter coefficients less than .1 MHz equal to 0
        ind_lowCutoff = int(lowCutoff/deltaFreq) + 1

        bscCurve = self.calculateBSC(freq[ind_lowCutoff:], d, sosm = 1540.,soss = 5570.,sosshear = 3374.7,rhom = 1020.,rhos = 2540.,maxang = 180,maxn = 25)
        #Add back in low frequency points by zero padding
        bscCurvePadded = numpy.zeros(len(freq))
        bscCurvePadded[ind_lowCutoff:] = bscCurve
        #Add freq step and number of freq points to backscatter array
        bscCurveExtraInfo = numpy.zeros( len(bscCurvePadded) + 2)
        bscCurveExtraInfo[0] = deltaFreq
        bscCurveExtraInfo[1] = len(freq)
        bscCurveExtraInfo[2:] = bscCurvePadded
        bscCurveExtraInfo.tofile(fname)



    def calculateBSC(self, freq, d, sosm = 1540.,soss = 5570.,sosshear = 3374.7,rhom = 1020.,rhos = 2540.,maxang = 180,maxn = 25):
        '''a function to calculate normalized BSC curve from Tony's Faran code. Assume Tony's
        %code only calculates 180 degree k3absscatlength values.
        %
        %Input:
        %freq:  The frequency over which I would like backscatter coefficients calculated in MHz
        %d: scatterer size, diameter, not radius (in unit of micro meter)
        %
        %output:
        BSCCurve: the BSC curve, in unit of 1/m/Sr
        %
        %Hairong Shi, 05/15/2006
        '''
        import numpy
        freq = freq*1.E6
        #Work out ka values for frequency range and scatterer diameter
        sos=1540.   #%speed of sound in surrounding media
        k=2*numpy.pi*freq/sos
        ka=k*d/2.*1e-6
        #get scattering length and convert to backscatter coefficient
        kabsscatlength = self.sphere(ka,sosm,soss,sosshear,rhom,rhos,maxang,maxn)
        BSCCurve=(kabsscatlength/k)**2

        return BSCCurve

    def sphere(self,ka,sosm,soss,sosshear,rhom,rhos,maxang,maxn):

        '''a function which calculates scattering cross sections for spheres
        according to Faran's theory; the output is a size of x3 vector giving
        the scattering length at 180 degrees
        and entries correspond to different values
        of ka (step size and max value of ka can be changed by manipulating the
        initialization of x3); note that the output is the absolute value of the
        scattering length (i.e. square root of the differential scattering cross
        section) times the wavenumber in the host medium (output is therefore
        unitless)

        ka = an array containing the (wavenumber)*(scatter radius) values
        over which scattering length will be calculated.
        sosm = speed of sound in host medium (fluid)
        soss = speed of sound in sphere
        sosshear = speed of sound for shear waves in sphere  (if poisson's ratio
        (sigma)is available, sosshear=soss*sqrt((1-2*sigma)/2/(1-sigma)))
        rhom = host medium density
        rhos = sphere density
        maxang = maximum angle (in degrees) for which to calculate a cross
        section
        maxn = number of terms to include in the summation (higher number =
        greater accuracy but more computation time)

        %UNITS ARE OF NO CONSEQUENCE AS LONG AS THEY ARE CONSISTENT, I.E.
        %SPEEDS OF SOUND MUST HAVE THE SAME UNITS AND DENSITIES MUST HAVE THE SAME
        %UNITS

        %9/20/04 Tony Gerig

        %Hairong Shi added comments:
        % for glass beads we used in lab, Poisson's ratio sigma=0.21,
        %sosm=1540m/s
        %soss=5570m/s
        %so sosshear=3374.7m/s
        %rhom=1020 kg/m^3
        %rhos=2540 kg/m^3
        %
        %so an example to run this code is:
        %   k3absscatlength = sphere(1540,5570,3374.7,1020,2540,180,100);

        '''
        import numpy
        from scipy.special import lpmn

        #initialization
        theta= 180
        x3=ka
        x2=x3*sosm/sosshear  #ka value for shear waves in sphere
        x1=x3*sosm/soss  #ka value for compressional waves in sphere

        #main loop over order n
        coeff = numpy.zeros((maxn + 1, len(x3)) ) + 1j*numpy.zeros((maxn + 1, len(x3)))
        for n in range(0,maxn + 1):
            #spherical bessel functions
            jx1=self.sphbess(n,x1)
            jx2=self.sphbess(n,x2)
            jx3=self.sphbess(n,x3)

            #spherical neumann functions
            nx3=self.sphneumm(n,x3)

            #derivatives of spherical bessel and neumann functions
            jpx3=n/(2*n+1.)*self.sphbess(n-1,x3)-(n+1)/(2*n+1.)*self.sphbess(n+1,x3)
            npx3=n/(2*n+1.)*self.sphneumm(n-1,x3)-(n+1)/(2*n+1.)*self.sphneumm(n+1,x3)
            jpx1=n/(2*n+1.)*self.sphbess(n-1,x1)-(n+1)/(2*n+1.)*self.sphbess(n+1,x1)
            jpx2=n/(2*n+1.)*self.sphbess(n-1,x2)-(n+1)/(2*n+1.)*self.sphbess(n+1,x2)

            #calculation of tandelta, tanalpha, tanbeta
            tandeltax3=-jx3/nx3
            tanalphax3=-x3*jpx3/jx3
            tanbetax3=-x3*npx3/nx3
            tanalphax1=-x1*jpx1/jx1
            tanalphax2=-x2*jpx2/jx2

            #calculation of tanxsi [eq. 30]
            term1=tanalphax1/(tanalphax1+1)
            term2=(n**2+n)/(n**2+n-1-0.5*x2**2+tanalphax2)
            term3=(n**2+n-0.5*x2**2+2*tanalphax1)/(tanalphax1+1)
            term4=((n**2+n)*(tanalphax2+1))/(n**2+n-1-0.5*x2**2+tanalphax2)
            tanxsi=-x2**2/2*(term1-term2)/(term3-term4)

            #calculation of tanphi [eq. 29]
            tanPhi = -rhom/rhos*tanxsi
            #calculation of coefficient (2n+1)*sin(eta)*cos(eta) part of [eqn. 31]
            taneta=tandeltax3*(tanPhi+tanalphax3)/(tanPhi+tanbetax3)
            coeff[n,:]=(2*n+1)*taneta/(1+taneta**2)+1j*(2*n+1)*taneta**2/(1+taneta**2)

        #legendre polynomials
        temp, deriv=lpmn(n,n,numpy.cos(numpy.pi/180*theta))
        #taking the first row is effectively taking m = 0
        legend=temp[0,:].reshape((1,maxn + 1))

        #matrix mult completes summation over n [eqn. 31]
        k3absscatlength=abs(numpy.dot(legend,coeff))
        output = k3absscatlength.reshape( len(x3) )
        output[numpy.isnan(output)] = 0
        return output

    def sphneumm(self,order,vect):

        #calculates the spherical neumann function of order 'order' for the passed
        #vector
        import numpy
        from scipy.special import yv
        sphneumm=numpy.sqrt(numpy.pi/2./vect)*yv(order+0.5,vect)
        return sphneumm

    def sphbess(self,order,vect):
        '''calculates the spherical bessel function of order 'order' for the passed
        vector'''
        import numpy
        from scipy.special import jv

        sphbess=numpy.sqrt(numpy.pi/2./vect)*jv(order+0.5,vect)
        return sphbess
