class faranBsc(object):


	 def SLtoBSC(self, x3, initFreq, stepFreq, endFreq, d, sosm = 1540.,soss = 5570.,sosshear = 3374.7,rhom = 1020.,rhos = 2540.,maxang = 180,maxn = 100)
		'''a function to calculate normalized BSC curve from Tony's Faran code. Assume Tony's 
		%code only calculates 180 degree k3absscatlength values.
		%
		%Input:
		%k3absscatlength: the result from Tony's Faran code
		%x3: ka values, it's an array, should be same as the parameters used in 
		%       Tony's result
		%initFreq: the initial frequency (in unit of MHz) to calculate BSC curve
		%stepFreq: the step frequency (MHz)
		%endFreq: the ending frequency (in unit of MHz) to calculate BSC curve
		%d: scatterer size, diameter, not radius (in unit of micro meter)
		%
		%output:
		%freq: frequency range (Hz)
		%BSCCurve: the BSC curve, in unit of 1/m/Sr
		%
		%Hairong Shi, 05/15/2006
		'''
		ksabsscatlength = self.sphere(sosm,soss,sosshear,rhom,rhos,maxang,maxn)
		import numpy
		sos=1540.;   #%speed of sound in surrounding media

		if x3.min()>(initFreq.*1e6)*2*numpy.pi*(d/2*1e-6)/sos or max(x3)<(endFreq*1e6)*2*numpy.pi*(d/2*1e-6)/sos 
		    
		    print 'Frequency out of range'
		    return

		freq=numpy.arange(initFreq, endFreq,stepFreq)*1e6;   #change unit to Hz
		k=2*numpy.pi*freq/sos;
		ka=k*d/2*1e-6;
		temp=interp1(x3, k3absscatlength, ka);  
		BSCCurve=(temp/k)**2;

		return freq, BSCCurve

	 def sphere(self,sosm,soss,sosshear,rhom,rhos,maxang,maxn)

		'''a function which calculates scattering cross sections for spheres
		#according to Faran's theory; the output is a (1+maxang) by (size of x3) matrix
		#where rows correspond to angles (0 to maxang) and columns correspond to different values
		#of ka (step size and max value of ka can be changed by manipulating the
		#initialization of x3); note that the output is the absolute value of the
		#scattering length (i.e. square root of the differential scattering cross
		#section) times the wavenumber in the host medium (output is therefore
		#unitless)

		%sosm = speed of sound in host medium (fluid)
		%soss = speed of sound in sphere
		%sosshear = speed of sound for shear waves in sphere  (if poisson's ratio
		%(sigma)is available, sosshear=soss*sqrt((1-2*sigma)/2/(1-sigma)))
		%rhom = host medium density
		%rhos = sphere density
		%maxang = maximum angle (in degrees) for which to calculate a cross
		%section
		%maxn = number of terms to include in the summation (higher number =
		%greater accuracy but more computation time)

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
		theta=maxang;
		x3=.1:.1:10;  #the range of ka over which cross sections are calculated
		x2=x3*sosm/sosshear;  #ka value for shear waves in sphere
		x1=x3*sosm/soss;  #ka value for compressional waves in sphere

		#main loop over order n
		for n in range(0,maxn)  	    
		    #spherical bessel functions
		    jx1=self.sphbess(n,x1);
		    jx2=self.sphbess(n,x2);
		    jx3=self.sphbess(n,x3);
		    
		    #spherical neumann functions
		    nx3=self.sphneumm(n,x3);
		    
		    #derivatives of spherical bessel and neumann functions
		    jpx3=n/(2*n+1)*self.sphbess(n-1,x3)-(n+1)/(2*n+1)*self.sphbess(n+1,x3);
		    npx3=n/(2*n+1)*self.sphneumm(n-1,x3)-(n+1)/(2*n+1)*self.sphneumm(n+1,x3);
		    jpx1=n/(2*n+1)*self.sphbess(n-1,x1)-(n+1)/(2*n+1)*self.sphbess(n+1,x1);
		    jpx2=n/(2*n+1)*self.sphbess(n-1,x2)-(n+1)/(2*n+1)*self.sphbess(n+1,x2);
		    
		    tandeltax3=-jx3/nx3;
		    tanalphax3=-x3*jpx3/jx3;
		    tanbetax3=-x3*npx3/nx3;
		    tanalphax1=-x1*jpx1/jx1;
		    tanalphax2=-x2*jpx2/jx2;
		    
		    #calculation of tanxsi and coefficient
		    term1=tanalphax1/(tanalphax1+1);
		    term2=(n**2+n)/(n**2+n-1-0.5*x2**2+tanalphax2);
		    term3=(n**2+n-0.5*x2**2+2*tanalphax1)./(tanalphax1+1);
		    term4=((n**2+n)*(tanalphax2+1))/(n**2+n-1-0.5*x2**2+tanalphax2);
		    tanxsi=-x2**2/2*(term1-term2)/(term3-term4);
		    
		    taneta=tandeltax3*(-rhom/rhos*tanxsi+tanalphax3)/(-rhom/rhos*tanxsi+tanbetax3);
		    coeff[n+1,:]=(2*n+1)*taneta/(1+taneta**2)+i*(2*n+1)*taneta**2/(1+taneta**2);
		    
		    #legendre polynomials
		    temp=legendre(n,n,numpy.cos(numpy.pi/180*theta));
		    legend[:,n+1]=temp(1,:);

		k3absscatlength=abs(legend*coeff);  %matrix mult completes summation over n
		return ksabsscatlength


	def sphneumm(self,order,vect)

		#calculates the spherical neumann function of order 'order' for the passed
		#vector
		import numpy
		from scipy.special import yv
		sphneumm=numpy.sqrt(numpy.pi/2./vect)*yv(order+0.5,vect);
		return sphneumm
			
	def sphbess(self,order,vect)
		'''calculates the spherical bessel function of order 'order' for the passed
		vector'''
		import numpy
		from scipy.special import jv

		sphbess=numpy.sqrt(numpy.pi/2./vect)*jv(order+0.5,vect);
