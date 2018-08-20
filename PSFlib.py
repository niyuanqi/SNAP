#################################################################
# Name:     PSFlib.py                                           #
# Author:   Yuan Qi Ni                                          #
# Version:  November 7, 2017                                    #
# Function: Program contains functions that define various PSF  #
#           shapes for Photometry.                              #
#################################################################

#essential modules
import numpy as np

#function: distance metric on images
def dist(x1, y1, x2, y2):
    #Euclidean distance
    return np.sqrt(np.square(x1-x2)+np.square(y1-y2))

###############################################
# 2D planar fit for sky background in annulus #
###############################################

#function: general 2D background sky plane
def D2plane((x, y), a, b, c):
    return (a*x + b*y + c).ravel()

###############################################
# Circular moffat function for PSF (legacy)   #
###############################################

#function: circular 2D moffat function
def S2moff((x, y), A, a, b, x0, y0):
    m = A*np.power(1+np.square(dist(x,y,x0,y0)/a),-b)
    return m.ravel()
#function: moffat a, b paramters -> fwhm
def S2moff_toFWHM(a,b):
    if b != 0 and np.power(2,1.0/b)-1 > 0:
        return np.abs(a)*(2*np.sqrt(np.power(2,1.0/b)-1))
    else:
        return 0
#function: integrate moffat function
def S2moff_integrate(A,a,b,A_err=0,a_err=None,b_err=None,f=0.9):
    if b != 0 and np.power(2,1.0/b)-1 > 0:
        Int = f*np.pi*np.square(a)*A/(b-1)
        opt_r = moff_aperture(a,b,f)
        if a_err is None or b_err is None:
            return Int, opt_r
        else:
            #calculate errors
            sig = abs(Int)*np.sqrt((2*a_err/a)**2+(A_err/A)**2+(b_err/(b-1))**2)
            return Int, sig, opt_r 
    else:
        return 0
#function: find aperture containing f fraction of moffat light
def S2moff_aperture(a,b,f=0.9):
    return a*np.sqrt(np.power(1 - f,1/(1-b)) - 1)

###############################################
# Elliptical moffat function for PSF          #
###############################################

#function: coordinate counterclockwise rotation by theta
def rot(x, y, theta):
    rad = theta*np.pi/180
    xp = x*np.cos(rad) - y*np.sin(rad)
    yp = x*np.sin(rad) + y*np.cos(rad)
    return xp, yp
#function: coordinate clockwise rotation by theta
def rotI(x, y, theta):
    rad = theta*np.pi/180
    xp = x*np.cos(rad) + y*np.sin(rad)
    yp = -x*np.sin(rad) + y*np.cos(rad)
    return xp, yp
#function: Mahalanobis distance on images
def Mdist(x, y, x0, y0, ax, ay, theta):
    xr, yr = rotI(x-x0, y-y0, theta)
    d = np.sqrt(np.square(xr/ax)+np.square(yr/ay))
    return d
#function: Get elliptical aperture using M distance
def E2ap_get(image, x0, y0, ax, ay, theta, r1, r2):
    rmax = max(f2*rx, f2*ry)
    xaxis = np.arange(max([0,x0-rmax]),min(image.shape[0],x0+rmax+1),dtype=int)
    yaxis = np.arange(max([0,y0-rmax]),min(image.shape[0],y0+rmax+1),dtype=int)
    api = np.array([image[y][x] for x in xaxis for y in yaxis if (Mdist(x,y,x0,y0,ax,ay,theta)<=r2 and Mdist(x,y,x0,y0,ax,ay,theta)>=r1)])
    apx = np.array([x for x in xaxis for y in yaxis if (Mdist(x,y,x0,y0,ax,ay,theta)<=r2 and Mdist(x,y,x0,y0,ax,ay,theta)>=r1)])
    apy = np.array([y for x in xaxis for y in yaxis if (Mdist(x,y,x0,y0,ax,ay,theta)<=r2 and Mdist(x,y,x0,y0,ax,ay,theta)>=r1)])
    return api, apx, apy
#function: elliptical 2D moffat function
def E2moff((x, y), A, ax, ay, b, theta, x0, y0):
    """
    This is the most frequently used function to model PSFs,
    but it is not perfect! Expect residuals!
    Residuals of fit is not a good estimator of noise in PSF.
    It's a systematic error, not a random one.
    In fact it cancels out when comparing integrated PSFs (to ref stars).
    
    Shape is elliptical Moffat with x-radius ax, y-radius ay
    rotated counterclockwise by angle theta in radians
    centered at position x0, y0, with sharpness b, scaling A.
    """
    m = A*np.power(1+np.square(Mdist(x,y,x0,y0,ax,ay,theta)),-b)
    return m.ravel()
#function: integrate elliptical moffat function
def E2moff_integrate(A, ax, ay, b, f=0.9):
    if b > 1:
        return f*A*np.pi*ax*ay/(b-1)
    else:
        return float('Inf')
#function: moffat a, b paramters -> fwhm
def E2moff_toFWHM(ax, ay, b):
    if b != 0 and np.power(2,1.0/b)-1 > 0:
        r = np.sqrt(np.power(2.0, 1/b) - 1)
    else:
        r = 0
    return ax*2*r, ay*2*r
#function: find aperture area containing f fraction of moffat light
def E2moff_apsize(ax,ay,b,f=0.9):
    if b > 1:
        return np.pi*ax*ay*(np.power(1 - f,1/(1-b)) - 1)
    else:
        return float('Inf')
#function: check veracity of PSF fit
def E2moff_verify(PSFpopt, x0=None, y0=None):
    if PSFpopt is None:
        return False
    elif len(PSFpopt) == 7: #moffat PSF
        #extract values from PSF
        A = PSFpopt[0]
        ax = abs(PSFpopt[1])
        ay = abs(PSFpopt[2])
        b = PSFpopt[3]
        theta = PSFpopt[4]
        X0 = PSFpopt[5]
        Y0 = PSFpopt[6]
        FWHMx, FWHMy = E2moff_toFWHM(ax, ay, b)
        est_area = E2moff_apsize(ax,ay,b,f=0.9)
        #check all manner of ridiculous scenarios
        if FWHMx > 20.0 or FWHMx < 1.0:
            #unphysical FWHM
            return False
        elif FWHMy > 20.0 or FWHMy < 1.0:
            #unphysical FWHM
            return False
        elif FWHMx > 2*FWHMy or FWHMx < 0.5*FWHMy:
            #FWHM too skewed
            return False
        elif est_area > np.pi*np.square(3.0*max(FWHMx, FWHMy)):
            #area under PSF domain too large
            return False
        elif dist(X0,Y0,x0,y0)>5:
            #not our source
            return False
        elif b <= 1.0:
            #divergent PSF
            return False
        else:
            #seems legit
            return True
    elif len(PSFpopt) == 4: #moffat shape
        #extract values from PSF
        ax = abs(PSFpopt[0])
        ay = abs(PSFpopt[1])
        b = PSFpopt[2]
        theta = PSFpopt[3]
        FWHMx, FWHMy = E2moff_toFWHM(ax, ay, b)
        est_area = E2moff_apsize(ax,ay,b,f=0.9)
        #check all manner of ridiculous scenarios
        if FWHMx > 20.0 or FWHMx < 1.0:
            #unphysical FWHM
            return False
        elif FWHMy > 20.0 or FWHMy < 1.0:
            #unphysical FWHM
            return False
        elif FWHMx > 2*FWHMy or FWHMx < 0.5*FWHMy:
            #unphysical FWHM
            return False
        elif est_area > np.pi*np.square(3.0*max(FWHMx, FWHMy)):
            #area under PSF domain too large
            return False
        elif b <= 1.0:
            #divergent PSF
            return False
        else:
            #seems legit
            return True

###############################################
# Sersic Profile                              #
###############################################

def b_integrand(r, re, n, C, b):
    """
    This is used in the integration to find b
    """
    return r*Sersic(r, re, n, C, b)

def Sersic(r, re, n, C, b):
    """
    Sersic Function for integration
    """
    ra=np.power((r/re), 1/n)
    return C*np.exp(-b*(ra-1))

def Sersic_integrate(Ie,re,n,e,f=0.9):
    from scipy import integrate
    from scipy.special import gamma

    if n > 0.35:
        #Approximate well using MacArthur 2003
        bn = 2*n - 1./3. + 4./(405.*n) + 46./(25515.*n**2) + 131./(1148175.*n**3) - 2194697./(30690717750.*n**4)
    else:
        #Numerically using Mohammad Akhlaghi 2012
        #This is the value of b for n=0.36
        b_test=0.426200378468
        #rate of increasing b
        b_drate=1.002
        #rate of decreasing b
        b_irate=1.001
        #Effective radius:
        re=2

        while True:
            a=integrate.quad(b_integrand, 0, float('inf'), 
                             args=(re, n, 1, b_test))[0]
            b=integrate.quad(b_integrand, 0, re, 
                             args=(re, n, 1, b_test))[0]
            I_diff=2*np.pi*(a-(2*b))
            if I_diff>-0.00001 and I_diff<0.00001:
                break
            elif I_diff<=-0.00001:
                b_test=b_test/b_drate
            elif I_diff>=0.00001:
                b_test=b_test*b_irate
            bn = b_test
        print "n={} --> b(n)={}".format(n, bn)
    return np.pi*re**2*Ie*2*n*gamma(2*n)*np.power(np.e, bn)/np.power(bn, 2*n)*f
    

###############################################
# Multi-object fit                            #
###############################################
        
#function: Composite moffat psf for multiple objects
def E2moff_multi((x,y), psftype, given, free):
    out = 0
    count = 0
    for i, psf in enumerate(psftype):
        #add moffat to output for each moffat
        if psf == '3':
            #given is empty, general psf params are all in free
            out+= E2moff((x, y),*free[count:count+7])
            count = count+7
        if psf == '2':
            #given contains [ax,ay,b,theta], free has [A, x0, y0]
            out+= E2moff((x,y),free[count],given[i][0],given[i][1],given[i][2],given[i][3],free[count+1],free[count+2])
            count = count+3
        if psf == '1':
            #given contains [ax,ay,b,theta,x0,y0], free has [A]
            out+= E2moff((x, y),free[count],*given[i])
            count = count+1
        if psf[0] == 's':
            #we need to use sersic profile
            from astropy.modeling.models import Sersic2D
            if psf[1] == 'n':
                #full sersic profile
                out+= Sersic2D(*free[count:count+7])(x,y)
                count = count+7
            else:
                #known sersic n
                out+= Sersic2D(amplitude=free[count],r_eff=free[count+1],n=float(psf[1:]),x_0=free[count+2],y_0=free[count+3],ellip=free[count+4],theta=free[count+5])(x,y)
                count = count+6
    return out
