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
        elif dist(X0,Y0,x0,y0)>10:
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
        theta = PSFpopt[2]
        b = PSFpopt[3]
        FWHMx, FWHMy = E2moff_toFWHM(ax, ay, b)
        est_area = E2_moff_apsize(ax,ay,b,f=0.9)
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

