#################################################################
# Name:     PSFlib.py                                           #
# Author:   Yuan Qi Ni                                          #
# Version:  November 7, 2017                                    #
# Function: Program contains functions that define various PSF  #
#           shapes for Photometry.                              #
#################################################################

#essential modules
import numpy as np
from scipy.signal import fftconvolve
from scipy.special import gamma, gammainc, gammaincinv
from scipy.optimize import fsolve

#function: distance metric on images
def dist(x1, y1, x2, y2):
    #Euclidean distance
    return np.sqrt(np.square(x1-x2)+np.square(y1-y2))

###############################################
# 2D planar fit for sky background in annulus #
###############################################

#function: general 2D background sky plane
def D2plane((x, y), a, b, c):
    return a*x + b*y + c

###############################################
# Circular moffat function for PSF (legacy)   #
###############################################

#function: circular 2D moffat function
def S2moff((x, y), A, a, b, x0, y0):
    m = A*np.power(1+np.square(dist(x,y,x0,y0)/a),-b)
    return m
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
    costh, sinth = np.cos(rad), np.sin(rad)
    xp = x*costh - y*sinth
    yp = x*sinth + y*costh
    return xp, yp
#function: coordinate clockwise rotation by theta
def rotI(x, y, theta):
    rad = theta*np.pi/180.0
    costh, sinth = np.cos(rad), np.sin(rad)
    xp = x*costh + y*sinth
    yp = -x*sinth + y*costh
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
def E2moff(pos, A, ax, ay, b, theta, x0, y0):
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
    m = A*np.power(1+np.square(Mdist(pos[0],pos[1],x0,y0,ax,ay,theta)),-b)
    return m
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

ker_in = dict()

def KerSersic(pos, PSF, nsamp=1, fker=1):
    global ker_in
    #oversampling factor
    n=nsamp #assume image PSF is well sampled, n=1
    #kernal radius as factor of fwhm
    f=fker #assume PSF is relatively sharp, f=1
    
    hashseq = hash(tuple([pos[0].tostring(),pos[1].tostring(),tuple(PSF)]))
    if hashseq in ker_in:
        #particular subimage and PSF already exists
        kz, sx, sy, ix, iy = ker_in[hashseq]
    else:   
        #psf parameters
        ax, ay, mb, mtheta = PSF
        mA = (mb-1)/(np.pi*ax*ay)
        #kernel radius f*fwhms
        mr = np.sqrt(np.power(2.0, 1/mb) - 1)
        kr = 2*f*max(ax,ay)*mr #in psf frame
        kr = int(np.ceil(np.absolute(kr)))
        #kernel size, nx oversampling
        kn = 2*n*kr+1
        #create kernel base
        kx = np.linspace(-kr, kr, kn, endpoint=True)
        ky = np.linspace(-kr, kr, kn, endpoint=True)
        kx, ky = np.meshgrid(kx, ky)
        #evaluate kernel
        kz = E2moff((kx, ky), mA, ax, ay, mb, mtheta, 0, 0)
        
        #sersic image size 2*fwhm over-extension and 4x oversampling
        xmin, xmax, ymin, ymax = np.min(pos[0]), np.max(pos[0]), np.min(pos[1]), np.max(pos[1])
        nx, ny = (xmax-xmin)*n + 2*n*kr+1, (ymax-ymin)*n + 2*n*kr+1
        sx = np.linspace(xmin-kr, xmax+kr, nx, endpoint=True)
        sy = np.linspace(ymin-kr, ymax+kr, ny, endpoint=True)
        sx, sy = np.meshgrid(sx, sy)
        #indices of original x,y
        ix = np.round((pos[0]-xmin)*n+n*kr).astype(int)
        iy = np.round((pos[1]-ymin)*n+n*kr).astype(int)
        
        #add to global variable
        ker_in[hashseq] = [kz, sx, sy, ix, iy]
    return kz, sx, sy, ix, iy

def Sersic2D((x,y), Ie, re, n, x0, y0, e, theta):
    bn = gammaincinv(2.0*n, 0.5)
    #effective radius scaling
    amax, amin = re, (1-e)*re
    xmaj, xmin = rotI(x-x0, y-y0, theta)
    z = np.sqrt(np.square(xmaj/amax) + np.square(xmin/amin))
    return Ie*np.exp(-bn*(np.power(z, 1./n)-1.))

def SersicK2D((x,y), PSF, Ie, re, n, x0, y0, e, theta, nsamp=1, fker=1):
    #get kernel and base
    kz, sx, sy, ix, iy = KerSersic((x,y), PSF, nsamp=nsamp, fker=fker)
    #evaluate sersic image
    sz = Sersic2D((sx,sy), Ie,re,n,x0,y0,e,theta)
    #convolve sersic image with kernel
    sz = fftconvolve(sz, kz, mode='same')/nsamp**2
    #obtain sz at x,y
    z = sz[tuple(np.array([iy, ix]))]
    #print "eval {0:07.3f}, {1:07.3f}, {2:07.3f}, {3:07.3f}".format(C, re, n, np.max(z))
    return z

def Sersic_integrate(Ie,re,n,e,f=0.9):
    bn = gammaincinv(2.0*n, 0.5)
    return f*(1-e)*2*np.pi*re**2*Ie*n*gamma(2.*n)*np.exp(bn)/np.power(bn, 2*n)

###############################################
# Core-Sersic Profile                         #
###############################################

csb_in = dict()

def CoreSersicb(re, n, gma, rbe):
    global csb_in
    hashseq = str([n,gma,rbe])
    if hashseq in csb_in:
        bn = csb_in[hashseq]
    else:
        rb = rbe*re
        rbn = np.power(rbe, 1./n)
        #func = lambda b: (n/np.power(b, 2*n))*np.exp(b*rbn)*gamma(2.*n)*(1 + gammainc(2.*n, b*rbn) - 2*gammainc(2.*n, b)) - (1./(2.-gma))*np.square(rbe)
        func = lambda b: 1 + gammainc(2.*n, b*rbn) - 2*gammainc(2.*n, b)
        bn = fsolve(func, gammaincinv(2.0*n, 0.5))[0]
        csb_in[hashseq] = bn
    return bn

def CoreSersic2D((x,y), x0, y0, Ib, e, theta, re, n, gma, rbe):
    #evaluate sersic image
    b = CoreSersicb(re, n, gma, rbe)
    a = 10. #sharp transition
    rb = rbe*re #rb as a fraction of re
    #get scaling factor
    C = Ib*np.power(2,-gma/a)*np.exp(b*np.power(2,1./(a*n))*np.power(rbe,1./n))
    #effective radius scaling
    amax, amin = re, (1-e)*re
    xmaj, xmin = rotI(x-x0, y-y0, theta)
    z = np.power(np.square(xmaj/amax) + np.square(xmin/amin), a/2.)
    zb = np.power(rbe, a)
    #spherical core sersic
    rs=np.power(z+zb, 1./(a*n))
    rc=np.power(1+zb/z, gma/a)
    return C*rc*np.exp(-b*rs)

def CoreSersicK2D((x,y), PSF, x0, y0, Ib, e, theta, re, n, gma, rbe, nsamp=1, fker=1):
    """
    Sersic Function for integration
    Trujillo et al. (2004)
    """
    #get kernel and base
    kz, sx, sy, ix, iy = KerSersic((x,y), PSF, nsamp=nsamp, fker=fker)
    #evaluate sersic image
    sz = CoreSersic2D((sx,sy), x0, y0, Ib, e, theta, re, n, gma, rbe)
    #convolve sersic image with kernel
    sz = fftconvolve(sz, kz, mode='same')/nsamp**2
    #obtain sz at x,y
    z = sz[tuple(np.array([iy, ix]))]
    #print "eval {0:07.3f}, {1:07.3f}, {2:07.3f}, {3:07.3f}, {4:07.3f}".format(C, re, n, gma, np.max(z))
    return z

def CoreSersic_integrate(Ib,re,n,gma,rbe,e,f=0.9):
    a = 10.
    rb = rbe*re
    rbn = np.power(rbe, 1./n)
    b = CoreSersicb(re, n, gma, rbe)
    return f*(1-e)*2*np.pi*Ib*(np.square(rb)/(2-gma) + np.exp(b*rbn)*n*(np.square(re)/np.power(b,2*n))*gamma(2.*n)*(1 - gammainc(2.*n, b*rbn)))

###############################################
# Multi-object fit                            #
###############################################

def PSFlen(psf):
    if psf == '3' or psf == '2' or psf =='1':
        return 7
    elif psf[0] == 's':
        if psf[1] == 'n':
            return 7
        else:
            return 3
    elif psf[0] == 'c':
        if psf[1] == 'n':
            return 9
        else:
            return 3
    else:
        print "No PSF in library"

def PSFparams(psf):
    if psf == '3':
        return ['C', 'ax', 'ay', 'b', 'theta', 'x0', 'y0']
    elif psf == '2':
        return ['C', 'x0', 'y0']
    elif psf == '1':
        return ['C']
    elif psf[0] == 's':
        if psf[1] == 'n':
            return ['Ie', 're', 'n', 'x0', 'y0', 'e', 'theta']
        elif psf[-1] == 'f':
            return ['Ie']
        else:
            return ['Ie', 'x0', 'y0']
    elif psf[0] == 'c':
        if psf[1] == 'n':
            return ['x0', 'y0', 'Ib', 'e', 'theta', 're', 'n', 'gma', 'rbe']
        elif psf[-1] == 'f':
            return ['Ib']
        else:
            return ['x0', 'y0', 'Ib']
    elif psf == 'p':
        #sky
        return ['a', 'b', 'c']

#function: Composite moffat psf for multiple objects
def E2moff_multi((x,y), psftype, PSF, given, free, skyflag=0, nsamp=4, fker=4):
    out = 0
    count = 0
    #get kernel and base
    outconv = 0
    kz, sx, sy, ix, iy = KerSersic((x,y), PSF, nsamp=nsamp, fker=fker)
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
            if psf[1] == 'n':
                #make ellipticity an independent parameter
                e, theta = np.absolute(free[count+5]), free[count+6]
                #print "e, theta", e,theta
                Ie, re = free[count], free[count+1]/np.sqrt(1-e)
                params = list([Ie,re])+list(free[count+2:count+5])+list([e, theta])
                #full sersic profile
                outconv+= Sersic2D((sx,sy), *params)
                count = count+7
            elif psf[-1] == 'f' and given != []:
                Ie = free[count]
                n=float(psf.split('-')[0][1:])
                e=float(psf.split('-')[-2])
                re=float(psf.split('-')[-3])/np.sqrt(1-e)
                theta=float(psf.split('-')[-1][:-1])
                
                #sersic with known n and position
                params = [Ie, re, n]+list(given[i])+[e, theta]
                outconv+= Sersic2D((sx,sy), *params)
                count = count+1
            else:
                Ie = free[count]
                n=float(psf.split('-')[0][1:])
                e=float(psf.split('-')[-2])
                re=float(psf.split('-')[-3])/np.sqrt(1-e)
                if psf[-1] == 'f':
                    theta=float(psf.split('-')[-1][:-1])
                else:
                    theta=float(psf.split('-')[-1])

                params = [Ie, re, n]+list(free[count+1:count+3])+[e, theta]
                #known sersic n
                outconv+= Sersic2D((sx,sy), *params)
                count = count+3
        if psf[0] == 'c':
            #we need to use core sersic profile
            if psf[1] == 'n':
                #make ellipticity an independent parameter
                e, theta = np.absolute(free[count+3]), free[count+4]
                Ib, re = free[count+2], free[count+5]/np.sqrt(1-e)
                params = list(free[count:count+2])+list([Ib, e, theta, re])+list(free[count+6:count+9])
                #full sersic profile
                outconv+= CoreSersic2D((sx,sy), *params)
                count = count+9
            elif psf[-1] == 'f' and given != []:
                n=float(psf.split('-')[0][1:])
                gma=float(psf.split('-')[1])
                rbe=float(psf.split('-')[2])
                e=float(psf.split('-')[-2])
                re=float(psf.split('-')[-3])/np.sqrt(1-e)
                theta=float(psf.split('-')[-1][:-1])
                #sersic with known n and position
                params = list(given[i])+list([free[count]])+[e, theta, re, n, gma, rbe]
                outconv+= CoreSersic2D((sx,sy), *params)
                count = count+1
            else:
                n=float(psf.split('-')[0][1:])
                gma=float(psf.split('-')[1])
                rbe=float(psf.split('-')[2])
                e=float(psf.split('-')[-2])
                re=float(psf.split('-')[-3])/np.sqrt(1-e)
                if psf[-1] == 'f':
                    theta=float(psf.split('-')[-1][:-1])
                else:
                    theta=float(psf.split('-')[-1])
                #known sersic n
                params = list(free[count:count+3])+[e, theta, re, n, gma, rbe]
                outconv+= CoreSersic2D((sx,sy), *params)
                count = count+3
    if not (outconv is 0):
        #do convolution
        #convolve sersic image with kernel
        outconv = fftconvolve(outconv, kz, mode='same')/nsamp**2
        out += outconv[tuple(np.array([iy, ix]))]
        #print np.max(outconv), free[15]
        #print np.max(outconv)
    if skyflag:
        out+= D2plane((x, y), *free[count:count+5])
        count = count+5
    return out
        
