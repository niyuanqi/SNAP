#################################################################
# Name:     Photometry.py                                       #
# Author:   Yuan Qi Ni                                          #
# Version:  August 25, 2016                                     #
# Function: Program contains functions that perform essential   #
#           photometric tasks.                                  #
#################################################################

#essential modules
import numpy as np

#maximum fev for curve_fit
maxfev = 10000

#function: distance metric on images
def dist(x1, y1, x2, y2):
    #Euclidean distance
    return np.sqrt(np.square(x1-x2)+np.square(y1-y2))

#function: photometric aperture at (x0,y0) from r1 to r2
def ap_get(image, x0, y0, r1, r2):
    xaxis = np.arange(max([0,x0-r2]), min(image.shape[1],x0+r2+1), dtype=int)
    yaxis = np.arange(max([0,y0-r2]), min(image.shape[0],y0+r2+1), dtype=int)
    api = np.array([image[y][x] for x in xaxis for y in yaxis if (dist(x0,y0,x,y)<=r2 and dist(x0,y0,x,y)>=r1)])
    apx = np.array([x for x in xaxis for y in yaxis if (dist(x0,y0,x,y)<=r2 and dist(x0,y0,x,y)>=r1)])
    apy = np.array([y for x in xaxis for y in yaxis if (dist(x0,y0,x,y)<=r2 and dist(x0,y0,x,y)>=r1)])
    return api, apx, apy

#function: photometric aperture around multiple sources from r1 to r2
def ap_multi(image, x0, y0, r1, r2):
    Nobj = len(x0)
    #Extract zone around all objects
    xaxis = np.arange(0,image.shape[1], dtype=int)
    yaxis = np.arange(0,image.shape[0], dtype=int)
    xmask = np.zeros(image.shape[1])
    ymask = np.zeros(image.shape[0])
    for i in range(Nobj):
        xap_single = np.absolute(xaxis-x0[i]) <= r2
        xmask = np.logical_or(xmask, xap_single)
        yap_single = np.absolute(yaxis-y0[i]) <= r2
        ymask = np.logical_or(ymask, yap_single)
    xaxis = xaxis[xmask]
    yaxis = yaxis[ymask]
    #find union of all apertures in zone
    apx, apy, api = [], [], []
    for x in xaxis:
        for y in yaxis:
            #is this pixel in an aperture?
            inap = False
            for i in range(Nobj):
                if dist(x0[i],y0[i],x,y)<=r2:
                    inap = True
            for i in range(Nobj):
                if dist(x0[i],y0[i],x,y)<r1:
                    inap = False
            if inap:
                apx.append(x)
                apy.append(y)
                api.append(image[y][x])
    return np.array(api), np.array(apx), np.array(apy)
    
#function: clean out cosmic rays and junk from PSF
def PSFclean(x,y,psf,ref,skyN=None,sat=40000,f=10):
    #remove saturated pixels
    mask = psf<sat
    if skyN is not None:
        #remove pixels that are 10sigma below or above fit (dead? hot?)
        mask1 = np.absolute(psf-ref)<f*np.sqrt(np.absolute(ref)+skyN**2)
        mask = np.logical_and(mask, mask1)
    #remove pixels that are a result of masking
    mask1 = np.absolute(psf) > 2e-30
    mask = np.logical_and(mask, mask1)
    return x[mask], y[mask], psf[mask]

#function: measure saturation level of CCD
def satpix(image):
    #Make histogram of pixel values
    ns, bins = np.histogram(image, bins=100)
    bins = (bins[1:]-bins[:-1])/2.0 + bins[:-1]
    #Half of max pixel value
    half = bins[-1]/2.0
    mask = bins > half
    #Determine pixel full well capacity
    full = bins[mask][np.argmax[ns[mask]]]
    #Return saturation level
    return full/2.0

#function: fits background sky plane and noise
def SkyFit(image, x0, y0, fwhm=5.0, sat=40000.0, verbosity=0):
    #update 180610: x0, y0 need to be lists (even if length is 1)
    
    from scipy.optimize import curve_fit
    from PSFlib import D2plane
    from MagCalc import PSFError
    
    #get background sky annulus
    #inner_annulus, inner_x, inner_y = ap_get(image, x0, y0, 4*fwhm, 5*fwhm)
    #outer_annulus, outer_x, outer_y = ap_get(image, x0, y0, 6*fwhm, 7*fwhm)
    inner_annulus, inner_x, inner_y = ap_multi(image, x0, y0, 4*fwhm, 5*fwhm)
    inner_x, inner_y, inner_annulus = PSFclean(inner_x,inner_y,inner_annulus,inner_annulus,sat=sat)
    outer_annulus, outer_x, outer_y = ap_multi(image, x0, y0, 6*fwhm, 7*fwhm)
    outer_x, outer_y, outer_annulus = PSFclean(outer_x,outer_y,outer_annulus,outer_annulus,sat=sat)
    #get first estimate mean background value
    innerB = np.mean(inner_annulus)
    outerB = np.mean(outer_annulus)
    #take outer annulus if inner annulus probably contains a star
    skyB = outerB if outerB<innerB else innerB
    skyi, skyx, skyy = (outer_annulus, outer_x, outer_y) if outerB<innerB else (inner_annulus, inner_x, inner_y)
    #fit sky background
    try:
        skypopt, skypcov = curve_fit(D2plane, (skyx, skyy), skyi, p0=[0,0,skyB], maxfev=maxfev, absolute_sigma=True)
        #Fit function
        skyTheo = D2plane((skyx,skyy),*skypopt)
        skyN = np.std(skyi-skyTheo)
        #filter out noisy pixels at 5sigma level (cos rays/hot pix)
        skyx, skyy, skyi = PSFclean(skyx,skyy,skyi,skyTheo,skyN,sat,10)
        
        #calculate better fit from cleaner data
        skypopt, skypcov = curve_fit(D2plane, (skyx, skyy), skyi, p0=skypopt, maxfev=maxfev, absolute_sigma=True)
        try:
            #try to calculate fit error
            skyperr = np.sqrt(np.diag(skypcov))
        except:
            #fit error uncalculable
            try:
                #take closer initial conditions, try again
                skypopt, skypcov = curve_fit(D2plane, (skyx, skyy), skyi, p0=skypopt, maxfev=maxfev, absolute_sigma=True)
                skyperr = np.sqrt(np.diag(skypcov))
            except:
                #fit error really is uncalculable, how???
                raise PSFError('Unable to fit sky.')
        #calculate sky noise near source
        skyTheo = D2plane((skyx,skyy),*skypopt)
        skyN = np.std(skyi-skyTheo)
        #calculate goodness of fit
        skyX2dof = np.square((skyi-skyTheo)/skyN)/(len(skyi)-3)
    except:
        #catastrophic failure of sky plane fitting, How???
        raise PSFError('Sky fitting catastrophic failure.')
    if verbosity > 0:
        print "sky plane fit parameters"
        print "[a, b, c] = "+str(skypopt)
        print "errors = "+str(skyperr)
    #return sky plane, sky noise, and goodness of fit
    return skypopt, skyperr, skyX2dof, skyN

#function: extracts PSF from source
def PSFextract(image, x0, y0, fwhm=5.0, fitsky=True, sat=40000.0, verbosity=0):
    
    from scipy.optimize import curve_fit
    from PSFlib import D2plane, E2moff, E2moff_toFWHM, E2moff_verify
    
    #fit sky background in an annulus
    skypopt, skyperr, skyX2dof, skyN = SkyFit(image, [x0], [y0], fwhm, sat, verbosity)
    
    #get fit box to fit psf
    fsize = 3
    intens, x, y = ap_get(image, x0, y0, 0, fsize*fwhm)
    #get an approximate fix on position
    x0 = np.sum(intens*x)/intens.sum()
    y0 = np.sum(intens*y)/intens.sum()
    #get centered fit box
    intens, x, y = ap_get(image, x0, y0, 0, fsize*fwhm)

    if fitsky:
        #get sky background
        sky = D2plane((x,y),*skypopt)
        #subtract sky background
        intens = intens - sky

    #filter out saturated pixels
    x, y, intens = PSFclean(x,y,intens,intens,skyN,sat,10)
    
    try:
        #fit 2d psf to background subtracted source light
        est = [image[int(y0)][int(x0)],fwhm/4.0,fwhm,3.0,0.0,x0,y0]
        bounds = ([-float("Inf"),0.01,0.01,1.01,0.0,0.0,0.0],[float("Inf"),2*fwhm,2*fwhm,float("Inf"),179.99,image.shape[0],image.shape[1]])
        PSFpopt, PSFpcov = curve_fit(E2moff, (x, y), intens, sigma=np.sqrt(np.absolute(intens)+skyN**2), p0=est, bounds=bounds, absolute_sigma=True, maxfev=maxfev)
        #Fit function
        I_theo = E2moff((x,y),*PSFpopt)
        #filter out noisy pixels at 5sigma level (cos rays/hot pix)
        x, y, intens = PSFclean(x,y,intens,I_theo,skyN,sat,10)

        #calculate better PSF from cleaner data
        PSFpopt, PSFpcov = curve_fit(E2moff, (x, y), intens, sigma=np.sqrt(np.absolute(intens)+skyN**2) , p0=PSFpopt, bounds=bounds, absolute_sigma=True, maxfev=maxfev)
        try:
            #try to calculate fit error
            PSFperr = np.sqrt(np.diag(PSFpcov))
        except:
            try:
                #take closer initial conditions
                PSFpopt, PSFpcov = curve_fit(E2moff, (x, y), intens, sigma=np.sqrt(np.absolute(intens)+skyN**2) , p0=PSFpopt, bounds=bounds, absolute_sigma=True, maxfev=maxfev)
                PSFperr = np.sqrt(np.diag(PSFpcov))
            except:
                PSFperr = [0]*5
        
        #calculate goodness of fit
        I_theo = E2moff((x, y),*PSFpopt)
        X2dof = np.sum(np.square((intens-I_theo)/np.sqrt(np.absolute(intens)+skyN**2)))/(len(intens)-len(PSFpopt))
    except:
        #catastrophic failure of PSF fitting
        print "PSF fitting catastrophic failure"
        PSFpopt = [0]*7
        PSFperr = [0]*7
        X2dof = 0
    
    #get values from fit
    A = PSFpopt[0]
    ax = abs(PSFpopt[1])
    ay = abs(PSFpopt[2])
    b = PSFpopt[3]
    theta = PSFpopt[4]
    X0 = PSFpopt[5]
    Y0 = PSFpopt[6]
    PSFpopt = [A, ax, ay, b, theta, X0, Y0]
    FWHMx, FWHMy = E2moff_toFWHM(ax, ay, b)
    
    if verbosity > 0:
        print "PSF moffat fit parameters"
        print "[A,ax,ay,b,theta,X0,Y0] = "+str(PSFpopt)
        print "parameter errors = "+str(PSFperr)
        print "Chi2 = "+str(X2dof)
        print "[FWHMx', FWHMy]' = "+str([FWHMx, FWHMy])
        
    #graph fits if verbosity is high enough
    if verbosity > 1:
        PSF_plot(image, x0, y0, PSFpopt, X2dof, skypopt, skyN, fitsky, fsize*fwhm)
    #check, if fit is ridiculous or noise object: give no fit
    if E2moff_verify(PSFpopt, x0, y0) and PSFpopt[0]>5*skyN:
        return PSFpopt, PSFperr, X2dof, skypopt, skyN
    else:
        return [0]*7, [0]*7, 0, [0]*3, skyN
    
#function: fits PSF to source
def PSFfit(image, PSF, PSFerr, x0, y0, fitsky=True, sat=40000.0, verbosity=0):

    from scipy.optimize import curve_fit
    from PSFlib import D2plane, E2moff, E2moff_toFWHM, E2moff_verify

    #get given fit parameters
    ax, axerr = PSF[0], PSFerr[0]
    ay, ayerr = PSF[1], PSFerr[1]
    b, berr = PSF[2], PSFerr[2]
    theta, thetaerr = PSF[3], PSFerr[3]
    FWHMx, FWHMy = E2moff_toFWHM(ax, ay, b)
    fwhm = max(FWHMx, FWHMy)

    #fit sky background in an annulus
    skypopt, skyperr, skyX2dof, skyN = SkyFit(image, [x0], [y0], fwhm, sat, verbosity=0)

    #get fit box to fit psf
    fsize = 3
    intens, x, y = ap_get(image, x0, y0, 0, fsize*fwhm)
    #get an approximate fix on position
    x0 = np.sum(intens*x)/intens.sum()
    y0 = np.sum(intens*y)/intens.sum()
    #get centered fit box
    intens, x, y = ap_get(image, x0, y0, 0, fsize*fwhm)
    
    if fitsky:
        #get sky background
        sky = D2plane((x,y),*skypopt)
        #subtract sky background
        intens = intens - sky

    #filter out saturated pixels
    x, y, intens = PSFclean(x,y,intens,intens,skyN,sat,10)
    
    try:
        #fit 2d fixed psf to background subtracted source light
        est = [image[int(y0)][int(x0)],x0,y0]
        bounds = ([-float("Inf"),0,0],[float("Inf"),image.shape[0],image.shape[1]])
        fitpopt, fitpcov = curve_fit(lambda (x, y),A,x0,y0: E2moff((x, y),A,ax,ay,b,theta,x0,y0), (x,y), intens, sigma=np.sqrt(np.absolute(intens)+skyN**2), p0=est, bounds=bounds, absolute_sigma=True, maxfev=maxfev)
        #parameters fitted to source
        A = fitpopt[0]
        X0 = fitpopt[1]
        Y0 = fitpopt[2]
        PSFpopt = [A,ax,ay,b,theta,X0,Y0]
        #Fit function
        I_theo = E2moff((x,y),*PSFpopt)
        #filter out noisy pixels at 5sigma level (cos rays/hot pix)
        x, y, intens = PSFclean(x,y,intens,I_theo,skyN,sat,10)

        #calculate better PSF from cleaner data
        fitpopt, fitpcov = curve_fit(lambda (x, y),A,x0,y0: E2moff((x, y),A,ax,ay,b,theta,X0,Y0), (x,y), intens, sigma=np.sqrt(np.absolute(intens)+skyN**2), p0=fitpopt, bounds=bounds, absolute_sigma=True, maxfev=maxfev)
        try:
            #try to calculate fit error
            fitperr = np.sqrt(np.diag(fitpcov))
        except:
            try:
                #take closer initial conditions
                fitpopt, fitpcov = curve_fit(lambda (x, y),A,x0,y0: E2moff((x, y),A,ax,ay,b,theta,x0,y0), (x,y), intens, sigma=np.sqrt(np.absolute(intens)+skyN**2), p0=fitpopt, bounds=bounds, absolute_sigma=True, maxfev=maxfev)
                fitperr = np.sqrt(np.diag(fitpcov))
            except:
                fitperr = [0]*3
        #parameters fitted to source
        A, Aerr = fitpopt[0], fitperr[0]
        X0, X0err = fitpopt[1], fitperr[1]
        Y0, Y0err = fitpopt[2], fitperr[2]
        PSFpopt = [A,ax,ay,b,theta,X0,Y0]
        PSFperr = [Aerr,axerr,ayerr,berr,thetaerr,X0err,Y0err]
        #calculate goodness of fit
        I_theo = E2moff((x, y),*PSFpopt)
        X2dof = np.sum(np.square((intens-I_theo)/np.sqrt(np.absolute(intens)+skyN**2)))/(len(intens)-len(fitpopt))
    except:
        #catastrophic failure of PSF fitting
        print "PSF fitting catastrophic failure"
        PSFpopt = [0]*7
        PSFperr = [0]*7
        X2dof = 0
    if verbosity > 0:
        print "PSF moffat fit parameters"
        print "[A,ax,ay,b,theta,X0,Y0] = "+str(PSFpopt)
        print "parameter errors = "+str(PSFperr)
        print "Chi2 = "+str(X2dof)
        print "[FWHMx', FWHMy]' = "+str([FWHMx, FWHMy])

    #graph fits if verbosity is high enough
    if verbosity > 1:
        PSF_plot(image, x0, y0, PSFpopt, X2dof, skypopt, skyN, fitsky, fsize*fwhm)
    #check if fit is ridiculous
    if E2moff_verify(PSFpopt, x0, y0):
        #not ridiculous, give back fit
        return PSFpopt, PSFperr, X2dof, skypopt, skyN
    else:
        return [0]*7, [0]*7, 0, [0]*3, skyN

#function: scales PSF to source location
def PSFscale(image, PSF, PSFerr, x0, y0, fitsky=True, sat=40000.0, verbosity=0):
    
    from scipy.optimize import curve_fit
    from PSFlib import D2plane, E2moff, E2moff_toFWHM, E2moff_verify

    #get given fit parameters
    ax, axerr = PSF[0], PSFerr[0]
    ay, ayerr = PSF[1], PSFerr[1]
    b, berr = PSF[2], PSFerr[2]
    theta, thetaerr = PSF[3], PSFerr[3]
    FWHMx, FWHMy = E2moff_toFWHM(ax, ay, b)
    fwhm = max(FWHMx, FWHMy)

    #fit sky background in an annulus
    skypopt, skyperr, skyX2dof, skyN = SkyFit(image, [x0], [y0], fwhm, sat, verbosity=0)

    #get fit box to fit psf
    fsize = 3
    intens, x, y = ap_get(image, x0, y0, 0, fsize*fwhm)
    if fitsky:
        #get sky background
        sky = D2plane((x,y),*skypopt)
        #subtract sky background
        intens = intens - sky

    #filter out saturated pixels
    x, y, intens = PSFclean(x,y,intens,intens,skyN,sat,10)
    
    try:
        #fit 2d fixed psf to background subtracted source light
        est = [image[int(y0)][int(x0)]]
        fitpopt, fitpcov = curve_fit(lambda (x, y),A: E2moff((x, y),A,ax,ay,b,theta,x0,y0), (x,y), intens, sigma=np.sqrt(np.absolute(intens)+skyN**2), p0=est, absolute_sigma=True, maxfev=maxfev)
        #parameters fitted to source
        PSFpopt = [fitpopt[0],ax,ay,b,theta,x0,y0]
        #Fit function
        I_theo = E2moff((x,y),*PSFpopt)
        #filter out noisy pixels at 5sigma level (cos rays/hot pix)
        x, y, intens = PSFclean(x,y,intens,I_theo,skyN,sat,10)

        #calculate better PSF from cleaner data
        fitpopt, fitpcov = curve_fit(lambda (x, y),A: E2moff((x, y),A,ax,ay,b,theta,x0,y0), (x,y), intens, sigma=np.sqrt(np.absolute(intens)+skyN**2), p0=fitpopt, absolute_sigma=True, maxfev=maxfev)
        try:
            #try to calculate fit error
            fitperr = np.sqrt(np.diag(fitpcov))
        except:
            try:
                #take closer initial conditions
                fitpopt, fitpcov = curve_fit(lambda (x, y),A: E2moff((x, y),A,ax,ay,b,theta,x0,y0), (x,y), intens, sigma=np.sqrt(np.absolute(intens)+skyN**2), p0=fitpopt, absolute_sigma=True, maxfev=maxfev)
                fitperr = np.sqrt(np.diag(fitpcov))
            except:
                fitperr = [0]
        #parameters fitted to source
        A, Aerr = fitpopt[0], fitperr[0]
        PSFpopt = [A,ax,ay,b,theta,x0,y0]
        PSFperr = [Aerr,axerr,ayerr,berr,thetaerr,0,0]
        #calculate goodness of fit
        I_theo = E2moff((x, y),*PSFpopt)
        X2dof = np.sum(np.square((intens-I_theo)/np.sqrt(np.absolute(intens)+skyN**2)))/(len(intens)-len(fitpopt))
    except:
        #catastrophic failure of PSF fitting
        print "PSF fitting catastrophic failure"
        PSFpopt = [0]*7
        PSFperr = [0]*7
        X2dof = 0
    if verbosity > 0:
        print "PSF moffat fit parameters"
        print "[A,ax,ay,b,theta,X0,Y0] = "+str(PSFpopt)
        print "parameter errors = "+str(PSFperr)
        print "Chi2 = "+str(X2dof)
        print "[FWHMx', FWHMy]' = "+str([FWHMx, FWHMy])

    #graph fits if verbosity is high enough
    if verbosity > 1:
        PSF_plot(image, x0, y0, PSFpopt, X2dof, skypopt, skyN, fitsky, fsize*fwhm)
    #check if fit is ridiculous, give back no fit
    if E2moff_verify(PSFpopt, x0, y0):
        return PSFpopt, PSFperr, X2dof, skypopt, skyN
    else:
        return [0]*7, [0]*7, 0, [0]*3, skyN

#function: fit multiple PSFs
def PSFmulti(image, PSF, PSFerr, psftype, x0, y0, fitsky=True, sat=40000.0, verbosity=0):

    from scipy.optimize import curve_fit
    from PSFlib import D2plane, E2moff_multi, E2moff_toFWHM, E2moff_verify

    maxfev = 100000

    #get given fit parameters
    ax, axerr = PSF[0], PSFerr[0]
    ay, ayerr = PSF[1], PSFerr[1]
    b, berr = PSF[2], PSFerr[2]
    theta, thetaerr = PSF[3], PSFerr[3]
    FWHMx, FWHMy = E2moff_toFWHM(ax, ay, b)
    fwhm = max(FWHMx, FWHMy)

    #number of objects
    Nobj = len(psftype)

    #fit sky background in an annulus
    skypopt, skyperr, skyX2dof, skyN = SkyFit(image, x0, y0, fwhm, sat, verbosity=0)
    
    #get fit box (around all sources) to multi-fit psf
    fsize = 3
    intens, x, y = ap_multi(image, x0, y0, 0, fsize*fwhm)
    if fitsky:
        #get sky background
        sky = D2plane((x,y),*skypopt)
        #subtract sky background
        intens = intens - sky

    #filter out saturated pixels
    x, y, intens = PSFclean(x,y,intens,intens,skyN,sat,10)

    #given parameters for fitting
    given = []
    for i in range(Nobj):
        if psftype[i] == 3:
            #given is empty, general psf params are all in free
            given.append([])
        if psftype[i] == 2:
            #given contains [ax,ay,b,theta], free has [A, x0, y0]
            given.append(PSF)
        if psftype[i] == 1:
            given.append([PSF[0],PSF[1],PSF[2],PSF[3],x0[i],y0[i]])
    given = np.array(given)
    #estimate free parameters for fitting
    est = []
    for i in range(Nobj):
        if psftype[i] == 3:
            #given is empty, general psf params are all in free
            est = np.concatenate((est,[image[int(y0[i])][int(x0[i])],fwhm/4.0,fwhm,3.0,0.0,x0[i],y0[i]]))
        if psftype[i] == 2:
            #given contains [ax,ay,b,theta], free has [A, x0, y0]
            est = np.concatenate((est,[image[int(y0[i])][int(x0[i])],x0[i],y0[i]]))
        if psftype[i] == 1:
            est = np.concatenate((est,[image[int(y0[i])][int(x0[i])]]))
    #estimate bounds on free parameters
    lbounds, ubounds = [], []
    for i in range(Nobj):
        if psftype[i] == 3:
            #given is empty, general psf params are all in free
            lbounds = np.concatenate((lbounds,[-float("Inf"),0.01,0.01,1.01,0.0,0.0,0.0]))
            ubounds = np.concatenate((ubounds,[float("Inf"),4*fwhm,4*fwhm,float("Inf"),179.99,image.shape[0],image.shape[1]]))
        if psftype[i] == 2:
            #given contains [ax,ay,b,theta], free has [A, x0, y0]
            lbounds = np.concatenate((lbounds,[-float("Inf"),0.0,0.0]))
            ubounds = np.concatenate((ubounds,[float("Inf"),image.shape[0],image.shape[1]]))
        if psftype[i] == 1:
            lbounds = np.concatenate((lbounds,[-float("Inf")]))
            ubounds = np.concatenate((ubounds,[float("Inf")]))
    bounds = (lbounds,ubounds)
    
    try:
        #fit 2d fixed psf to background subtracted source light
        fitpopt, fitpcov = curve_fit(lambda (x, y),*free: E2moff_multi((x, y),psftype, given, free), (x,y), intens, sigma=np.sqrt(np.absolute(intens)+skyN**2), p0=est, bounds=bounds, absolute_sigma=True, maxfev=maxfev)
        #Fit function
        I_theo = E2moff_multi((x, y),psftype, given, fitpopt)
        #filter out noisy pixels at 5sigma level (cos rays/hot pix)
        x, y, intens = PSFclean(x,y,intens,I_theo,skyN,sat,10)

        #calculate better PSF from cleaner data
        fitpopt, fitpcov = curve_fit(lambda (x, y),*free: E2moff_multi((x, y),psftype, given, free), (x,y), intens, sigma=np.sqrt(np.absolute(intens)+skyN**2), p0=fitpopt, bounds=bounds, absolute_sigma=True, maxfev=maxfev)
        try:
            #try to calculate fit error
            fitperr = np.sqrt(np.diag(fitpcov))
        except:
            try:
                #take closer initial conditions
                fitpopt, fitpcov = curve_fit(lambda (x, y),*free: E2moff_multi((x, y),psftype, given, free), (x,y), intens, sigma=np.sqrt(np.absolute(intens)+skyN**2), p0=fitpopt, bounds=bounds, absolute_sigma=True, maxfev=maxfev)
                fitperr = np.sqrt(np.diag(fitpcov))
            except:
                fitperr = [0]*len(est)
        #parameters fitted to source
        PSFpopt, PSFperr = [], []
        count = 0
        for i in range(Nobj):
            if psftype[i] == 3:
                #given is empty, general psf params are all in free
                PSFpopt.append(fitpopt[count:count+7])
                PSFperr.append(fitperr[count:count+7])
                count = count+7
            if psftype[i] == 2:
                #given contains [ax,ay,b,theta], free has [A, x0, y0]
                PSFpopt.append([fitpopt[count],given[i][0],given[i][1],given[i][2],given[i][3],fitpopt[count+1],fitpopt[count+2]])
                PSFperr.append([fitperr[count],axerr,ayerr,berr,thetaerr,fitperr[count+1],fitperr[count+2]])
                count = count+3
            if psftype[i] == 1:
                #given contains [ax,ay,b,theta,x0,y0], free has [A]
                PSFpopt.append([fitpopt[count],given[i][0],given[i][1],given[i][2],given[i][3],given[i][4],given[i][5]])
                PSFperr.append([fitperr[count],axerr,ayerr,berr,thetaerr,0,0])
                count = count+1
        #calculate goodness of fit
        I_theo = E2moff_multi((x, y),psftype, given, fitpopt)
        X2dof = np.sum(np.square((intens-I_theo)/np.sqrt(np.absolute(intens)+skyN**2)))/(len(intens)-len(fitpopt))
        
        #Graph residual if verbosity is high enough
        if verbosity > 1:
            import matplotlib.pyplot as plt
            I_t = np.copy(image)
            x1, x2 = int(min(x0)-3*fwhm),int(max(x0)+3*fwhm)
            y1, y2 = int(min(y0)-3*fwhm),int(max(y0)+3*fwhm)
            for i in np.arange(x1, x2):
                for j in np.arange(y1, y2):
                    I_t[j][i] = E2moff_multi((i, j),psftype, given, fitpopt) + D2plane((i, j),*skypopt)
            sb = image[y1:y2,x1:x2]-I_t[y1:y2,x1:x2]
            sbmax = 5*skyN
            sbmin = -5*skyN
            plt.title("Multi-object fit residual")
            plt.imshow(sb, cmap='Greys',vmax=sbmax,vmin=sbmin)
            plt.colorbar()
            plt.scatter(np.array(x0, dtype=int)-x1, np.array(y0, dtype=int)-y1, color='r', marker='.')
            plt.show()
    except:
        #catastrophic failure of PSF fitting
        print "PSF fitting catastrophic failure"
        PSFpopt = [[0]*7]*Nobj
        PSFperr = [[0]*7]*Nobj
        X2dof = 0
    if verbosity > 0:
        print "PSF moffat fit parameters"
        for i in range(Nobj):
            print "Object "+str(i+1)+":"
            print "[A,ax,ay,b,theta,X0,Y0] = "+str(PSFpopt[i])
            print "parameter errors = "+str(PSFperr[i])
            print "Chi2 = "+str(X2dof)
            FWHMx, FWHMy = E2moff_toFWHM(PSFpopt[i][1], PSFpopt[i][2], PSFpopt[i][3])
            print "[FWHMx', FWHMy]' = "+str([FWHMx, FWHMy])

    #graph fits if verbosity is high enough
    if verbosity > 1:
        PSFmulti_plot(image, x0, y0, PSFpopt, X2dof, skypopt, skyN, fitsky, window=15)
        
    #check if fit is ridiculous, give back no fit
    ridic = False
    checks = []
    for i in range(Nobj):
        if psftype[i] != 3:
            if not E2moff_verify(PSFpopt[i], x0[i], y0[i]):
                ridic = True
                if verbosity > 0:
                    print "Bad PSF: Object "+str(i+1)
                    print dist(PSFpopt[i][5],PSFpopt[i][6],x0[i],y0[i])
    if not ridic:
        #None of the fits were ridiculous
        return PSFpopt, PSFperr, X2dof, skypopt, skyN
    else:
        return [[0]*7 for i in range(Nobj)], [[0]*7 for i in range(Nobj)], 0, [0]*3, skyN

#function: plot PSF fitting
def PSF_plot(image, x0, y0, PSFpopt, X2dof, skypopt, skyN, fitsky, window=15):

    import matplotlib.pyplot as plt
    from PSFlib import D2plane, E2moff, E2moff_toFWHM
    
    #get fit parameters
    A = PSFpopt[0]
    ax = abs(PSFpopt[1])
    ay = abs(PSFpopt[2])
    b = PSFpopt[3]
    theta = PSFpopt[4]
    X0 = PSFpopt[5]
    Y0 = PSFpopt[6]
    FWHMx, FWHMy = E2moff_toFWHM(ax, ay, b)
    
    x = np.arange(x0-window,x0+window+1,dtype=int)
    xt = np.arange(x0-window,x0+window+1,0.1)
    Ix_im = np.array([image[int(y0)][i] for i in x])
    y = np.arange(y0-window,y0+window+1,dtype=int)
    yt = np.arange(y0-window,y0+window+1,0.1)
    Iy_im = np.array([image[i][int(x0)] for i in y])
    if FWHMx*FWHMy != 0:
        #compute PSF fit
        Ix_theo = E2moff((xt,np.array([int(y0)]*len(xt))),*PSFpopt)
        Ix_res = Ix_im - E2moff((x,np.array([int(y0)]*len(x))),*PSFpopt)
        Iy_theo = E2moff((np.array([int(x0)]*len(yt)),yt),*PSFpopt)
        Iy_res = Iy_im - E2moff((np.array([int(x0)]*len(y)),y),*PSFpopt)
    else:
        Ix_theo = 0
        Iy_theo = 0
        Ix_res = Ix_im
        Iy_res = Iy_im
    
    if fitsky:
        Ix_theo = Ix_theo+D2plane((xt,np.array([int(y0)]*len(xt))),*skypopt)
        Ix_res = Ix_res-D2plane((x,np.array([int(y0)]*len(x))),*skypopt)
        Iy_theo = Iy_theo+D2plane((np.array([int(x0)]*len(yt)),yt),*skypopt)
        Iy_res = Iy_res-D2plane((np.array([int(x0)]*len(y)),y),*skypopt)
        
    #plot psf moffat in X and Y slice side by side
    f, ax = plt.subplots(2, 2, sharex="col", sharey="row")
    ax[0][0].errorbar(x,Ix_im,yerr=np.sqrt(np.absolute(Ix_im)+skyN**2),fmt='r+', label='slice at Y='+str(int(y0)))
    ax[0][0].plot(xt, Ix_theo, label='fit X2/dof='+str(X2dof)[:6])
    ax[0][0].set_xlabel('x-pixel')
    ax[0][0].set_ylabel('data #')
    ax[0][0].set_xlim(min(x),max(x))
    ax[0][0].legend()
    ax[1][0].errorbar(x,Ix_res,yerr=np.sqrt(np.absolute(Ix_im)+skyN**2),fmt='r+', label='residuals')
    ax[1][0].plot(x,np.zeros(len(x)))
    ax[1][0].legend()
    ax[0][1].errorbar(y,Iy_im,yerr=np.sqrt(np.absolute(Iy_im)+skyN**2),fmt='r+', label='slice at X='+str(int(x0)))
    ax[0][1].plot(yt, Iy_theo, label='fit X2/dof='+str(X2dof)[:6])
    ax[0][1].set_xlabel('y-pixel')
    ax[0][1].set_ylabel('data #')
    ax[0][1].set_xlim(min(y),max(y))
    ax[0][1].legend()
    ax[1][1].errorbar(y,Iy_res,yerr=np.sqrt(np.absolute(Iy_im)+skyN**2),fmt='r+', label='residuals')
    ax[1][1].plot(y,np.zeros(len(y)))
    ax[1][1].legend()
    f.subplots_adjust(wspace=0)
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels()[0] for a in ax[1, :]], visible=False)
    plt.setp([a.get_yticklabels()[0] for a in ax[:, 0]], visible=False)
    plt.setp([a.get_xticklabels()[-1] for a in ax[1, :]], visible=False)
    plt.setp([a.get_yticklabels()[-1] for a in ax[:, 0]], visible=False)
    plt.suptitle("PSF Moffat Fit")
    plt.show()

#function: plot multi-object PSF fitting
def PSFmulti_plot(image, x0, y0, PSFpopt, X2dof, skypopt, skyN, fitsky, window=15):

    import matplotlib.pyplot as plt
    from PSFlib import D2plane, E2moff_multi, E2moff_toFWHM

    psfpopt = np.concatenate(PSFpopt)
    psftype = [3]*len(PSFpopt)

    for i in range(len(PSFpopt)):
        ax = abs(PSFpopt[i][1])
        ay = abs(PSFpopt[i][2])
        b = PSFpopt[i][3]
        FWHMx, FWHMy = E2moff_toFWHM(ax, ay, b)
        
        x = np.arange(x0[i]-window,x0[i]+window+1,dtype=int)
        xt = np.arange(x0[i]-window,x0[i]+window+1,0.1)
        Ix_im = np.array([image[int(y0[i])][j] for j in x])
        y = np.arange(y0[i]-window,y0[i]+window+1,dtype=int)
        yt = np.arange(y0[i]-window,y0[i]+window+1,0.1)
        Iy_im = np.array([image[j][int(x0[i])] for j in y])
        if FWHMx*FWHMy != 0:
            #compute PSF fit
            Ix_theo = E2moff_multi((xt,np.array([int(y0[i])]*len(xt))), psftype, [], psfpopt)
            Ix_res = Ix_im - E2moff_multi((x,np.array([int(y0[i])]*len(x))), psftype, [], psfpopt)
            Iy_theo = E2moff_multi((np.array([int(x0[i])]*len(yt)),yt), psftype, [], psfpopt)
            Iy_res = Iy_im - E2moff_multi((np.array([int(x0[i])]*len(y)),y), psftype, [], psfpopt)
        else:
            Ix_theo = 0
            Iy_theo = 0
            Ix_res = Ix_im
            Iy_res = Iy_im
    
        if fitsky:
            Ix_theo = Ix_theo+D2plane((xt,np.array([int(y0[i])]*len(xt))),*skypopt)
            Ix_res = Ix_res-D2plane((x,np.array([int(y0[i])]*len(x))),*skypopt)
            Iy_theo = Iy_theo+D2plane((np.array([int(x0[i])]*len(yt)),yt),*skypopt)
            Iy_res = Iy_res-D2plane((np.array([int(x0[i])]*len(y)),y),*skypopt)
        
        #plot psf moffat in X and Y slice side by side
        f, ax = plt.subplots(2, 2, sharex="col", sharey="row")
        ax[0][0].errorbar(x,Ix_im,yerr=np.sqrt(np.absolute(Ix_im)+skyN**2),fmt='r+', label='slice at Y='+str(int(y0[i])))
        ax[0][0].plot(xt, Ix_theo, label='fit X2/dof='+str(X2dof)[:6])
        ax[0][0].set_xlabel('x-pixel')
        ax[0][0].set_ylabel('data #')
        ax[0][0].set_xlim(min(x),max(x))
        ax[0][0].legend()
        ax[1][0].errorbar(x,Ix_res,yerr=np.sqrt(np.absolute(Ix_im)+skyN**2),fmt='r+', label='residuals')
        ax[1][0].plot(x,np.zeros(len(x)))
        ax[1][0].legend()
        ax[0][1].errorbar(y,Iy_im,yerr=np.sqrt(np.absolute(Iy_im)+skyN**2),fmt='r+', label='slice at X='+str(int(x0[i])))
        ax[0][1].plot(yt, Iy_theo, label='fit X2/dof='+str(X2dof)[:6])
        ax[0][1].set_xlabel('y-pixel')
        ax[0][1].set_ylabel('data #')
        ax[0][1].set_xlim(min(y),max(y))
        ax[0][1].legend()
        ax[1][1].errorbar(y,Iy_res,yerr=np.sqrt(np.absolute(Iy_im)+skyN**2),fmt='r+', label='residuals')
        ax[1][1].plot(y,np.zeros(len(y)))
        ax[1][1].legend()
        f.subplots_adjust(wspace=0)
        f.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels()[0] for a in ax[1, :]], visible=False)
        plt.setp([a.get_yticklabels()[0] for a in ax[:, 0]], visible=False)
        plt.setp([a.get_xticklabels()[-1] for a in ax[1, :]], visible=False)
        plt.setp([a.get_yticklabels()[-1] for a in ax[:, 0]], visible=False)
        plt.suptitle("PSF Moffat Fit")
        plt.show()
    
#function: calculates photometric intensity and sky background using PSF
def PSF_photometry(image, x0, y0, PSFpopt, PSFperr, skypopt, skyN, verbosity=0):

    from PSFlib import D2plane, E2moff, E2moff_verify, E2moff_integrate, E2moff_apsize
    
    #get some critical values
    A, Aerr = PSFpopt[0], PSFperr[0]
    ax, axerr = abs(PSFpopt[1]), PSFperr[1]
    ay, ayerr = abs(PSFpopt[2]), PSFperr[2]
    b, berr = PSFpopt[3], PSFperr[3]
    theta, thetaerr = PSFpopt[4], PSFperr[4]
    X0 = PSFpopt[5]
    Y0 = PSFpopt[6]
    
    #fraction of light to include (Kron)
    frac = 0.9
    if b > 1:
        #integrate source PSF
        Io = E2moff_integrate(A,ax,ay,b,frac)
        #compute optimal aperture size (90% source light)
        ap_size = E2moff_apsize(ax,ay,b,frac)
        sigmao = np.sqrt(np.absolute(Io) + (skyN**2)*ap_size)
        SNo = Io/sigmao
    else:
        if verbosity > 0:
            print "Divergent integral, PSF not viable for PSF photometry"
        Io = 0
        SNo = 0
        SNr = 0
        opt_r = 0
    #output information
    if verbosity > 0:
        print "\n"
        print "signal to noise: "+str(SNo)
        print "\n"
    
    #return intensity, and signal to noise
    return Io, SNo

#function: calculates photometric intensity and sky background using aperture
def Ap_photometry(image, x0, y0, skypopt, skyN, radius=None, PSF=None, fitsky=True, verbosity=0):
    
    from PSFlib import D2plane, E2moff_verify
    
    if radius is None and E2moff_verify(PSF, x0, y0):
        #get some critical values to calculate radius
        frac = 0.9 #Kron aperture light fraction
        ax = abs(PSF[0])
        ay = abs(PSF[1])
        b = PSF[2]
        theta = PSF[3]
        radius = 0.5*(ax+ay)*np.sqrt(np.power(1 - frac,1/(1-b)) - 1)
    else:
        radius = 0

    if radius != 0:
        #extract PSF aperture
        PSF_extract, x, y = ap_get(image, x0, y0, 0, radius)
        if fitsky:
            #subtract sky background
            PSF_sky = D2plane((x,y),*skypopt)
            PSF_extract = PSF_extract - PSF_sky
        #integrate aperture
        Io = np.sum(PSF_extract)
        #proper signal to noise calculation for unscaled intensities
        #noise is sqrt(intensity) is the best we can do
        sigmar = np.sqrt(np.absolute(Io) + (skyN**2)*PSF_extract.size)
        SNo = Io/sigmar
    else:
        print "Unable to integrate, aperture=0."
        Io = 0
        SNo = 0
    
    if verbosity > 0:
        print "\n"
        print "Signal to noise: "+str(SNo)
        print "Pixel aperture radius: "+str(radius)
        print "\n"

    #graph photometry if verbosity is high enough
    if verbosity > 1 and radius > 0:

        import matplotlib.pyplot as plt

        #Extract PSF
        x = np.arange(x0-radius,x0+radius+1,dtype=int)
        Ix_im = np.array([image[int(y0)][i] for i in x])  
        y = np.arange(y0-radius,y0+radius+1,dtype=int)
        Iy_im = np.array([image[i][int(x0)] for i in y])
        if fitsky:
            Ix_im = Ix_im-D2plane((x,np.array([int(y0)]*len(x))),*skypopt)
            Iy_im = Iy_im-D2plane((np.array([int(x0)]*len(y)),y),*skypopt)
        
        #plot psf moffat in X and Y slice side by side
        f, ax = plt.subplots(1, 2, sharey="row")
        ax[0].errorbar(x,Ix_im,yerr=np.sqrt(np.absolute(Ix_im)+skyN**2),fmt='r+', label='slice at Y='+str(int(y0)))
        ax[0].set_xlabel('x-pixel')
        ax[0].set_ylabel('data #')
        ax[0].set_xlim(min(x),max(x))
        ax[0].legend()
        ax[1].errorbar(y,Iy_im,yerr=np.sqrt(np.absolute(Iy_im)+skyN**2),fmt='r+', label='slice at X='+str(int(x0)))
        ax[1].set_xlabel('y-pixel')
        ax[1].set_ylabel('data #')
        ax[1].set_xlim(min(y),max(y))
        ax[1].legend()
        f.subplots_adjust(wspace=0)
        plt.suptitle("Aperture Section")
        plt.show()
        plt.setp([a.get_yticklabels()[0] for a in ax], visible=False)
        plt.setp([a.get_yticklabels()[-1] for a in ax], visible=False)
    
        #plot psf aperture to snr, intensity
        ap_r = np.linspace(1,2*radius,100)
        Is = np.zeros(len(ap_r))
        sigmas = np.zeros(len(ap_r))
        for i, rad in enumerate(ap_r):
            #extract PSF aperture
            PSF_extract, x, y = ap_get(image, x0, y0, 0, rad)
            if fitsky:
                #subtract sky background
                PSF_sky = D2plane((x,y),*skypopt)
                PSF_extract = PSF_extract - PSF_sky
            #integrate aperture
            Is[i] = np.sum(PSF_extract)
            sigmas[i] = np.sqrt(np.absolute(Is[i]) + (skyN**2)*PSF_extract.size)
        SNs = Is/sigmas

        f, ax = plt.subplots(2, sharex=True)
        ax[0].set_title("Aperture Signal to Noise")
        ax[0].plot(ap_r,SNs,label='SNR')
        ax[0].plot([radius, radius], [min(SNs), max(SNs)], 'g')
        ax[0].set_ylabel("SN")
        ax[0].legend()
        ax[1].errorbar(ap_r,Is,yerr=sigmas,fmt='r+',label='summed intens')
        ax[1].plot([radius, radius], [min(Is), max(Is)], 'g')
        ax[1].set_ylabel("I")
        ax[1].set_xlabel("Aperture R (pix)")
        ax[1].legend()
        plt.setp([a.get_yticklabels()[0] for a in ax], visible=False)
        plt.setp([a.get_yticklabels()[-1] for a in ax], visible=False)
        f.subplots_adjust(hspace=0)
        plt.show()
    elif verbosity > 1:
        print "Unable to plot, aperture=0."

    #return intensity, and signal to noise
    return Io, SNo
