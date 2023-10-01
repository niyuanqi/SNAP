#################################################################
# Name:     FitsSandbox.py                                      #
# Author:   Yuan Qi Ni                                          #
# Version:  November 15, 2016                                   #
# Function: Program contains routines for generating various    #
#           test files in fits format.                          #
#################################################################

#essential modules
import numpy as np

#class: exception, data can not be handled
class DataError(Exception):
    """
    Exception raised for data that FitsSandbox cannot handle
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return "Data incompatible: must be "+repr(self.value) 

#function: plot fits images
def plotFits(image, rmax=1., rmin=0., title=None, invert=True, cmap='gray'):
    """
    Plot images in HDUList
    Input------------------------------
    image  : image (ndarray, ImageHDU or HDUlist)
    rmax   : ratio of colorbar vmax wrt max pixel
    rmin   : ratio of colorbar vmin wrt max pixel
    title  : string title of image
    invert : boolean True to invert colormap when plotting
    cmap   : matplotlib colormap
    """

    from astropy.io import fits
    import matplotlib.pyplot as plt
    
    if type(image) is fits.HDUList:
        #if HDUlist, break up into ImageHDUs
        for i, hdu in enumerate(image):
            if title is None:
                plotFits(hdu, rmax, rmin, title, cmap, invert)
            else:
                plotFits(hdu, rmax, rmin, title+" "+str(i), cmap, invert)
    elif type(image) is fits.ImageHDU or type(image) is fits.PrimaryHDU:
        #if ImageHDU, plot data image
        plotFits(image.data, rmax, rmin, title, cmap, invert)
    elif type(image) is np.ndarray:
        #if data image, plot
        vmax = rmax*np.amax(image)
        vmin = rmin*np.amin(image)
        if invert:
            cmap = cmap+'_r'
        if title is not None:
            plt.title(title)
        plt.imshow(image, interpolation='nearest',
                   vmin=vmin, vmax=vmax, cmap=cmap, origin='lower')
        plt.colorbar()
        plt.xlabel("x pixel")
        plt.ylabel("y pixel")
        plt.show()
    else:
        #not an image hdu object
        raise DataError([fits.HDUList,fits.PrimaryHDU,
                         fits.ImageHDU,np.ndarray])

#function: plot fits images
def genFits(filename):
    """
    Generate a blank HDUlist
    Input------------------------------
    filename  : name of HDUlist
    """
    print "gen"



#function: replace fits hdu image with noise
def map_Noise(hdulist, val, scale):
    """
    Map noise over HDUlist images for testing
    Input------------------------------
    hdulist : Fits Header Data Unit (ImageHDU or HDUlist)
    val     : noise central mean
    scale   : noise standard deviation
    Output-----------------------------
    out     : Image HDU or HDUlist with noise mapped
    """

    from astropy.io import fits

    if type(hdulist) is fits.HDUList:
        #if HDUlist, break up into ImageHDUs
        out = []
        for i, hdu in enumerate(hdulist):
            #recursion to map ImageHDU to noise
            out.append(map_Noise(hdu, val, scale))
        out = fits.HDUList(out)
    elif type(hdulist) is fits.ImageHDU or type(hdulist) is fits.PrimaryHDU:
        #if ImageHDU, map data to noise
        noise = map_Noise(hdulist.data, val, scale)
        #generate new hdu with noise
        out = hdulist.copy()
        out.data = noise
    elif type(hdulist) is np.ndarray:
        template = np.copy(hdulist)
        if scale != 0.0:
            noise = np.random.normal(val, scale, template.size)
            noise = np.reshape(noise, template.shape)
        else:
            noise = np.ones(template.size)*val
            noise = np.reshape(noise, template.shape)
        out = noise
    else:
        #not an image hdu object
        raise DataError([fits.HDUList,fits.PrimaryHDU,
                         fits.ImageHDU,np.ndarray])
    #return noise mapped hdu
    return out

