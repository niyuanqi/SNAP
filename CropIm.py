#################################################################
# Name:     CropIm.py                                           #
# Author:   Yuan Qi Ni                                          #
# Date:     July 25, 2017                                       #
# Function: Program crops fits file to small, source centered   #
#           square with dimension D.                            #
#################################################################

#function: crop image given
def make_crop_image(filename, outname, radius):

    from astropy.io import fits
    from astropy.wcs import WCS
    import numpy as np
    import os

    #load HDU image
    print "loading hdu"
    hdulist = fits.open("../raw/"+filename)
    info = hdulist.info()
    image = hdulist[0].data
    header = hdulist[0].header
    wcs = WCS(header)
    
    #convert position of source to world coordinates
    X, Y = wcs.all_world2pix(ra, dec, 0)
    print "Source located at: " + str(X) + ", " + str(Y)
    
    crop = fits.PrimaryHDU()
    xlims = [int(max(X-radius,0)), int(min(X+radius,image.shape[0]))]
    ylims = [int(max(Y-radius,0)), int(min(Y+radius,image.shape[1]))]
    #print xlims[0], xlims[1], ylims[0], ylims[1]
    crop.data = image[ylims[0]:ylims[1], xlims[0]:xlims[1]]
    crop.header = header
    crop.header.update(wcs[ylims[0]:ylims[1], xlims[0]:xlims[1],].to_header())
    
    print "writing hdu"
    print "cropped "+outname
    if os.path.exists(outname) : os.remove(outname)
    crop.writeto(outname)

#main function
if __name__ == "__main__":
    
    import argparse
    
    #receive arguments
    parser = argparse.ArgumentParser(description='make diff images.')
    parser.add_argument('src_name', type=str, help='input file name')
    parser.add_argument('out_name', type=str, help='output file name')
    parser.add_argument('-r', '--radius', type=float, help='output file size')
    args = parser.parse_args()

    #create crop image
    make_crop_image(args.src_name, args.out_name, args.radius)
