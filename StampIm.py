#################################################################
# Name:     StampIm.py                                          #
# Author:   Yuan Qi Ni                                          #
# Date:     August 3, 2017                                      #
# Function: Program makes png stamp images from fits files.     #
#################################################################

#function: make fits file into stamp image array
def make_stamp_image(filename, ra, dec, radius=100):
    """
    Makes stamp image png from fits file
    ------------------------------------
    Inputs:
    ------------------------------------
    filename = string fits file to be made into stamp
    ra, dec = float degree position of center
    radius = float arcsec radius of image
    ------------------------------------
    Outputs:
    ------------------------------------
    stamp image array
    """

    import numpy as np
    from astropy.io import fits
    from astropy.wcs import WCS
    import matplotlib.pyplot as plt

    #load HDU image
    hdulist = fits.open(filename)
    wcs = WCS(filename)
    image = hdulist[0].data
    #close HDU image
    hdulist.close()

    cX, cY = wcs.all_world2pix(ra, dec, 0)
    cX, cY = int(round(cX)), int(round(cY))
    image = image.T[cX-radius:cX+radius+1].T
    image = image[cY-radius:cY+radius+1]
    return image

#function: create A4 collage of stamp images with filename as subtext
def make_image_collage(files, names, outname, ra, dec, radius=100, scale=0.001, ppi=100, spacing=50 ):
    """
    Makes stamp image pdf from fits file
    ------------------------------------
    Inputs:
    ------------------------------------
    files = stamp images to be incorporated into collage
    outname = string output name of collage pdf
    ra, dec = float degree position of center
    radius = float arcsec radius of image
    scale = float fraction of max intensity to set image scale 
    ppi = int pixels per inch of A4 paper
    spacing = float pixel minimum spacing between stamps
    ------------------------------------
    Outputs:
    ------------------------------------
    pdf file with outname
    """

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    #size of A4 page is 8.5x11 inches
    #calculate pixel dimensions
    width = int(8.5*ppi)
    length = int(11*ppi)

    #a4 sheets
    papers = []
    texts = []
    cents = []

    #make blank array of A4 paper
    paper = np.zeros([length, width])
    marker = [length-1,0]
    textloc = []
    textname = []
    centloc = []

    with PdfPages(outname) as pdf:
        #for each image, place onto A4 paper
        for i, filename in enumerate(files):
            print "Stamping",filename
            #read image
            image = make_stamp_image(filename, ra, dec, radius)
            print "Image size"
            print image.shape
            #leave some space to place image from marker
            print "Deciding marker placement"
            corner1 = [marker[0]-spacing-image.shape[0], marker[1]+spacing]
            corner2 = [corner1[0]+image.shape[0], corner1[1]+image.shape[1]]
            print corner1, corner2
            if corner2[1] > width:
                print "New row"
                #go to new row
                marker = [corner1[0], 0]
                corner1 = [marker[0]-spacing-image.shape[0], marker[1]+spacing]
                corner2 = [corner1[0]+image.shape[0], corner1[1]+image.shape[1]]
                print corner1, corner2
            if corner1[0] < 0:
                print "New page"
                #save current page
                papers.append(paper)
                texts.append([textloc, textname])
                cents.append(centloc)
                #make new page
                paper = np.zeros([length, width])
                marker = [length-1,0]
                textloc = []
                textname = []
                centloc = []
                corner1 = [marker[0]-spacing-image.shape[0], marker[1]+spacing]
                corner2 = [corner1[0]+image.shape[0], corner1[1]+image.shape[1]]
                print corner1, corner2
            print "Image stamped"
            paper[corner1[0]:corner2[0], corner1[1]:corner2[1]] = image
            textloc.append([corner1[1], corner1[0]+image.shape[0]+5])
            centloc.append([corner1[1]+image.shape[0]/2.0, corner1[0]+image.shape[0]/2.0])
            textname.append(names[i])
            marker = [corner2[0]+spacing, corner2[1]]
        #all done? save current page
        papers.append(paper)
        texts.append([textloc, textname])
        cents.append(centloc)

        vmax = scale*max([np.amax(paper) for paper in papers])
        vmin = scale*min([np.amin(paper) for paper in papers])
        for i in range(len(papers)):
            print "plotting page",str(i+1)
            fig = plt.figure(figsize=(8.5, 11.0))
            plt.imshow(papers[i], interpolation='nearest', vmin=vmin, vmax=vmax, cmap='Greys', origin='lower', extent=[0, width-1, 0, length-1])
            textloc = texts[i][0]
            textname = texts[i][1]
            centloc = cents[i]
            for j in range(len(textloc)):
                plt.text(textloc[j][0], textloc[j][1], textname[j])
                cX, cY = centloc[j][0], centloc[j][1]
                plt.plot([cX-0.8*spacing,cX-0.3*spacing],[cY,cY],c='k',lw=3)
                plt.plot([cX+0.3*spacing,cX+0.8*spacing],[cY,cY],c='k',lw=3)
                plt.plot([cX,cX],[cY-0.8*spacing,cY-0.3*spacing],c='k',lw=3)
                plt.plot([cX,cX],[cY+0.3*spacing,cY+0.8*spacing],c='k',lw=3)
            plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)
            pdf.savefig(fig)
            plt.close()
