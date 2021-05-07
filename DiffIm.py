#################################################################
# Name:     DiffIm.py                                           #
# Author:   Yuan Qi Ni                                          #
# Version:  January 10, 2016                                    #
# Function: Program contains routines used for fits image       #
#           subtraction using wcsremap astrometric alignment    #
#           and hotpants photometric alignment software by      #
#           Andrew Becker. Based on make_image_diff.py          #
#################################################################

#Sample usage
#python make_image_diff_n.py N300-1.Q0.B.161030_0504.C.034140.005604N3646.0060.nh.fits N300-1.Q0.B.150626_1842-151019_1231.XCSA.005605N3646.00081.00081.FM30.BS0512.ALL.coadd.NEW.REF.fits diff.fits conv.fits

#remove incompatible header information
def remove_tan_from_header(inname, outdir, extnum=0):
    """
    update header for astrometry.net produced header file
    """
    
    from astropy.io import fits
    import os

    basename = os.path.basename(inname)
    outname = os.path.join(outdir, basename)

    f = fits.open(inname)
    header = f[extnum].header
    for ctype in ["CTYPE1", "CTYPE2"]:
        if  header[ctype].endswith("-SIP"):
            header[ctype] = header[ctype][:-4]

    f.writeto(outname, clobber=True, output_verify="fix")
    return outname

#use wcsremap to match images astrometrically
def run_wcsremap(ref_name, src_name, outdir):

    import subprocess
    import os

    basename = os.path.basename(ref_name)
    outname = os.path.join(outdir, basename)
    
    #call wcsremap
    subprocess.call(['wcsremap', '-template', src_name, '-source',
                     ref_name, '-outIm', outname])

    return outname

#make difference image
def basic_diff_image(src_name, ref_name, out_name, conv_name, tmpdir="DITemp", delete_temp=True):
    try:
        
        import os
        import subprocess
        
        #make temporary directory
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)
            
        #fix header info
        src_name2 = remove_tan_from_header(src_name, tmpdir)
        
        #remap reference file to source file coordinates
        remapped_ref = run_wcsremap(ref_name, src_name2, tmpdir)
        
        #hotpants arguments
        default_flags = ['hotpants', '-inim', src_name, '-tmplim',
                         remapped_ref, '-outim', out_name, '-oci',
                         conv_name, '-tl', '-100', '-il', '-100',
                         '-n', 'i', '-c', 't', '-ko', '2',
                         '-nrx', '2', '-nry', '2',
                         '-nsx', '15', '-nsy', '15',
                         '-ng','4','7','0.70','6','1.50','4','3.00','3','6.0']
            
        #subtract remapped reference file from source file
        subprocess.call(default_flags)

        print "SUBTRACTION COMPLETE"
        print "output:",out_name
        
    finally:
        if delete_temp:
            
            import shutil
            
            shutil.rmtree(tmpdir)

#make difference image
def make_diff_image(src_name, ref_name, out_name, conv_name, tmp_fwhm=None, src_fwhm=None, imx=9216, imy=9232, tmpdir="DITemp", delete_temp=True):
    try:
        
        import os
        import subprocess
        import numpy as np

        #figure out band
        fo = src_name.split('/')[2]
        fo = '.'.join(fo.split('.')[2:5])
        band = fo[0]
        print ""
        print "Band =", band

        #check which image is better
        if src_fwhm is not None and tmp_fwhm is not None:
            if src_fwhm < tmp_fwhm:
                print "Image is better than template."
                sigma_match = np.sqrt(tmp_fwhm**2 - src_fwhm**2)/2.
                fwhm_c = tmp_fwhm
                fwhm_f = src_fwhm
                better = 'i'
            else:
                print "Template is better than image."
                sigma_match = np.sqrt(src_fwhm**2 - tmp_fwhm**2)/2.
                fwhm_c = src_fwhm
                fwhm_f = tmp_fwhm
                better = 't'
        else:
            print "Assuming template is better."
            sigma_match = None
            fwhm_c = None
            fwhm_f = None
            better = 't'
        print ""

        #make temporary directory
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)

        #fix header info
        src_name2 = remove_tan_from_header(src_name, tmpdir)
        
        #remap reference file to source file coordinates
        remapped_ref = run_wcsremap(ref_name, src_name2, tmpdir)

        #set hotpants flags depending on which is better
        flags = ['-inim', src_name, '-tmplim',
                 remapped_ref, '-outim', out_name, '-oci',
                 conv_name, '-tl', '-100', '-il', '-100',
                 '-n', 'i', '-ko', '2', '-c', better]
        
        #set kernel extraction depending on fwhm
        if (src_fwhm > 0.95*tmp_fwhm) and fwhm_c is not None:
            #parameters that work well for poor science images
            flags += ['-r', str(2.5*fwhm_c/2.0)]
            flags += ['-rss', str(3.0*fwhm_c)]
            flags += ['-nsx', str(round(imx/(30.0*fwhm_c)))]
            flags += ['-nsy', str(round(imy/(30.0*fwhm_c)))]
            flags += ['-ng','3','6',str(fwhm_f/4.0),'4',str(fwhm_f/2.0),'2',str(fwhm_f)]    
        else:
            #parameters that work well for excellent science images
            flags += ['-nsx', '30', '-nsy', '30']
            flags += ['-ng','4','7','0.70','6','1.50','4','3.00','3','6.0']
    
        #call hotpants
        subprocess.call(['hotpants'] + flags)

        print "SUBTRACTION COMPLETE"
        print "output:",out_name
        
    finally:
        if delete_temp:
            
            import shutil
            
            shutil.rmtree(tmpdir)
    
#main function
if __name__ == "__main__":
    
    import argparse
    
    #receive arguments
    parser = argparse.ArgumentParser(description='make diff images.')
    parser.add_argument('src_name', type=str, help='input file name')
    parser.add_argument('ref_name', type=str, help='reference file name')
    parser.add_argument('out_name', type=str, help='output file name')
    parser.add_argument('conv_name', type=str, help='convolved file name')
    parser.add_argument('--debug', dest='debug', action='store_const',
                        const=True, default=False,
                        help='give this flag to not delete temp files')
    args = parser.parse_args()
    #whether to delete temporary directory for intermediate files created
    delete_temp = not args.debug

    #create difference image
    basic_diff_image(args.src_name, args.ref_name,
                     args.out_name, args.conv_name,
                     delete_temp=delete_temp)
