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

#use hotpants to match images photometrically and subtract
def run_hotpants(src_name, tmp_name, out_name, conv_name, fwhm=None, imsize=3692.8, tmp_sat=40000, src_sat=45000, tmp_neg=-100, src_neg=-100, tmp_mask=None, src_mask=None):

    import numpy as np
    import subprocess

    #default_flags = ['hotpants', '-inim', src_name, '-tmplim',
    #                 tmp_name, '-outim', out_name, '-oci',
    #                 conv_name, '-tl', '-100', '-il', '-100',
    #                 '-nrx', '2', '-nry', '2',
    #                 '-nsx', '15', '-nsy', '15',
    #                 '-ng','4','7','0.70','6','1.50','4','3.00','3','6.0']

    #filenames, normalize to src and force convolution on tmp
    flags = ['-n', 'i', '-inim', src_name, '-tmplim', tmp_name, '-outim',
             out_name, '-oci', conv_name, '-c', 't', '-ko', '2']
    #artifact mask files
    if tmp_mask is not None:
        flags += ['-tmi', tmp_mask]
    if src_mask is not None:
        flags += ['-imi', src_mask]
    flags += ['-mins', '1']
    #negative and saturation limits
    flags += ['-tl', str(tmp_neg), '-il', str(src_neg)]
    flags += ['-tu', str(tmp_sat), '-iu', str(src_sat)]
    #seeing based convolution kernels
    if fwhm is not None:
        flags += ['-r', str(2.5*fwhm/2.0)]
        flags += ['-rss', str(3.0*fwhm)]
    #make sure each stamp is ~2.5', imsize given in arcsec
    flags += ['-nsx', str(np.floor(imsize/(2.5*60)))]
    flags += ['-nsy', str(np.floor(imsize/(2.5*60)))]
    #gaussians with which to compose kernel
    if fwhm is not None:
        flags += ['-ng','3','6',str(fwhm/2.0),'4',str(fwhm),'2',str(2*fwhm)]
    else:
        flags += ['-ng','3','6','3.00','4','6.00','2','12.0']
        
    #call hotpants
    subprocess.call(['hotpants'] + flags)

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
        remaped_ref = run_wcsremap(ref_name, src_name2, tmpdir)

        #hotpants arguments
        default_flags = ['hotpants', '-inim', src_name, '-tmplim',
                         ref_name, '-outim', out_name, '-oci',
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
def make_diff_image(src_name, ref_name, out_name, conv_name, fwhm=None, imsize=3692.8, tmp_sat=40000, src_sat=45000, tmp_neg=-100, src_neg=-100, tmp_mask=None, src_mask=None, tmpdir="DITemp", delete_temp=True):
    try:
        
        import os
        
        #make temporary directory
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)
            
        #fix header info
        src_name2 = remove_tan_from_header(src_name, tmpdir)
        
        #remap reference file to source file coordinates
        remaped_ref = run_wcsremap(ref_name, src_name2, tmpdir)

        #hotpants arguments
        args = [src_name, remaped_ref, out_name, conv_name]
        if fwhm is not None:
            args += [fwhm]
        args += [imsize, tmp_sat, src_sat, tmp_neg, src_neg]
        if tmp_mask is not None:
            args += [tmp_mask]
        if src_mask is not None:
            args += [src_mask]
            
        #subtract remapped reference file from source file
        run_hotpants(*args)

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
    make_diff_image(args.src_name, args.ref_name, args.out_name, args.conv_name,
                    delete_temp=delete_temp)
