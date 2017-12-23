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
#python make_image_diff_n.py N300-1.Q0.B.161030_0504.C.034140.005604N3646.0060.nh.fits N300-1.Q0.B.150626_1842-151019_1231.XCSA.005605N3646.00081.00081.FM30.BS0512.ALL.coadd.NEW.REF.fits diff.fits

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
def run_hotpants(src_name, template_name, out_name, conv_name):

    import subprocess

    #call hotpants
    subprocess.call(['hotpants', '-inim', src_name, '-tmplim',
                     template_name, '-outim', out_name, '-oci',
                     conv_name, '-tl', '-100', '-il', '-100',
                     '-nrx', '2', '-nry', '2',
                     '-nsx', '15', '-nsy', '15',
                     '-ng','4','7','0.70','6','1.50','4','3.00','3','6.0'])

#make difference image
def make_diff_image(src_name, ref_name, out_name, conv_name, delete_temp=True):
    try:
        
        import os
        
        #make temporary directory
        tmpdir = 'DITemp'
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)
            
        #fix header info
        src_name2 = remove_tan_from_header(src_name, tmpdir)
        
        #remap reference file to source file coordinates
        remaped_ref = run_wcsremap(ref_name, src_name2, tmpdir)
        
        #subtract remapped reference file from source file
        run_hotpants(src_name, remaped_ref, out_name, conv_name)

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
