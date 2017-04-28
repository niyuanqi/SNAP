#################################################################
# Name:     LCgen.py                                            #
# Author:   Yuan Qi Ni                                          #
# Date:     April 28, 2017                                      #
# Function: Program applies Source Extractor to fits images.    #
#           Automates process of locating parameter files,      #
#           one line application.                               #
#################################################################

#function: apply SExtractor
def AutoSEx(src_name, out_name='none', config_name='none',
            check_out='NONE', check_name='none', xml_out='N',
            xml_name='none'):

    import subprocess

    #did user specify output name?
    if out_name == 'none':
        out_name = src_name[:-5]+".cat"
    #did user specify check image name?
    if args.check_name == 'none':
        check_name = src_name[:-5]+"."+check_out+".fits"
    #did user specify xml file name?
    if args.xml_name == 'none':
        xml_name = src_name[:-5]+".xml"
        
    #did user specify config file path?
    if args.config_name == 'none':
        #find path in which SExtractor is located
        sex_path = subprocess.check_output("which sex", shell=True)[:-5]
        sex_path = sex_path + "/../share/sextractor/"
        config_name = sex_path + "default.sex"
        params_name = sex_path + "default.param"
        filter_name = sex_path + "default.conv"
        subprocess.call(['sex', src_name, '-c', config_name,'-CATALOG_NAME', out_name, '-CHECKIMAGE_TYPE', check_out, '-CHECKIMAGE_NAME', check_name, '-WRITE_XML', xml_out, '-XML_NAME', xml_name, '-PARAMETERS_NAME', params_name, "-FILTER_NAME", filter_name])
    else:
        subprocess.call(['sex', src_name, '-c', config_name,'-CATALOG_NAME', out_name, '-CHECKIMAGE_TYPE', check_out, '-CHECKIMAGE_NAME', check_name, '-WRITE_XML', xml_out, '-XML_NAME', xml_name])

#main function
if __name__ == "__main__":
    
    import argparse
    import subprocess
    
    #receive arguments
    parser = argparse.ArgumentParser(description='Applies source extractor. Try sex -d to see what parameters can be additionally given to configure SExtractor.')
    parser.add_argument('src_name', type=str, help='input fits file name')
    parser.add_argument("--config_name", type=str, default='none', help="Configuration file name which specifies all other configuration related files. Routine by default takes config files in SExtractor path.")
    parser.add_argument("--out_name", type=str, default='none', help="Output file name, default is input_name.cat")
    parser.add_argument("--check_out", type=str, default='none', help="SExtractor outputs check image of given type")
    parser.add_argument("--check_name", type=str, default='none', help="Output check image name, default is input_name.<check_type>.fits")
    parser.add_argument("--xml_out", action='store_true', help="SExtractor outputs xml file.")
    parser.add_argument("--xml_name", type=str, default='none', help="Output xml file name, default is input_name.xml")
    args = parser.parse_args()

    #did user specify check image?
    if args.check_out == 'none':
        check_out = 'NONE'
    else:
        check_out = args.check_out
    #did user specify xml file?
    if args.xml_out:
        xml_out = 'Y'
    else:
        xml_out = 'N'
        
    #run SExtractor routine
    AutoSEx(args.src_name, args.out_name, args.config_name, check_out, args.check_name,
            xml_out, args.xml_name)
