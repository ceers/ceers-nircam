# Wrapper for TweakReg in jwst v1.8+
#

__author__ = "Anton M. Koekemoer, STScI"
__version__ = "1.0.0"
__license__ = "BSD3"

# Version history
# 1.0.0 -- M Bagley, UT Austin, edited be compatible with jwst v1.8+ TweakReg
# 0.5.0 -- M Bagley, UT Austin, moving source detection parameters to config 
#          file adding routine to save wcs models
# 0.4.0 -- Micaela Bagley, UT Austin, edited to use configuration files for
#          each filter, reducing unused imports
# 0.3.0 -- Switching to using user-supplied source catalogs for each input
#          image, to use Source Extractor's windowed funtions for centroids
# 0.2.0 -- Tweaking the source detection using Photutils



import os
import sys
from glob import glob
import argparse
from configparser import ConfigParser

import numpy as np
from astropy.io import fits
from astropy.table import Table

from math import *

from asdf import AsdfFile
from jwst.datamodels import ImageModel
# jwst TweakReg imports
from jwst.tweakreg import TweakRegStep
from jwst.associations import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base


#################################
### SET I/O AND OTHER SETUP HERE 
IODIR = 'calibrated'
# Directory containing absolute reference catalog
CATDIR = 'catalogs'
# Filename of Absolute reference catalog
ABS_REFCAT = 'CEERS_EGS_HST_v1.9_cat_radecmag.ecsv'
# Suffix of input files to be processed
FILE_SUFFIX = 'cal'
# Directory to save tweaked WCS models
TWEAKDIR = 'tweakreg_wcs'
#################################


def save_wcs(calimage):
    """save the tweakreged WCS for assignment to other reduction versions"""
    model = ImageModel(calimage)
    wcs = model.meta.wcs
    f = AsdfFile()
    f.tree['filename'] = model.meta.filename
    f.tree['wcsinfo'] = model.meta.wcsinfo.instance
    f.tree['wcs'] = wcs
    # check that output dir exists
    try:
        os.makedirs(TWEAKDIR)
    except FileExistsError:
        pass
    f.write_to(os.path.join(TWEAKDIR,
               os.path.basename(calimage.replace('fits','asdf'))))


def make_image_catalogs(imgfile_list, params, rerun_catalogs=False):
    """Run Source Extractor on individual images to get input source catalogs

    First splits each image into sci and rms extensions. 

    Runs SE with hard coded source detection parameters. These may need 
    to be tweaked for shallower or deeper observations, but they work well
    for CEERS Epoch 1 imaging - detecting plenty of real sources and avoiding
    too many spurious sources.

    Converts the catalog files into the ecsv format needed by TweakReg,
    using the XWIN_IMAGE,YWIN_IMAGE for the centroid positions.

    Args:
        imgfile_list (list, str): list of image file names
        params (dict): dictionary of parameters from config file, used for 
            SE source detection
        rerun_catalogs (Optional [bool]): If True, rerun SE to generate 
            source catalogs even if the catalog files already exist
    Outputs:
        *sci.fits: science extension for each input image
        *rms.fits: RMS (ERR) extension for each input image
        *cat.txt: output catalog file from SE for each input image
        *cat.idxy.ecsv: simplified output catalog for each input image, 
            including source ID, windowed centroid positions, and AUTO flux,
            and saved as an ecsv for TweakReg
    """

    sefiles = ['se.config', 'se.conv', 'se.nnw', 'se.outputs']

    for filename in sefiles:
        if (not os.path.exists(filename)):
            print('Missing %s, copy to working directory'%filename)
            exit()

    for imgfile in imgfile_list:
        imgfile_sci  = imgfile[:-5] + '_sci.fits'
        imgfile_rms  = imgfile[:-5] + '_rms.fits'
        imgfile_cat  = imgfile[:-5] + '_cat.txt'
        catfile_idxy = imgfile_cat = imgfile[:-5] + '_cat.idxy.ecsv'

        if (not os.path.exists(catfile_idxy)):
            data_sci = fits.getdata(imgfile, 'SCI')
            data_rms = fits.getdata(imgfile, 'ERR')
            # for now, best to just use a flat value for rms
            data_rms[:,:] = np.median(data_rms[np.where(data_rms != 0)])
   
            fits.writeto(imgfile_sci, data_sci, overwrite=True)
            fits.writeto(imgfile_rms, data_rms, overwrite=True)

        if (os.path.exists(catfile_idxy)) & (not rerun_catalogs):
            print('%s exists, skipping'%catfile_idxy)
            continue

        # SE command
        run_str = 'sex  '+imgfile_sci + \
                  '  -c  se.config -WEIGHT_IMAGE '+imgfile_rms + \
                  ' -WEIGHT_TYPE MAP_RMS' + \
                  ' -CATALOG_NAME '+imgfile_cat+ \
                  ' -PIXEL_SCALE 0.03  '+ \
                  ' -SEEING_FWHM 0.1 '  + \
                  ' -GAIN 1 ' + \
                  ' -MAG_ZEROPOINT    21. '  + \
                  ' -BACK_SIZE        128'   + \
                  ' -BACK_FILTERSIZE  4'     + \
                  ' -BACKPHOTO_TYPE   LOCAL' + \
                  ' -BACKPHOTO_THICK  64'    + \
                  ' -DETECT_THRESH '  +params['detect_thresh'] + \
                  ' -ANALYSIS_THRESH '+params['detect_thresh'] + \
                  ' -DETECT_MINAREA ' +params['detect_minarea'] + \
                  ' -DEBLEND_NTHRESH '+params['deblend']

        a = os.system('/bin/rm -f '+imgfile_cat)
        a = os.system(run_str)
      
        # Convert this to the format needed by tweakreg
        catalog = Table.read(imgfile_cat, format='ascii.sextractor')
        catxy = Table([catalog['NUMBER'], catalog['XWIN_IMAGE'], 
                       catalog['YWIN_IMAGE'], catalog['FLUX_AUTO']],  
                      names=['id', 'xcentroid', 'ycentroid', 'flux'])
        catxy['xcentroid']  -=  1.
        catxy['ycentroid']  -=  1.
        catxy.write(catfile_idxy, format='ascii.ecsv', overwrite=True)


def extract_results(tweakreg_logfile, imgfile_list, sca_n, results):
    """Extract the fit results from the output tweakreg logfile

    """
    f = open(tweakreg_logfile,'r')
    lines = f.readlines()
    f.close()
    nlines = len(lines)

    # First look for relative alignment results 
    # this is hardcoded to just look for x,y shifts, no rotation
    for imgfile in imgfile_list[1:]:
        #
        i = -1
        while (i < nlines-1):
            #
            i += 1
            line = lines[i]
            #
            if (line.find('stpipe.tweakreg - INFO - Aligning image catalog \'GROUP ID: '+os.path.basename(os.path.splitext(imgfile)[0])) >= 0):
                #
                result_found = False
                nmatch = -1
                #
                while ((not result_found) and (i < nlines-1)):
                    #
                    i += 1
                    line = lines[i]
                    line_split = str.split(line)
                    #
                    if (line.find('stpipe.tweakreg - INFO - XSH') >= 0):
                        xsh, ysh = float(line_split[8]), float(line_split[10])
                    #
                    if (line.find('stpipe.tweakreg - INFO - FIT RMSE') >= 0):
                        rmse, mae = float(line_split[9]), float(line_split[12])
                    #
                    if (line.find('stpipe.tweakreg - INFO - Final solution based on') >= 0):
                        nmatch = int(float(line_split[11]))
                    #
                    if (line.find('stpipe.tweakreg - WARNING - Not enough matches') >= 0):
                        xsh, ysh, rmse, mae, nmatch = 0, 0, 0, 0, 0
                    #
                    # Finally, at this point a "result" has been found, either way
                    #
                    if (nmatch > -1):
                        result_found = True
                        results.append('Relative alignment:   %46s  %6i  %12.8f  %12.8f  %12.8f  %12.8f\n'  % (imgfile, nmatch, rmse, mae,  xsh, ysh))


    # Now look for the absolute alignment results 
    # this is hardcoded to look for x,y shifts, rotation, and scale
    #
    # Note: unfortunately this logfile is cumulative, so each time this 
    # subroutine is called, it's necessary to skip the first several 
    # instances of "GROUP ID: 987654" until we reach the value of the 
    # counter "sca_n"
    #
    # It's rather kluddy but is just a consequence of the logfile being 
    # cumulative for all the SCA runs.
    #
    n = 0
    i = -1
    while (i < nlines-1):
        #
        i += 1
        line = lines[i]
        #
        if (line.find('stpipe.tweakreg - INFO - Aligning image catalog \'GROUP ID: 987654') >= 0):
            #
            n += 1
            if (n == sca_n):	# only do it for the current SCA, ie skip it for all the preceding instances in this logfile.
                #
                result_found = False
                nmatch = -1
                #
                while ((not result_found) and (i < nlines-1)):
                    #
                    i += 1
                    line = lines[i]
                    line_split = str.split(line)
                    #
                    if (line.find('stpipe.tweakreg - INFO - XSH') >= 0):
                        xsh, ysh, rot, scale = float(line_split[8]), float(line_split[10]), float(line_split[12]), float(line_split[14])
                    # 
                    if (line.find('stpipe.tweakreg - INFO - FIT RMSE') >= 0):
                        rmse, mae = float(line_split[9]), float(line_split[12])
                    #
                    if (line.find('stpipe.tweakreg - INFO - Final solution based on') >= 0):
                        nmatch = int(float(line_split[11]))
                    #
                    if (line.find('stpipe.tweakreg - WARNING - Not enough matches') >= 0):
                        xsh, ysh, rot, scale, rmse, mae, nmatch = 0, 0, 0, 0, 0, 0, 0
                    #
                    # Finally, at this point a "result" has been found, either way
                    #
                    if (nmatch > -1):
                        result_found = True
                        results.append('Absolute alignment:   %-46s  %6i  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f\n'  % (tweakreg_logfile[9:-4], nmatch, rmse, mae,  xsh, ysh, rot, scale))



def run_tweakregstep(visitstr, filt, save_results):
    """A wrapper to prepare input catalogs, run TweakReg, and extract results

    Passing each input catalog as part of the SCA asn file. This worked better
    than creating a file listing the user-provided catalogs and passing this
    with the `catfile` kwarg.

    Args:
        visitstr (str): String of program+observation+visit ID of images to run through TweakReg (for example: jw01345001001)
        filt (str): Filter to run through TweakReg. All calibrated images with visitstr and filt will be grouped together
        save_results (bool): If True, save results of fits by writing output *tweakreg.fits files. Set to False to save time while iterating 

    Outputs:
        [visitstr]_[filt]_nrc*.json - association files grouping input images
            by NIRCam detector
        tweakreg_[visitstr]_[filt]_nrc*.log - output log file from the relative
            and absolute astrometric fits for this group of images
        tweakreg_results_all.[viststr]_[filt].txt - overall output file 
            summarizing the relative and absolute astrometric fits for all 
            groups of images
    """
    # Set up the absolute astrometric reference catalog
    abs_refcat = os.path.join(CATDIR, ABS_REFCAT)

    # Set up other parameters for tweakreg, pull from config file
    cfg = ConfigParser()
    cfg.read('%s.cfg'%filt.lower())
    options = cfg.options(visitstr)
    params = {}
    for option in options:
        params[option] = cfg.get(visitstr, option)

    # start logfile fresh
    if (os.path.exists('tweakreg_out.log')):
        os.system('/bin/rm -f tweakreg_out.log')

    # list of images
    images = glob(os.path.join(IODIR,'%s*_%s.fits'%(visitstr,FILE_SUFFIX)))
    imgfile_list = sorted([x for x in images if fits.getheader(x)['FILTER'] == filt.upper()])

    # get list of detectors
    sca_list = []
    for imgfile in imgfile_list:
        sca = str.split(imgfile,'_')[3]
        if (not sca in sca_list): sca_list.append(sca)
    sca_list.sort()

    tweakreg_results = 'tweakreg_results_all.%s_%s.txt' % (visitstr, filt)

    results = []
    results.append('Absolute/Relative     Filename                                        NMATCH    RMSE(arcsec)  MAE(arcsec)   XSH(arcsec)   YSH(arcsec)   ROT(deg)      SCALE\n')
    results.append('------------------    ----------------------------------------------  ------    ----------    ----------    ----------    ----------    ----------    ----------\n')

    sca_n = 1
    for sca in sca_list:
        # making sca_n always 1 because the logfile is being reset each time
        # we call tweakreg, different from before, but this is easier than 
        # editing the absolute astrometry from the extract_results function
        #sca_n += 1
        images = glob(os.path.join(IODIR,'%s*_%s*_%s.fits'%(visitstr,sca,
                                                               FILE_SUFFIX)))
        imgfile_list = sorted([x for x in images if fits.getheader(x)['FILTER'] == filt.upper()])

        # Make SE catalogs for each image using windowed coords and 
        # save in format expected by TweakReg
        make_image_catalogs(imgfile_list, params)
       
        # Create asn file for visit
        asn_name = '_'.join([visitstr, filt, sca])
        asn_file = '.'.join([asn_name, 'json'])
        if (not os.path.exists(asn_file)):
            asn = asn_from_list.asn_from_list(imgfile_list, 
                                              rule=DMS_Level3_Base, 
                                              product_name=asn_name)
 
            # add in SE catalogs for tweakreg
            for m,member in enumerate(asn['products'][0]['members']):
                cat = asn['products'][0]['members'][m]['expname'].replace('.fits', '_cat.idxy.ecsv')
                asn['products'][0]['members'][m]['tweakreg_catalog'] = cat
       
            with open(asn_file, 'w') as outfile:
                name, serialized = asn.dump(format='json')
                outfile.write(serialized)

        tweakreg_logfile = 'tweakreg_%s.log'%asn_name

        # Calling tweakreg
        # harcoded parameters are ones that are unlikely to change
        tweakreg_output = TweakRegStep.call(asn_file,
                        ### Source finding parameters
                        # save_catalogs should be ignored for user-provided cats
                        save_catalogs       = False,
                        use_custom_catalogs = True,
                        catalog_format      = 'ecsv',
#                        catfile             = catfile,
                        kernel_fwhm         = float(params['kernel_fwhm']),
                        snr_threshold       = float(params['snr_threshold']),
                        # ignoring parameters related to interal 
                        # object detection
#                        sharplo
#                        sharphi
#                        roundlo
#                        roundhi
#                        brightest
#                        peakmax
#                        bkg_boxsize         = int(params['bkg_boxsize']),
                        ### Relative astrometry, optimize alignment order
                        enforce_user_order  = True,
                        ### Relative astrometry, reference catalog parameters
                        # This doesn't always work when set to True
                        expand_refcat       = False,
                        ### Relative astrometry, object matching parameters
                        minobj              = int(params['minobj']),
                        searchrad           = float(params['searchrad']),
                        # don't use 2dhst for relative alignment, should 
                        # be pretty good for back-to-back exposures
                        use2dhist           = False,	
                        separation          = float(params['separation']),
                        tolerance           = float(params['tolerance']),
#                        xoffset
#                        yoffset
                        ### Relative astrometry, catalog fitting parameters
                        fitgeometry         = 'shift',
                        nclip               = int(params['nclip']),
                        sigma               = float(params['sigma']),
                        ### Absolute astrometry fitting parameters
                        abs_refcat          = abs_refcat,
                        abs_minobj          = int(params['abs_minobj']),
                        abs_searchrad       = float(params['abs_searchrad']),
                        abs_use2dhist       = True,
                        # reducing abs_separation can bring in more
                        # more objects to match
                        abs_separation      = float(params['abs_separation']),
                        abs_tolerance       = float(params['abs_tolerance']),
                        abs_fitgeometry     = params['abs_fitgeom'],
                        abs_nclip           = int(params['abs_nclip']),
                        abs_sigma           = float(params['abs_sigma']),
                        save_abs_catalog    = False,
                        save_results        = save_results,
                        logcfg              = 'tweakreg_log.cfg')

        # Find log file to parse output
        # log filename is set in tweakreg-log.cfg
        if (os.path.exists('tweakreg_out.log')):
            os.system('/bin/cp -p tweakreg_out.log  '+tweakreg_logfile)
            extract_results(tweakreg_logfile, imgfile_list, sca_n, results)

        if save_results:
            # save the tweaked WCS for assignment to other reduction versions
            for image in imgfile_list:
                tweakimage = os.path.basename(image).replace(
                                          '%s.fits'%FILE_SUFFIX,'tweakreg.fits')
                save_wcs(tweakimage)

    outfile = open(tweakreg_results,'w')
    outfile.writelines(results)
    outfile.close()


def main():
    parser = argparse.ArgumentParser(description='A wrapper for TweakReg')
    parser.add_argument('visitstr', type=str,
                        help='String of program+observation+visit ID of images to run through TweakReg (for example: jw01345001001)')
    parser.add_argument('filt', type=str,
                        help='Filter to run through TweakReg. All calibrated images in IODIR with visitstr, FILE_SUFFIX and filt will be grouped together. Wrapper looks for config file called [filt].cfg.')
    parser.add_argument('--save_results', action='store_true',
                        help='Save results of fit. This takes extra time and is not recommended while iterating on parameters.')
    args = parser.parse_args()

    run_tweakregstep(args.visitstr, args.filt, args.save_results)


if __name__ == '__main__':
    main()


