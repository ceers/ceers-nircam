__author__ = "Micaela B. Bagley, UT Austin"
__version__ = "0.1.0"
__license__ = "BSD3"

import os
import sys
from glob import glob
import argparse
import numpy as np

from jwst.datamodels import ImageModel
import asdf

#################################
### SET I/O AND OTHER SETUP HERE 
INPUTDIR = 'calibrated'
OUTPUTDIR = 'calibrated'
# Directory containing Tweakreg WCS solutions
TWEAKDIR = 'tweakreg_wcs'
# Suffix of input files to be processed
FILE_SUFFIX = 'cal'
# Suffix of output images with updated WCS
OUTPUT_SUFFIX = 'tweakreg'
#################################


def update_wcs(cal):
    """Update header with tweaked WCS saved from previous TweakReg run

    Args: 
        cal (str): filename of input cal file

    Outputs:
        - Output image with updated WCS model with suffix set by 
            OUTPUT_SUFFIX, ready for outlier detection
    """
    # open with jwst datamodels
    model = ImageModel(os.path.join(INPUTDIR,cal))
    
    # find tweakregged asdf file for this image
    base = cal.split('_%s.fits'%FILE_SUFFIX)[0]
    print('%s base: %s'%(cal, base))
    # directory of asdf wcs from tweakregged images
    tweakwcs = os.path.join(TWEAKDIR, '%s_tweakreg.asdf'%base)

    tweakreg = asdf.open(tweakwcs)

    wcs = tweakreg['wcs']
    wcsinfo = tweakreg['wcsinfo']

    print('%s updating wcs from %s'%(cal,tweakwcs))
    model.meta.wcs = wcs
    model.meta.wcsinfo = wcsinfo
    model.meta.cal_step.tweakreg = 'COMPLETE'

    # save output
    model.save(os.path.join(OUTPUTDIR, cal.replace('_%s.fits'%FILE_SUFFIX,
                                                   '_%s.fits'%OUTPUT_SUFFIX)))


def main():
    parser = argparse.ArgumentParser(description='Update WCS from a previous run of TweakReg')
    parser.add_argument('--image', type=str,
                        help='Filename of a single calibrated image')
    parser.add_argument('--all_images', action='store_true',
                        help='Run on all images in INPUTDIR with suffix FILE_SUFFIX')
    args = parser.parse_args()

    if args.all_images:
        images = glob(os.path.join(INPUTDIR, '*_%s.fits'%FILE_SUFFIX))

        for image in images:
            update_wcs(os.path.basename(image))

    else:
        update_wcs(args.image)


if __name__== '__main__':
    main()



