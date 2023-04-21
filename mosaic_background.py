__author__ = "Micaela B. Bagley, UT Austin"
__version__ = "0.2.0"
__license__ = "BSD3"

# History
# 0.2.0 -- use config files when setting different parameters for each filter

import os
from glob import glob
import argparse
from configparser import ConfigParser
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import SigmaClip
from photutils import Background2D, BkgZoomInterpolator
from photutils import BiweightLocationBackground

from jwst.datamodels import ImageModel
import stdatamodels

import background_subtraction

#################################
### SET I/O AND OTHER SETUP HERE 
IODIR = 'calibrated'
# Filename of merged source mask
OUTMASK = 'merged_mask.fits'
# Suffix of 2D background subtracted images using individual masks
BKG_SUFFIX = 'bkgsub1' 
# Suffix of 2D background subtracted images using merged masks
MERGED_BKG_SUFFIX = 'mbkgsub1'
#################################


def run_background_and_tiermask(image):
    """Create source mask and do initial round of background subtraction

    """
    # get background parameters for image filter
    filt = fits.getheader(image, 0)['FILTER']

    if filt in ['F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F410M', 'F444W']:
        section = 'nircam'
    else:
        section = filt.lower()

    cfg = ConfigParser()
    cfg.read('mosaic_background.cfg')
    options = cfg.options(section)
    params = {}
    for option in options:
        params[option] = cfg.get(section, option)

    # Insantiate the SubtractBackground object. Set the output suffix
    bs = background_subtraction.SubtractBackground()
    bs.suffix = BKG_SUFFIX
    bs.replace_sci = True
    bs.bg_box_size = int(params['bg_box_size'])
    bs.bg_filter_size = int(params['bg_filter_size'])
    bs.ring_radius_in = int(params['ring_radius_in'])
    bs.ring_width = int(params['ring_width'])

    tier_dilate_size = [int(x) for x in params['tier_dilate_size'].split(',')]
    bs.tier_dilate_size = tier_dilate_size

    tier_nsigma = [float(x) for x in params['tier_nsigma'].split(',')]
    bs.tier_nsigma = tier_nsigma
 
    tier_npixels = [int(x) for x in params['tier_npixels'].split(',')]
    bs.tier_npixels = tier_npixels

    # Print out the parameters being used
    pprint.pprint(bs.__dict__)
    print("")

    bs.do_background_subtraction(IODIR, image)


def merge_masks(bkgfiles):
    """ """
    mask = None

    for bkgfile in bkgfiles:
        with fits.open(bkgfile) as hdu:

            # HST and NIRCam file extensions are different
            if hdu[0].header['TELESCOP'] == 'HST':
                sci_extension = 0
            else:
                sci_extension = 'SCI'

            if hdu[sci_extension].header['FILTER'] == 'F277W':
                wcs = WCS(hdu[sci_extension].header)
            input_tiermask = hdu['TIERMASK'].data
            # Bordermask is bit 1...clear that
            this_source_mask = np.left_shift(np.right_shift(input_tiermask,1),1)
            if mask is None:
                mask = this_source_mask
            else:
               mask = mask | this_source_mask 

    merged_mask = mask
    mask = mask.astype(np.int32)

    hduout = fits.PrimaryHDU(mask, header=wcs.to_header())
    hduout.writeto(os.path.join(IODIR, OUTMASK), overwrite=True)


def run_final_background_subtraction(image):
    """ """
    # get tiermask from bgk-subtracted image to get bordermask specific 
    # to this filter image
    bkgimage = image.replace('i2d',BKG_SUFFIX)
    with fits.open(bkgimage) as hdumask:
        bordermask = hdumask['TIERMASK'].data == 1 

    merged_mask = fits.getdata(os.path.join(IODIR,OUTMASK))
    sourcemask = merged_mask | bordermask
    mask = sourcemask != 0

    hdu = fits.open(image)

    # HST and NIRCam file extensions are different
    if hdu[0].header['TELESCOP'] == 'HST':
        sci_extension = 0
    else:
        sci_extension = 'SCI'

    sci = hdu[sci_extension].data

    bkg = Background2D(sci,
                       box_size = 10,
                       sigma_clip = SigmaClip(sigma=3),
                       filter_size = 5,
                       bkg_estimator = BiweightLocationBackground(),
                       exclude_percentile = 90,
                       mask = mask,
                       interpolator = BkgZoomInterpolator())

    bkgsub = sci - bkg.background
    bkgsub = np.choose(bordermask, (bkgsub,0.))
    # Replace sci extension with background-subtracted version
    hdu[sci_extension].data = bkgsub

    # Append background and mask as extensions
    bkgd = fits.ImageHDU(bkg.background, header=wcs.to_header(), name='BKGD')
    hdu.append(bkgd)
    bkgmask = fits.ImageHDU(sourcemask.astype(np.int8), header=wcs.to_header(),
                            name='BKGMASK')
    hdu.append(bkgmask)

    # Write out the file
    outfile = image.replace('i2d', MERGED_BKG_SUFFIX)
    hdu.writeto(os.path.join(IODIR, outfile), overwrite=True)
    hdu.close()


def main():
    parser = argparse.ArgumentParser(description='Mask sources in mosaics across all filters and subtract a 2D background model from mosaics')
    parser.add_argument('pointing', type=str,
                        help='CEERS pointing ID, background subtract all mosaics for this pointing')
    parser.add_argument('--add_hst', action='store_true',
                        help='Include HST images in background subtraction and source masking')
    args = parser.parse_args()

    # get list of mosaics for background subtraction
    images = glob(os.path.join(IODIR, 'ceers_%s_*i2d.fits'%args.pointing))
    images.sort()
    # Include HST images
    if args.add_hst:
        hstimages = glob(os.path.join(HSTDIR, 
                         'egs_all_*_030mas_v1.9_%s_*;'*args.pointing))
        hstimages.sort()
        images = images + hstimages

    # first run background subtraction on each mosaic individually to 
    # mask sources
    for image in images:
        run_background_and_tiermask(os.path.join(IODIR, image))
            
    # now merge source masks from all available filters, using the BKG_SUFFIX
    # images that include the individual source masks as an extension
    bkgimages = glob(os.path.join(IODIR, '*_%s.fits'%BKG_SUFFIX))
    bkgimages.sort()
    merge_masks(bkgimages)

    # run final background subtraction on each image using merged mask
    for image in images:
        run_final_background_subtraction(image)



if __name__ == '__main__':
    main()


