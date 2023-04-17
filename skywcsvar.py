__author__ = "Micaela B. Bagley, UT Austin"
__version__ = "0.5.1"
__license__ = "BSD3"

# Version history
# 0.5.1 -- Using global variables for image suffix 
# 0.5.0 -- Incorporating 2D background model determination into process(). 
#          These background-subtracted images are used for the variance 
#          rescaling
# 0.4.0 -- Adding VAR_RNOISE rescaling to include sky rms
# 0.3.0 -- Adding step to fill holes in variance maps
# 0.2.0 -- Adding WCS update 

import os
import sys
from glob import glob
import pprint
import argparse
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip, biweight_location
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

from jwst.datamodels import ImageModel
import asdf
import background_subtraction
import compute_cal_sky_variance


#################################
### SET I/O AND OTHER SETUP HERE 
INPUTDIR = 'calibrated'
OUTPUTDIR = 'calibrated'
# Directory containing Tweakreg WCS solutions
TWEAKDIR = 'tweakreg_wcs'
# Suffix of input files to be processed
FILE_SUFFIX = 'crf'
# Suffix of 2D background subtracted image
BKG_SUFFIX = 'bkgsub1'
# Suffix of output images ready for mosaicking
OUTPUT_SUFFIX = 'match'
#################################


def gaussian(x, a, mu, sig):
    return a * np.exp(-(x-mu)**2/(2*sig**2))


def fit_sky(data, plot_sky=False, ax=None, color='C1', **kwargs):
    """Fit distribution of sky fluxes with a Gaussian"""
    bins = np.arange(-1.5, 2.0, 0.001)
    h,b = np.histogram(data, bins=bins)
    bc = 0.5 * (b[1:] + b[:-1])
    binsize = b[1] - b[0]

    p0 = [10, bc[np.argmax(h)], 0.01]
    popt,pcov = curve_fit(gaussian, bc, h, p0=p0)

    if plot_sky:
        ax.plot(bc, gaussian(bc, *popt), color, **kwargs)

    return popt[1]


def run_background_subtraction(fitsfile):
    """A wrapper to call background subtraction for an individual image.

    This will estimate a full background model for an image, more than 
    just a single pedestal. This background subtraction is used to calculate
    the sky variance for rescaling the variance maps. We do subtract this 
    background model from images before mosaicking because it can oversubtract
    flux in regions around extended/low surface brightness sources.

    Args:
        fitsfile (str): filename of image to be background-subtracted
    """
    # Insantiate the SubtractBackground object. Set the output suffix
    bs = background_subtraction.SubtractBackground()
    bs.suffix = BKG_SUFFIX
    bs.replace_sci = True

    # Print out the parameters being used
    pprint.pprint(bs.__dict__)
    print("")

    bs.do_background_subtraction(INPUTDIR, fitsfile)


def process(cal, plot_sky=False):
    """Performs 5 additional corrections to cal files before mosaicking

    process() performs 5 steps:
      1. Run a background subtraction routine that creates a 2D model of 
         the background in an image. This will be used to calculate the 
         sky variance. The routine saves a background-subtracted version of 
         the cal file with he suffix set by BKG_SUFFIX. This does not alter 
         the original cal file because we only want to subtract a pedestal 
         from the images before stacking. This background subtraction 
         routine can oversubtract in regions around extended/low surface 
         brightness sources

      2. Calculate and subtract a single pedestal value from the cal file

      3. Update the WCS if necessary. Copies the wcs model from the asdf
         files saved following the Tweakred step to the cal file. This 
         allows one use the updated/corrected WCS in many future reductions 
         without rerunning Tweakreg

      4. Rescale the variance maps. Determines a robust sky variance in 
         each of the *_[BKG_SUFFIX].fits image and scales the VAR_RNOISE
         array to reproduce this value. The VAR_RNOISE arrays are used for 
         inverse variance weighting during drizzling, so this step ensures 
         that the resulting error arrays will include the rms sky fluctuations

      5. Fix holes in the variance maps. In v1.7.2 and previous versions of 
         the jwst Pipeline, we found that known bad pixels had values of 
         exactly zero in the variance arrays. When the dithered images were 
         coadded, these areas in the output error array had relatively low 
         rms compared with the average. In these areas, the input error array 
         with the missing data or bad pixel did not contribute to the rms of 
         the affected pixel, creating 'holes' in the output error arrays. We 
         set these to infinity in the individual variance maps to correctly
         down-weight bad pixels during drizzling.

    Args:
        cal (str): filename of input cal file
        plot_sky (Optional [bool]): if True, plots binned pixel fluxes and
            the Gaussian fit when determining the sky pedestal

    Outputs:
        - 2D Background subtracted image with suffix set by BKG_SUFFIX
        - Output image will all corrections performed with suffix set by 
            OUTPUT_SUFFIX, ready to be combined in a mosaic
    """

    ### Run background subtraction, which will be used for sky variance calc
    run_background_subtraction(cal)


    ### Now calculate the sky pedestal value in the image for subtraction
    # open with jwst datamodels
    model = ImageModel(os.path.join(INPUTDIR,cal))
    dq = model.dq
    sci = model.data
    
    # read in segmentation map created during 1/f correction
    # image3_part1.asdf adds the association id a3001 to the filename
    segmap = os.path.join(INPUTDIR, 
                cal.replace('_a3001_%s.fits'%FILE_SUFFIX,'_rate_1fmask.fits'))
    seg = fits.getdata(segmap)
    w = np.where((dq == 0) & (seg == 0))
    data = sci[w]
    data = data.flatten()

    if plot_sky:
        bins = np.arange(-1.5, 2.0, 0.001)
        fig,ax = plt.subplots(1, 1, tight_layout=True, figsize=(15,8))
        ax.hist(data, bins=bins, color='k', alpha=0.3)
    else:
        ax = None

    try:
        sky = fit_sky(data, plot_sky=plot_sky, ax=ax, color='C1', 
                      alpha=0.5, lw=2)
    except RuntimeError as e:
        print('!!! Error %s !!!'%cal)
        err = open('errors.list', 'a')
        err.write('{}\t{}\n'.format(cal,e))
        err.close()
        sky = 0

    meddata = np.median(data)

    # iterate on sigma clipping
    clipped = sigma_clip(data, sigma=5, sigma_upper=0, sigma_lower=10, 
                         maxiters=5, masked=False)
    medclip = np.median(clipped)
    biweight = biweight_location(clipped)

    if plot_sky:
        ax.hist(clipped, bins=bins, color='C0', alpha=0.3)
        ax.axvline(meddata, color='C1', lw=1)
        ax.axvline(medclip, color='C2', lw=1)
        ax.axvline(biweight, color='C3', lw=1)

    try:
        sky = fit_sky(clipped, plot_sky=plot_sky, ax=ax, color='C2', 
                      alpha=0.6, lw=1)
    except RuntimeError as e:
        print('!!! Error %s !!!'%cal)
        err = open('errors.list', 'a')
        err.write('{}\t{}\n'.format(cal,e))
        err.close()
        sky = 0

    if plot_sky:
        plt.show()

    print('%s'%cal)
    print('  clipped median: %f'%medclip)
    print('  biweight background: %f'%biweight_location(clipped))
    print('  gaussian-fit background: %f'%sky)
    print('%s subtracting sky: %f'%(cal,sky))
    # subtract off sky    
    model.data -= sky

    model.meta.background.level = sky
    model.meta.background.subtracted = True
    model.meta.background.method = 'local'


    ### update WCS if necessary 
    # this should already have been done before outlier detection, but
    # do it here just in case 
    # find tweakregged asdf file for this image
    if model.meta.cal_step.tweakreg != 'COMPLETE':  
        # image3_part1.asdf adds the association id a3001 to the filename
        base = cal.split('_a3001_%s.fits'%FILE_SUFFIX)[0]
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


    ### rescale variance maps
    print('%s rescaling readnoise variance'%cal)
    # Insantiate the SubtractBackground object. Set the output suffix
    sv = compute_cal_sky_variance.ScaledVariance()

    # Print out the parameters being used
    pprint.pprint(sv.__dict__)

    # use the 2D background subtracted image
    fitsfile = cal.replace('%s.fits'%FILE_SUFFIX, '%s.fits'%BKG_SUFFIX)
    
    sv.read_file(INPUTDIR, fitsfile)
    # directly pull corrected readnoise, rather than writing to file
    sv.correct_the_variance()
    varcorr = sv.predicted_skyvar
    model.var_rnoise = varcorr


    ### fix holes in variance maps
    print('%s fixing variance map holes'%cal)
    rnoise = model.var_rnoise
    poisson = model.var_poisson
    flat = model.var_flat

    w = np.where(rnoise == 0)
    rnoise[w] = np.inf

    w = np.where(poisson == 0)
    poisson[w] = np.inf

    w = np.where(flat == 0)
    flat[w] = np.inf

    model.var_rnoise = rnoise
    model.var_poisson = poisson
    model.flat = flat
    print('success %s'%cal)


    # save output
    model.save(os.path.join(OUTPUTDIR, cal.replace('_%s.fits'%FILE_SUFFIX, 
                                                   '_%s.fits'%OUTPUT_SUFFIX)))
   
    print('finished: %s'%cal)


def main():
    parser = argparse.ArgumentParser(description='Subtrat sky pedestal from image, update WCS if necessary, rescale and fill holes in VAR_RNOISE maps')
    parser.add_argument('--image', type=str,
                        help='Filename of a single calibrated image')
    parser.add_argument('--all_images', action='store_true', 
                        help='Run on all images in INPUTDIR with suffix FILE_SUFFIX')
    args = parser.parse_args()

    if args.all_images:
        images = glob(os.path.join(INPUTDIR, '*_%s.fits'%FILE_SUFFIX))

        for image in images:
            measure_sky(image, plot_sky=False)

    else:
        measure_sky(args.image, plot_sky=False)


if __name__== '__main__':
    main()



