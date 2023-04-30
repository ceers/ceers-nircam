# Measure and subtract the horizontal and vertical striping patterns 
# (1/f noise) from rate file

__author__ = "Micaela B. Bagley, UT Austin"
__version__ = "0.3.0"
__license__ = "BSD3"

# Version history
# 0.3.0 -- Adding option to provide a manual masking threshold instead of 
#          using default value
# 0.2.5 -- Reverting to full-row median if too many pixels in a given 
#          amp-row are masked
# 0.2.0 -- Fitting amplfier-by-amplifier instead of using the full row
# 0.1.6 -- Dilating the first source mask tier to get more of the wings of
#          bright sources
# 0.1.5 -- Using Harry's tiered source masking approach to more aggressively
#          mask source flux

import os
import shutil
import logging
from datetime import datetime
import argparse
import numpy as np
from astropy.io import fits
import astropy.stats as astrostats
from astropy.convolution import Ring2DKernel, Gaussian2DKernel, convolve_fft
from scipy.optimize import curve_fit
from scipy.ndimage import binary_dilation
from scipy.ndimage import median_filter
from photutils import  make_source_mask

# jwst-related imports
from jwst.datamodels import ImageModel, FlatModel, dqflags
from jwst.flatfield.flat_field import do_correction
from stdatamodels import util
import crds
# Pipeline 
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
# After first call to a jwst module, all logging will appear to come from 
# stpipe, and use the stpipe configuration/format.


#################################
### SET I/O AND OTHER SETUP HERE 
INPUTDIR = 'calibrated'
OUTPUTDIR = 'calibrated'
# Fraction of masked pixels in an amp-row that triggers switch to 
# full-row median
MASKTHRESH = 0.8
#################################


### from jwst/refpix/reference_pixels.py:
# NIR Reference section dictionaries are zero indexed and specify the values
# to be used in the following slice: (rowstart: rowstop, colstart:colstop)
# The 'stop' values are one more than the actual final row or column, in
# accordance with how Python slices work
NIR_reference_sections = {'A': {'top': (2044, 2048, 0, 512),
                                'bottom': (0, 4, 0, 512),
                                'side': (0, 2048, 0, 4),
                                'data': (0, 2048, 0, 512)},
                          'B': {'top': (2044, 2048, 512, 1024),
                                'bottom': (0, 4, 512, 1024),
                                'data': (0, 2048, 512, 1024)},
                          'C': {'top': (2044, 2048, 1024, 1536),
                                'bottom': (0, 4, 1024, 1536),
                                'data': (0, 2048, 1024, 1536)},
                          'D': {'top': (2044, 2048, 1536, 2048),
                                'bottom': (0, 4, 1536, 2048),
                                'side': (0, 2048, 2044, 2048),
                                'data': (0, 2048, 1536, 2048)}
                         }

### taking the reference rows/columns into account
NIR_amps = {'A': {'data': (4, 2044, 4, 512)},
            'B': {'data': (4, 2044, 512, 1024)},
            'C': {'data': (4, 2044, 1024, 1536)},
            'D': {'data': (4, 2044, 1536, 2044)}
            }


def gaussian(x, a, mu, sig):
    return a * np.exp(-(x-mu)**2/(2*sig**2))


def fit_sky(data):
    """Fit distribution of sky fluxes with a Gaussian"""
    bins = np.arange(-1, 1.5, 0.001)
    h,b = np.histogram(data, bins=bins)
    bc = 0.5 * (b[1:] + b[:-1])
    binsize = b[1] - b[0]

    p0 = [10, bc[np.argmax(h)], 0.01]
    popt,pcov = curve_fit(gaussian, bc, h, p0=p0)

    return popt[1]


def collapse_image(im, mask, dimension='y', sig=2.):
    """collapse an image along one dimension to check for striping.

    By default, collapse columns to show horizontal striping (collapsing
    along columns). Switch to vertical striping (collapsing along rows)
    with dimension='x' 

    Striping is measured as a sigma-clipped median of all unmasked pixels 
    in the row or column.

    Args:
        im (float array): image data array
        mask (bool array): image mask array, True where pixels should be 
            masked from the fit (where DQ>0, source flux has been masked, etc.)
        dimension (Optional [str]): specifies which dimension along which 
            to collapse the image. If 'y', collapses along columns to 
            measure horizontal striping. If 'x', collapses along rows to 
            measure vertical striping. Default is 'y'
        sig (Optional [float]): sigma to use in sigma clipping
    """
    # axis=1 results in array along y
    # axis=0 results in array along x
    if dimension == 'y':
#        collapsed = np.median(im, axis=1)
        res = astrostats.sigma_clipped_stats(im, mask=mask, sigma=sig, 
                                             cenfunc='median',
                                             stdfunc='std', axis=1)
    elif dimension == 'x':
#        collapsed = np.median(im, axis=0)
        res = astrostats.sigma_clipped_stats(im, mask=mask, sigma=sig, 
                                             cenfunc='median',
                                             stdfunc='std', axis=0)

    return res[1]
    

def masksources(image):
    """Detect sources in an image using a tiered approach for different 
       source sizes 

    image (str): filename of image

    """
    model = ImageModel(image)
    sci = model.data
    err = model.err
    wht = model.wht
    dq = model.dq

    # bad pixel mask for make_source_mask
    bpflag = dqflags.pixel['DO_NOT_USE']
    bp = np.bitwise_and(dq, bpflag)
    bpmask = np.logical_not(bp == 0)

    log.info('masking, estimating background')
    # make a robust estimate of the mean background and replace blank areas
    sci_nan = np.choose(np.isnan(err),(sci,err))
    # Use the biweight estimator as a robust estimate of the mean background
    robust_mean_background = astrostats.biweight_location(sci_nan, c=6.,
                                                          ignore_nan=True)
    sci_filled = np.choose(np.isnan(err),(sci,robust_mean_background))
    
    log.info('masking, initial source mask')
    # make an initial source mask
    ring = Ring2DKernel(40, 3)
    filtered = median_filter(sci, footprint=ring.array)

    log.info('masking, mask tier 1')
    # mask out sources iteratively
    # Try a reasonably big filter for masking the bright stuff
    convolved_difference = convolve_fft(sci-filtered,Gaussian2DKernel(25))
    mask1 = make_source_mask(convolved_difference, nsigma=3., npixels=15,
                             mask=np.isnan(err))
    # grow the largest mask 
    temp = np.zeros(sci.shape)
    temp[mask1] = 1
    sources = np.logical_not(temp == 0)
    dilation_sigma = 10
    dilation_window = 11
    dilation_kernel = Gaussian2DKernel(dilation_sigma, x_size=dilation_window,
                                       y_size=dilation_window)

    source_wings = binary_dilation(sources, dilation_kernel)
    temp[source_wings] = 1
    mask1 = np.logical_not(temp == 0)

    log.info('masksources: mask tier 2')
    # A smaller smoothing for the next tier
    convolved_difference = convolve_fft(sci-filtered,Gaussian2DKernel(15))
    mask2 = make_source_mask(convolved_difference, nsigma=3., npixels=15,
                             mask=mask1) | mask1

    log.info('masksources: mask tier 3')
    # Still smaller 
    convolved_difference = convolve_fft(sci-filtered,Gaussian2DKernel(5))
    mask3 = make_source_mask(convolved_difference, nsigma=3., npixels=5,
                             mask=mask2) | mask2

    log.info('masksources: mask tier 4')
    # Smallest 
    convolved_difference = convolve_fft(sci-filtered,Gaussian2DKernel(2))
    mask4 = make_source_mask(convolved_difference, nsigma=3., npixels=3,
                             mask=mask3,dilate_size=3) | mask3

    finalmask = mask4

    # save output mask
    outputbase = os.path.join(OUTPUTDIR, os.path.basename(image))
    maskname = outputbase.replace('.fits', '_1fmask.fits')
    log.info('masksources: saving mask to %s'%maskname)
    outmask = np.zeros(finalmask.shape, dtype=int)
    outmask[finalmask] = 1
    fits.writeto(maskname, outmask, overwrite=True)
    return outmask


def measure_fullimage_striping(fitdata, mask):
    """Measures striping in countrate images using the full rows.

    Measures the horizontal & vertical striping present across the 
    full image. The full image median will be used for amp-rows that
    are entirely or mostly masked out.

    Args:
        fitdata (float array): image data array for fitting
        mask (bool array): image mask array, True where pixels should be 
            masked from the fit (where DQ>0, source flux has been masked, etc.)

    Returns:
        (horizontal_striping, vertical_striping): 
    """

    # fit horizontal striping, collapsing along columns
    horizontal_striping = collapse_image(fitdata, mask, dimension='y')
    # remove horizontal striping, requires taking transpose of image
    temp_image = fitdata.T - horizontal_striping
    # transpose back
    temp_image2 = temp_image.T

    # fit vertical striping, collapsing along rows
    vertical_striping = collapse_image(temp_image2, mask, dimension='x')

    return horizontal_striping, vertical_striping


def measure_striping(image, origfilename, thresh=None, apply_flat=True, mask_sources=True, save_patterns=False):
    """Removes striping in rate.fits files before flat fielding.

    Measures and subtracts the horizontal & vertical striping present in 
    countrate images. The striping is most likely due to 1/f noise, and 
    the RefPixStep with odd_even_columns=True and use_side_ref_pixels=True
    does not fully remove the pattern, no matter what value is chosen for 
    side_smoothing_length. There is also residual vertical striping in NIRCam 
    images simulated with Mirage.

    The measurement/subtraction is done along one axis at a time, since 
    the measurement along x will depend on what has been subtracted from y.

    Note: 
        The original rate image file is copied to *_rate_pre1f.fits, and 
        the rate image with the striping patterns removed is saved to 
        *_rate.fits, overwriting the input filename

    Args:
        image (str): rate image filename, including full relative path
        origfilename (str): filename to rename original rate file
        thresh (Optional [float]): fraction of masked amp-row pixels above 
            which full row fit is used
        apply_flat (Optional [bool]): if True, identifies and applies the 
            corresponding flat field before measuring striping pattern. 
            Applying the flat first allows for a cleaner measure of the 
            striping, especially for the long wavelength detectors. 
            Default is True.
        mask_sources (Optional [bool]): If True, masks out sources in image
            before measuring the striping pattern so that source flux is 
            not included in the calculation of the sigma-clipped median.
            Sources are identified using the Mirage seed images.
            Default is True.
        save_patterns (Optional [bool]): if True, saves the horizontal and
            vertical striping patterns to files called *horiz.fits and 
            *vert.fits, respectively
    """
    try:
        crds_context = os.environ['CRDS_CONTEXT']
    except KeyError:
        crds_context = crds.get_default_context()
    
    # if thresh is not defined by user, use global default
    if thresh is None:
        thresh = MASKTHRESH

    # set up output filename, this will also be used for saving 
    # other outputs like the source mask and striping patterns
    outputbase = os.path.join(OUTPUTDIR, os.path.basename(image))

    model = ImageModel(image)
    log.info('Measuring image striping')
    log.info('Working on %s'%os.path.basename(image))

    # apply the flat to get a cleaner meausurement of the striping
    if apply_flat:
        log.info('Applying flat for cleaner measurement of striping patterns')
        # pull flat from CRDS using the current context
        crds_dict = {'INSTRUME':'NIRCAM', 
                     'DETECTOR':model.meta.instrument.detector, 
                     'FILTER':model.meta.instrument.filter, 
                     'PUPIL':model.meta.instrument.pupil, 
                     'DATE-OBS':model.meta.observation.date,
                     'TIME-OBS':model.meta.observation.time}
        flats = crds.getreferences(crds_dict, reftypes=['flat'], 
                                   context=crds_context)
        # if the CRDS loopup fails, should return a CrdsLookupError, but 
        # just in case:
        try:
            flatfile = flats['flat']
        except KeyError:
            log.error('Flat was not found in CRDS with the parameters: {}'.format(crds_dict))
            exit()

        log.info('Using flat: %s'%(os.path.basename(flatfile)))
        with FlatModel(flatfile) as flat:
            # use the JWST Calibration Pipeline flat fielding Step 
            model,applied_flat = do_correction(model, flat)
            
    # construct mask for median calculation
    mask = np.zeros(model.data.shape, dtype=bool)
    mask[model.dq > 0] = True
    
    # mask out sources
    if mask_sources:
        # first look for a source mask in OUTPUTDIR that already exists
        srcmask = outputbase.replace('.fits', '_1fmask.fits')
        if os.path.exists(srcmask):
            log.info('Using existing source mask %s'%srcmask)
            seg = fits.getdata(srcmask)
        else:
            log.info('Detecting sources to mask out source flux')
            seg = masksources(image)
        
        wobj = np.where(seg > 0)
        mask[wobj] = True

    # measure the pedestal in the unmasked parts of the image
    log.info('Measuring the pedestal in the image')
    pedestal_data = model.data[~mask]
    pedestal_data = pedestal_data.flatten()
    median_image = np.median(pedestal_data)
    log.info('Image median (unmasked and DQ==0): %f'%(median_image))
    try:
        pedestal = fit_sky(pedestal_data)
    except RuntimeError as e:
        log.error("Can't fit sky, using median value instead")
        pedestal = median_image
    else:
        log.info('Fit pedestal: %f'%pedestal)

    # subtract off pedestal so it's not included in fit  
    model.data -= pedestal

    # measure full pattern across image
    full_horizontal, vertical_striping = measure_fullimage_striping(model.data, 
                                                                    mask)

    horizontal_striping = np.zeros(model.data.shape)
    vertical_striping = np.zeros(model.data.shape)

    # keep track of number of number of times the number of masked pixels 
    # in an amp-row exceeds thersh and a full-row median is used instead
    ampcounts = []
    for amp in ['A','B','C','D']:
        ampcount = 0
        rowstart, rowstop, colstart, colstop = NIR_amps[amp]['data']
        ampdata = model.data[:, colstart:colstop]
        ampmask = mask[:, colstart:colstop]
        # fit horizontal striping in amp, collapsing along columns
        hstriping_amp = collapse_image(ampdata, ampmask, dimension='y')
        # check that at least 1/4 of pixels in each row are unmasked
        nmask = np.sum(ampmask, axis=1)
        for i,row in enumerate(ampmask):
            if nmask[i] > (ampmask.shape[1]*thresh):
                # use median from full row
                horizontal_striping[i,colstart:colstop] = full_horizontal[i]
                ampcount += 1
            else:
                # use the amp fit 
                horizontal_striping[i,colstart:colstop] = hstriping_amp[i]
        ampcounts.append('%s-%i'%(amp,ampcount))   

    ampinfo = ', '.join(ampcounts)
    log.info('%s, full row medians used: %s /%i'%(os.path.basename(image), 
                                                  ampinfo, rowstop-rowstart))

    # remove horizontal striping    
    temp_sub = model.data - horizontal_striping

    # fit vertical striping, collapsing along rows
    vstriping = collapse_image(temp_sub, mask, dimension='x')
    vertical_striping[:,:] = vstriping

    # save horizontal and vertical patterns 
    if save_patterns:
        fits.writeto(outputbase.replace('.fits', '_horiz.fits'), 
                     horizontal_striping, overwrite=True)
        fits.writeto(outputbase.replace('.fits', '_vert.fits'), 
                     vertical_striping, overwrite=True)

    model.close()
    
    # copy image 
    log.info('Copying input to %s'%origfilename)
    shutil.copy2(image, origfilename)

    # remove striping from science image
    with ImageModel(image) as immodel:
        sci = immodel.data
        # to replace zeros
        wzero = np.where(sci == 0)
        temp_sci = sci - horizontal_striping
        # transpose back
        outsci = temp_sci - vertical_striping
        outsci[wzero] = 0
        # replace NaNs with zeros and update DQ array
        # the image has NaNs where an entire row/column has been masked out
        # so no median could be calculated.
        # All of the NaNs on LW detectors and most of them on SW detectors
        # are the reference pixels around the image edges. But there is one
        # additional row on some SW detectors 
#        refpixflag = dqflags.pixel['REFERENCE_PIXEL']
#        wref = np.bitwise_and(immodel.dq, refpixflag)
#        outsci[np.where(wref)] = 0
        wnan = np.isnan(outsci)
        bpflag = dqflags.pixel['DO_NOT_USE']
        outsci[wnan] = 0
        immodel.dq[wnan] = np.bitwise_or(immodel.dq[wnan], bpflag)

        # write output
        immodel.data = outsci
        # add history entry
        time = datetime.now()
        stepdescription = 'Removed horizontal,vertical striping; remstriping.py %s'%time.strftime('%Y-%m-%d %H:%M:%S')
        # writing to file doesn't save the time stamp or software dictionary 
        # with the History object, but left here for completeness
        software_dict = {'name':'remstriping.py',
                         'author':'Micaela Bagley',
                         'version':'1.0',
                         'homepage':'ceers.github.io'}
        substr = util.create_history_entry(stepdescription, 
                                              software=software_dict)
        immodel.history.append(substr)
        log.info('Saving cleaned image to %s'%outputbase)
        immodel.save(outputbase)


def main():

    parser = argparse.ArgumentParser(description='Measure and remove horizontal and vertical striping pattern (1/f noise) from rate file')
    parser.add_argument('image', type=str,
                        help='Filename of rate image for pattern subtraction')

    parser.add_argument('--thresh', type=float, 
                        help='The threshold (fraction of masked pixels in an amp-row) above which to switch to a full-row median')
    parser.add_argument('--save_patterns', action='store_true',
                        help='Save the horizontal and vertical striping patterns as FITS files')
    args = parser.parse_args()
    image = os.path.join(INPUTDIR, args.image)

    # Original rate will be copied to INPUTDIR with suffix pre1f
    pre1f = image.replace('rate.fits', 'rate_pre1f.fits')

    measure_striping(image, pre1f, thresh=args.thresh, apply_flat=True, 
                     mask_sources=True, save_patterns=args.save_patterns)


if __name__ == '__main__':
    main()

    
