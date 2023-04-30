__author__ = "Micaela B. Bagley, UT Austin"
__version__ = "0.2.0"
__license__ = "BSD3"

# Version history
# 0.2.0 -- Adding option to provide a manual scaling factor instead of 
#          fitting for scaling factor. This becomes useful to override the 
#          scaling for images where bright sources, artifacts, etc. cause 
#          fit_wisp_feature to over or under fit 
# 0.1.5 -- Smoothing wisp template before fitting, and using the affected 
#          subregion of the detector to fit the scaling for wisp template
# 0.1.1 -- Using median absolute deviation instead of sigma-clipped variance

import os
import shutil
import argparse
import logging
from datetime import datetime
import numpy as np
from astropy.io import fits
from astropy.stats import median_absolute_deviation
from scipy.ndimage import binary_dilation, gaussian_filter
import matplotlib.pyplot as plt
# photutils
from photutils.segmentation import (make_source_mask, detect_threshold, 
                                    detect_sources)
from photutils.background import Background2D
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma

# jwst-related imports
from jwst.datamodels import ImageModel, FlatModel
from jwst.datamodels import SegmentationMapModel
from jwst.flatfield.flat_field import do_correction
from stdatamodels import util
import crds
# Pipeline 
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

#################################
### SET I/O AND OTHER SETUP HERE 
INPUTDIR = 'calibrated'
OUTPUTDIR = 'calibrated'
# Directory containing wisp templates
WISPDIR = 'wisps/wisps_2022_08_26'
#################################


def create_seg_map(data, segname='segm.fits', saveseg=False):
    """Detects sources in the calibrated files to create a source mask.

    Source detection parameters are set for a relatively cold run. The 
    goal is to detect and mask the largest, brightest sources. If the 
    source detection is too aggressive, the wisp features themselves may
    be detected and masked.

    Args:
        data (jwst.datamodel): jwst image datamodel
        saveseg (Optional [bool]): Set to True to save the segmentation map
        segname (Optional [str]): Filename for output segmentation map

    Returns:
        segm.data (int array): segmentation map
    """
    ## source detection parameters
    # masking sources before background determination
    mask_nsigma = 5.
    mask_npixels = 5.
    mask_filter_fwhm = 1.
    mask_filter_size = 3
    # background 
    back_box_size = 256
    back_filter_size = (3,3)
    # detection parameters
    detect_thresh = 5.5 
    kernel_sigma = 1.5
    kernel_x_size = 3
    kernel_y_size = 3
    connectivity = 8
    detect_npixels = 15 
    # dilation kernel
    dilation_sigma = 7
    dilation_window = 15

    # prep image and err
    im = data.data
    err = data.err
    wcsinfo = data.meta.wcsinfo
    err[np.isnan(err)] = 1.e5
    im = np.asarray(im)

    # make a source mask
    log.info('creating source mask')
    mask = make_source_mask(im, mask_nsigma, mask_npixels,
                            filter_fwhm=mask_filter_fwhm,
                            filter_size=mask_filter_size)

    # use a coverage mask to identify missing data
    covmask = np.zeros(mask.shape, dtype=bool)
    covmask[err == 1.e5] = True

    # background
    log.info('background determination')
    bkg = Background2D(im, back_box_size, mask=mask, coverage_mask=covmask,
                       filter_size=back_filter_size)

    # detection threshold
    # either subtract off the background or include it in the threshold def
    threshold = bkg.background + (detect_thresh * err)

    # convolution kernel
    sigma = kernel_sigma * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=kernel_x_size, y_size=kernel_y_size)
    kernel.normalize()

    # detecting sources (convolution not used for cold run)
    log.info('detecting sources')
    segm = detect_sources(im, threshold, npixels=detect_npixels,
                          connectivity=connectivity)
#                          kernel=kernel, connectivity=connectivity)

    # grow segmentation map
    sources = np.logical_not(segm.data == 0)
    dilation_kernel = Gaussian2DKernel(dilation_sigma, x_size=dilation_window,
                                       y_size=dilation_window)

    source_wings = binary_dilation(sources, dilation_kernel)
    segm.data[source_wings] = 1

    if saveseg:
        segm_model = SegmentationMapModel(segm.data)
        segm_model.update(data, only='PRIMARY')
#        segm_model.meta.wcs = data.meta.wcs
        segm_model.meta.wcsinfo = wcsinfo
        segm_model.save(segname)

    return segm.data


def calc_variance(data, template, coeff):
    """Calculates the absolute median deviation of wisp subtracted image.

    Determines the variance of the function: image - coefficient * template.
    Using the median absolute deviation squarred. This is not scaled to 
    represent the standard deviation of normally distributed data, as would 
    be appropriate for an error estimator. However, fit_wisp_feature() will 
    find the coefficient that minimizes this variance, and so the relative 
    values are what matter. 

    Args:
        data (float): image array of masked data values
        template (float): image array of wisp template
        coeff (float): coefficient for scaling wisp template

    Returns:
        var_mad (float): median absolute deviation squarred for given coeff
    """
    func = data - coeff * template
    sigma_mad = median_absolute_deviation(func, ignore_nan=True)
    var_mad = sigma_mad**2
    return var_mad


def fit_wisp_feature(image, origfilename, apply_flat=True, fit_scaling=True, scale=1.0):
    """Scale wisp template and subtract from rate file.

    The wisp template is scaled to match the strength of the feature 
    present in the image. This is done by minimizing the variance of 
    (image - coefficient * template). 

    Original rate file is renamed to avoid being overwritten. 

    Args:
        image (str): Filename of input rate file
        origfilename (str): Filename to rename original rate file 
        apply_flat (Optional [bool]): set to True to apply flat to image first
        fit_scaling (Optional [bool]): set to True to determine a scaling for
            wisp template by minimizing the variance of image - coeff*template
        scale (Optional [float]): set a scaling value to use instead of 
            fitting. This is ignored if fit_scaling is True

    Outputs:
        Wisp-subtracted rate file is saved to OUTPUTDIR/image. Original
        rate file is renamed to origfilename. Additionally, a figure showing
        the variance determined for an array of cofficients is saved as 
        OUTPUTDIR/*_wisp.pdf.
    """
    # plotting parameters
    plt.rcParams['axes.linewidth'] = 2.
    plt.rcParams['axes.labelsize'] = 20
    plt.rcParams['axes.formatter.useoffset'] = False
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    try:
        crds_context = os.environ['CRDS_CONTEXT']
    except KeyError:
        crds_context = crds.get_default_context()

    # only need to do wisp subtraction for subset of filters and detectors
    detectors = ['NRCA3','NRCA4','NRCB3','NRCB4']

    model = ImageModel(os.path.join(INPUTDIR,image))
    filt = model.meta.instrument.filter
    detector = model.meta.instrument.detector
    if (filt not in ['F150W','F200W']) | (detector not in detectors):
        log.info('No correction for %s %s (%s)'%(filt,detector,image))
        return

    # check that image has not already been corrected
    for entry in model.history:
        if 'Wisps: subtracted' in entry['description']:
            log.info('%s already corrected for wisps, exiting'%image)
            return

    log.info('Fitting wisps feature')
    log.info('Working on %s'%image)

    # apply the flat to get a cleaner meausurement of the striping
    if apply_flat:
        log.info('Applying flat to match wisp templates')
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

    # read in template and mask nans
    wisptemplate = os.path.join(WISPDIR,'wisps_%s_%s.fits'%(detector.lower(),
                                                            filt.upper()))
    wisptemp = fits.getdata(wisptemplate)
    wisptemp[np.isnan(wisptemp)] = 0
    wisptemp[model.data == 0] = 0

    # source and bad pixel masking is only necessary if scaling template
    if fit_scaling:
        # construct mask for median calculation
        mask = np.zeros(model.data.shape, dtype=bool)
        mask[model.dq > 0] = True

        # source detection
        segfilename = os.path.join(OUTPUTDIR, image.replace('_rate.fits',
                                                            '_segm_cold.fits'))
        segm = create_seg_map(model, segfilename, saveseg=True)
        wobj = np.where(segm > 0)
        mask[wobj] = True

        maskedim = model.data.copy()
        maskedim[mask] = 0

        # smooth template before variance scaling, seems to help
        wisptemp = gaussian_filter(wisptemp, 2.0, mode='constant', truncate=3.5)

        # consider subsets of image focused around wisps for variance scaling
        if detector.lower() == 'nrca3':
            x1 = 300
            x2 = 1300
            y1 = 1100
            y2 = 2046
        elif detector.lower() == 'nrca4': 
            x1 = 450
            x2 = 1450
            y1 = 0
            y2 = 600
        elif detector.lower() == 'nrcb3': 
            x1 = 700
            x2 = 1450
            y1 = 0
            y2 = 800
        elif detector.lower() == 'nrcb4': 
            x1 = 900
            x2 = 1600
            y1 = 850
            y2 = 2046
   
        imseg = maskedim[y1:y2,x1:x2]
        wispseg = wisptemp[y1:y2,x1:x2]

        log.info('fitting coefficients')    
        coeffs = np.arange(0.01, 1.2, 0.02)
        variance_mad = np.zeros(coeffs.shape[0])
        for i,c in enumerate(coeffs):
            variance_mad[i] = calc_variance(imseg, wispseg, c)

        fig,(ax1,ax2) = plt.subplots(2, 1, figsize=(12,8), tight_layout=True)
    
        # fit with a curve to base scaling off of trend, rather than scatter
        fit_mad = np.polyfit(coeffs, variance_mad, deg=2)
        pfit_mad = np.poly1d(fit_mad)
        ax1.plot(coeffs, pfit_mad(coeffs)*1e4, 'C0', lw=1.5)
        ax1.plot(coeffs, variance_mad*1e4, 'C0o', lw=1.5)

        # show difference between curve and measured variances
        diff = variance_mad - pfit_mad(coeffs)
        ax2.plot(coeffs, diff*1e6, 'C0', lw=1)

        m = np.argmin(pfit_mad(coeffs))
        minval = coeffs[m]
        for ax in [ax1,ax2]:
            ax.axvline(minval, color='k', ls='--', lw=1.5)

        log.info('%s - fit coefficient for %s = %.2f'%(image, 
                    os.path.basename(wisptemplate), minval))

        for ax in [ax1,ax2]:
            ax.axvline(minval, color='C0', ls='--', lw=1.5)

        log.info('%s - using coefficient for %s = %.2f'%(image, 
                    os.path.basename(wisptemplate), minval))

        ax2.set_xlabel('coefficient')
        ax2.set_ylabel(r'residuals (10$^{-6}$)')
        ax1.set_ylabel(r'var (from MAD, 10$^{-4}$)', labelpad=10)

        outplot = os.path.join(OUTPUTDIR, 
                              image.replace('_rate.fits','_wisp.pdf'))
        fig.savefig(outplot)
#        plt.show()

    else:
        minval = scale
        log.info('%s - using coefficient for %s = %.2f'%(image, 
                    os.path.basename(wisptemplate), minval))

    # close model and open a clean version to clear anything we've done
    # to it (ie, flat fielding)
    model.close()
    del wisptemp

    # copy original
    log.info('Copying input to %s'%origfilename)
    shutil.copy2(os.path.join(INPUTDIR,image), 
                 os.path.join(OUTPUTDIR,origfilename))

    wisptemp = fits.getdata(wisptemplate)
    wisptemp[np.isnan(wisptemp)] = 0
    wisptemp[model.data == 0] = 0

    model = ImageModel(os.path.join(INPUTDIR,image))
    # subtract out wisp
    temp = minval * wisptemp
    corr = model.data - temp
    model.data = corr
   
    # add history entry
    time = datetime.now()
    stepdescription = 'Wisps: subtracted %s, scale=%.2f; wispsub.py %s'%(wisptemplate,minval,time.strftime('%Y-%m-%d %H:%M:%S'))
    substr = util.create_history_entry(stepdescription)
    model.history.append(substr)

    output = os.path.join(OUTPUTDIR, image)
    model.save(output)
    log.info('cleaned image saved to %s'%output)
    model.close()


def main():
    # affected detectors / filters
    # nrca3, nrca4, nrcb3, nrcb4
    # F150W, F200W

    parser = argparse.ArgumentParser(description='Scale wisp template and subtract from rate file')
    parser.add_argument('image', type=str, 
                        help='Filename of rate image for wisp subtraction')
    parser.add_argument('--fit_scaling', action='store_true',
                        help='Determine the optimal scaling for wisp template? Default is True')
    parser.add_argument('--scale', type=float, default=0.0,
                        help='Scaling factor to use instead of fitting (ignored if --fit_scaling=True). Default is 0.')
    args = parser.parse_args()

    origfilename = os.path.basename(args.image).replace('.fits','_prewisp.fits')

    fit_wisp_feature(args.image, origfilename, apply_flat=True, 
                     fit_scaling=args.fit_scaling, scale=args.scale)



if __name__ == '__main__':
    main()


