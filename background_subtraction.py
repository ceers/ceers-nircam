#  Iterative source removal and background subtraction
#

__author__ = "Henry C. Ferguson, STScI"
__version__ = "1.4.0"
__license__ = "BSD3"

# History
# 1.1.0 -- output the tier masks as a bitmask
#       -- Fixed bug where the number of masks was implicitly set to 4
# 1.1.1 -- Fixed bug where estimate of background rms in tier_mask was ignoring the mask
# 1.2.1 -- Appends the background-subtracted image as a separate extension rather than replacing the SCI
# 1.2.2 -- Makes replace_sci an option for the location of sky-subtracted image 
# 1.3.0 -- Optionally mask bits that are set in the DQ array, if it is there
# 1.3.1 -- Only pass tier_mask the mask for bad pixels and off-detector (not previous_mask)
# 1.4.0 -- Added clipped_ring_median to try to have the ring-median have less suppression in the outskirts of galaxies

import numpy as np
from astropy.io import fits
from astropy import stats as astrostats
from os import path
from dataclasses import dataclass

# Imports for background estimation
from photutils import (
    Background2D,  # For estimating the background
    MedianBackground, BiweightLocationBackground, SExtractorBackground, # Not all are used
    BkgIDWInterpolator, BkgZoomInterpolator)  # For interpolating background
from photutils.segmentation import detect_sources, make_source_mask
from photutils.utils import circular_footprint
from photutils.utils import ShepardIDWInterpolator as idw
from astropy.convolution import (
    convolve, convolve_fft, Box2DKernel, Tophat2DKernel,
    Ring2DKernel, Gaussian2DKernel)
from scipy.ndimage import median_filter
from astropy.wcs import WCS

#import dill # Just for debugging

# If datamodels are used to allow selecting DQ flags by name
from jwst.datamodels import dqflags 

@dataclass
class SubtractBackground:
    tier_nsigma: list = (3., 3., 3., 3.)
    tier_npixels = list = (15, 15, 5, 2)
    tier_kernel_size: list = (25, 15, 5, 2)
    tier_dilate_size: list = (0, 0, 0, 3) 
    faint_tiers_for_evaluation: list = (3,4) # numbering starts at 1
    ring_radius_in:  float = 40
    ring_width: float = 3
    ring_clip_max_sigma: float = 5.
    ring_clip_box_size: int = 100  
    ring_clip_filter_size: int = 3
    bg_box_size: int = 5
    bg_filter_size: int = 3
    bg_exclude_percentile: int = 90
    bg_sigma: float = 3
    plot_smooth: int = 0
    suffix: str = "bkgsub"
    interpolator: str = 'zoom'
    replace_sci: bool = False
    dq_flags_to_mask: list = ('SATURATED',)

    def open_file(self,directory,fitsfile):
        with fits.open(path.join(directory,fitsfile)) as hdu:
            sci = hdu[1].data
            err = hdu[2].data
            self.has_dq = False
            for h in hdu:
                if 'EXTNAME' in h.header:
                    if h.header['EXTNAME'] == 'DQ':
                        self.has_dq = True
                        self.dq = h.data
                        print(f"{fitsfile} has a DQ array")
        return sci,err
    
    # Convenience routine for inspecting the background and mask
    def plot_mask(self, scene, bkgd, mask, zmin, zmax, smooth=0, slices=None):
        '''Make a three-panel plot of:
             * the mask for the whole image,
             * the scene times the mask
             * a zoomed-in region, with the mask shown as contours
        '''
        if slices:
            rows = slices[0]
            cols = slices[1]
            _mask = mask[rows,cols]
            _scene = scene[rows,cols]
            _bkgd = bkgd[rows,cols]
        else:
            _mask = mask
            _scene = scene
            _bkgd = bkgd
        plt.figure(figsize=(20, 10))
        plt.subplot(131)
        plt.imshow(_mask, vmin=0, vmax=1, cmap=plt.cm.gray, origin='lower')
        plt.subplot(132)
        smooth = self.plot_smooth
        if smooth == 0:
            plt.imshow((_scene-_bkgd)*(1-_mask), vmin=zmin, vmax=zmax, origin='lower')
        else:
            smoothed = convolve((_scene-_bkgd)*(1-_mask), Gaussian2DKernel(smooth))
            plt.imshow(smoothed*(1-_mask), vmin=zmin/smooth, vmax=zmax/smooth,
                       origin='lower')
        plt.subplot(133)
        plt.imshow(_scene-_bkgd, vmin=zmin, vmax=zmax,origin='lower')
        plt.contour(_mask, colors='red', alpha=0.2)
    
    def replace_masked(self, sci, mask):
        sci_nan = np.choose(mask,(sci,np.nan))
        robust_mean_background = astrostats.biweight_location(sci_nan,c=6.,ignore_nan=True)
        sci_filled = np.choose(mask,(sci,robust_mean_background))
        return sci_filled
    
    def off_detector(self, sci, err):
        return np.isnan(err) # True if OFF detector, False if on detector

    def mask_by_dq(self):
        self.dqmask = np.zeros(len(self.dq),bool)
        for flag_name in self.dq_flags_to_mask:
            flagbit = dqflags.pixel[flag_name]
            self.dqmask = self.dqmask | (np.bitwise_and(self.dq,flagbit) != 0)
    
    def ring_median_filter(self, sci, mask):
        print(f"Ring median filtering with radius, width = ",end='')
        print(f"{self.ring_radius_in}, {self.ring_width}")
        sci_filled = self.replace_masked(sci,mask)
        ring = Ring2DKernel(self.ring_radius_in, self.ring_width)
        filtered = median_filter(sci, footprint=ring.array)
        return sci-filtered

    def clipped_ring_median_filter(self, sci, mask):
        # First make a smooth background (clip_box_size should be big)
        bkg = Background2D(sci,
              box_size = self.ring_clip_box_size,
              sigma_clip = astrostats.SigmaClip(sigma=self.bg_sigma),
              filter_size = self.ring_clip_filter_size,
              bkg_estimator = BiweightLocationBackground(),
              exclude_percentile = 90,
              mask = mask,
              interpolator = BkgZoomInterpolator())
        # Estimate the rms after subtracting this
        background_rms = astrostats.biweight_scale((sci-bkg.background)[~mask]) 
        # Apply a floating ceiling to the original image
        ceiling = self.ring_clip_max_sigma * background_rms + bkg.background
        # Pixels above the ceiling are masked before doing the ring-median filtering
        ceiling_mask = sci > ceiling
        print(f"Ring median filtering with radius, width = ",end='')
        print(f"{self.ring_radius_in}, {self.ring_width}")
        sci_filled = self.replace_masked(sci,mask | ceiling_mask)
        ring = Ring2DKernel(self.ring_radius_in, self.ring_width)
        filtered = median_filter(sci_filled, footprint=ring.array)
        return sci-filtered
    
    def tier_mask(self, img, mask, tiernum = 0):
        background_rms = astrostats.biweight_scale(img[~mask]) # Already has been ring-median subtracted
        # Replace the masked pixels by the robust background level so the convolution doesn't smear them
        background_level = astrostats.biweight_location(img[~mask]) # Already has been ring-median subtracted
        replaced_img = np.choose(mask,(img,background_level))
        convolved_difference = convolve_fft(replaced_img,Gaussian2DKernel(self.tier_kernel_size[tiernum]),allow_huge=True)
        # First detect the sources, then make masks from the SegmentationImage
        seg_detect = detect_sources(convolved_difference, 
                    threshold=self.tier_nsigma[tiernum] * background_rms,
                    npixels=self.tier_npixels[tiernum], 
                     mask=mask)
        if self.tier_dilate_size[tiernum] == 0:
            mask = seg_detect.make_source_mask()
        else:
            footprint = circular_footprint(radius=self.tier_dilate_size[tiernum])
            mask = seg_detect.make_source_mask(footprint)
        print(f"Tier #{tiernum}:")
        print(f"  kernel_size = {self.tier_kernel_size[tiernum]}")
        print(f"  tier_nsigma = {self.tier_nsigma[tiernum]}")
        print(f"  tier_npixels = {self.tier_npixels[tiernum]}")
        print(f"  tier_dilate_size = {self.tier_dilate_size[tiernum]}")
        print(f"  median of ring-median-filtered image = {np.median(img)}")
        print(f"  biweight rms of ring-median-filtered image  = {background_rms}")
        # For debugging #####################################################################
        # dill.dump(convolved_difference,open(f"convolved_difference{tiernum}.pkl","wb"))
        # dill.dump(img,open(f"bkgsub_img{tiernum}.pkl","wb"))
        # dill.dump(mask,open(f"bkgsub_masktier{tiernum}.pkl","wb"))
        #####################################################################################
        return mask

    def mask_sources(self, img, bitmask, starting_bit=1): 
        ''' Iteratively mask sources 
            Wtarting_bit lets you add bits for these masks to an existing bitmask
        '''
        print(f"ring-filtered background median: {np.median(img)}")
        first_mask = bitmask != 0
        for tiernum in range(len(self.tier_nsigma)):
            mask = self.tier_mask(img, first_mask, tiernum=tiernum)
            bitmask = np.bitwise_or(bitmask,np.left_shift(mask,tiernum+starting_bit))
        return bitmask
    
    def estimate_background(self, img, mask):
        bkg = Background2D(img, 
                    box_size = self.bg_box_size,
                    sigma_clip = astrostats.SigmaClip(sigma=self.bg_sigma),
                    filter_size = self.bg_filter_size,
                    bkg_estimator = BiweightLocationBackground(),
                    exclude_percentile = self.bg_exclude_percentile,
                    mask = mask,
                    interpolator = BkgZoomInterpolator())
        return bkg
    
    def estimate_background_IDW(self, img, mask):
        bkg = Background2D(img, 
                    box_size = self.bg_box_size,
                    sigma_clip = astrostats.SigmaClip(sigma=self.bg_sigma),
                    filter_size = self.bg_filter_size,
                    bkg_estimator = BiweightLocationBackground(),
                    exclude_percentile = self.bg_exclude_percentile,
                    mask = mask,
                    interpolator = BkgIDWInterpolator())
        return bkg

    def evaluate_bias(self, bkgd, err, mask):
        on_detector = np.logical_not(np.isnan(err)) # True if on detector, False if not
    
        mean_masked = bkgd[mask & on_detector].mean()
        std_masked = bkgd[mask & on_detector].std()
        stderr_masked = mean_masked/(np.sqrt(len(bkgd[mask]))*std_masked)
    
        mean_unmasked = bkgd[~mask & on_detector].mean()
        std_unmasked = bkgd[~mask & on_detector].std()
        stderr_unmasked = mean_unmasked/(np.sqrt(len(bkgd[~mask]))*std_unmasked)
        
        diff = mean_masked - mean_unmasked
        significance = diff/np.sqrt(stderr_masked**2 + stderr_unmasked**2)
        
        print(f"Mean under masked pixels   = {mean_masked:.4f} +- {stderr_masked:.4f}")
        print(f"Mean under unmasked pixels = "
              f"{mean_unmasked:.4f} +- {stderr_unmasked:.4f}")
        print(f"Difference = {diff:.4f} at {significance:.2f} sigma significance")
    
    # Customize the parameters for the different steps here
    def do_background_subtraction(self,datadir,fitsfile):
        # Background subtract all the bands
        print(fitsfile)
        sci, err = self.open_file(datadir,fitsfile)

        # Set up a bitmask
        bitmask = np.zeros(sci.shape,np.uint32) # Enough for 32 tiers

        # First level is for masking pixels off the detector
        off_detector_mask = self.off_detector(sci,err)
        #bitmask = np.bitwise_or(bitmask,np.left_shift(off_detector_mask,0))
        
        # Mask by DQ bits if desired and DQ file exists
        if self.has_dq:
            self.mask_by_dq()
            mask = off_detector_mask | self.dqmask
        else:
            mask = off_detector_mask 
        bitmask = np.bitwise_or(bitmask,np.left_shift(mask,0))

        # Ring-median filter 
        #filtered = self.ring_median_filter(sci, mask)
        filtered = self.clipped_ring_median_filter(sci, mask)
        
        # Mask sources iteratively in tiers
        bitmask = self.mask_sources(filtered, bitmask, starting_bit=1)
        mask = (bitmask != 0) 

        # Estimate the background using just unmasked regions
        if self.interpolator == 'IDW':
            bkg = self.estimate_background_IDW(sci, mask)
        else:
            bkg = self.estimate_background(sci, mask)
        bkgd = bkg.background

        # Subtract the background
        bkgd_subtracted = sci-bkgd

        # Evaluate
        print("Bias under bright sources:")
        self.evaluate_bias(bkgd,err,mask) # Under all the sources
        print("\nBias under fainter sources")
        faintmask = np.zeros(sci.shape,bool)
        for t in self.faint_tiers_for_evaluation:
            faintmask = faintmask | (np.bitwise_and(bitmask,2**t) != 0)
        self.evaluate_bias(bkgd,err,faintmask) # Just under the fainter sources
        
        # Write out the results
        prefix = fitsfile[:fitsfile.rfind('_')]
        outfile = f"{prefix}_{self.suffix}.fits"
        outpath = path.join(datadir,outfile)
        hdu = fits.open(path.join(datadir,fitsfile))
        wcs = WCS(hdu[1].header) # Attach WCS to it
        # Replace or append the background-subtracted image
        # Replace
        if self.replace_sci:
            hdu[1].data = bkgd_subtracted
        # Append 
        else:
            newhdu = fits.ImageHDU(bkgd_subtracted,header=wcs.to_header(),name='BKGSUB')
            hdu.append(newhdu)
        # Append an extension with the bitmask from the tiers of source rejection
        newhdu = fits.ImageHDU(bitmask,header=wcs.to_header(),name='TIERMASK')
        hdu.append(newhdu)
        # Write out the new FITS file
        hdu.writeto(outpath,overwrite=True)
        print(f"Writing out {outpath}")
        print("")
