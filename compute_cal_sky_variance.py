#!/usr/bin/env python
# Compute a rough inverse-variance map for drizzling by rescaling the VAR_RDNOISE
# by the measured sky variance in source-free 7x7 pixel regions
# Currently masks anything that is non-zero in the background mask. This could be tweaked...

__author__ = "Henry C. Ferguson, STScI"
__version__ = "0.1.0"
__license__ = "BSD3"

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits
from scipy import ndimage
import sys
from astropy.nddata import block_reduce
from astropy.stats import biweight_location, biweight_midvariance
from dataclasses import dataclass
from os import path

@dataclass
class ScaledVariance:
    block_size: int = 7
    mask_extension: int = 9

    def read_file(self,datadir,fitsfile):
        self.fitsfile = fitsfile
        self.datadir = datadir
        self.hdu = fits.open(path.join(datadir,fitsfile))
        self.sci = self.hdu[1].data
        self.var_rdnoise = self.hdu[6].data
        self.mask = self.hdu[self.mask_extension].data

    def compute_variance(self,img):
        blk_img = block_reduce(img,self.block_size)
        # Mask bins with any originally-masked pixels
        blk_mask = block_reduce(self.mask,self.block_size) != 0 
        unmasked_bins = blk_img[blk_mask == 0]
        variance = biweight_midvariance(unmasked_bins)
        return variance / self.block_size**2 # Return equivalent 1-pixel variance

    def masked_mean(self,img):
        blk_img = block_reduce(img,self.block_size)
        # Mask bins with any originally-masked pixels
        blk_mask = block_reduce(self.mask,self.block_size) != 0 
        unmasked_img = blk_img[blk_mask == 0]
        mean = biweight_location(unmasked_img)
        return mean / self.block_size**2 # because block_reduce sums by default

    def what_fraction_unmasked(self):
        # Mask bins with any originally-masked pixels
        blk_mask = block_reduce(self.mask,self.block_size) 
        on_detector = (~(blk_mask == 1)).sum()
        unmasked_bins = (blk_mask == 0).sum()
        return unmasked_bins/on_detector

    def correct_the_variance(self):
        self.unmasked_frac = self.what_fraction_unmasked()
        self.skyvar = self.compute_variance(self.sci)
        self.masked_mean_var_rdnoise = self.masked_mean(self.var_rdnoise)
        self.correction_factor = self.skyvar / self.masked_mean_var_rdnoise
        self.predicted_skyvar = self.correction_factor * self.var_rdnoise
        print(f"{self.fitsfile}")
        print(f"Robust masked mean VAR_RDNOISE: {self.masked_mean_var_rdnoise}")
        print(f"Robust masked mean SKY_VARIANCE: {self.skyvar}")
        print(f"Correction factor: {self.correction_factor}")
        print(f"Fraction of pixels unmasked: {self.unmasked_frac}")
        
    def write_file(self):
        prefix = self.fitsfile[:self.fitsfile.rfind('_')] # replace last suffix 
        outfile = prefix+'_skycor.fits'
        self.hdu[6].data = self.predicted_skyvar
        self.hdu.writeto(path.join(self.datadir,outfile))
        self.hdu.close()
