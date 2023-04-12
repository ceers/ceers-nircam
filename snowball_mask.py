#!/usr/bin/env python
# Create a mask for snowballs from a DQ file

__author__ = "Henry C. Ferguson, STScI"
__version__ = "0.2.0"
__license__ = "BSD3"

import numpy as np
from astropy.io import fits
from scipy import ndimage
from astropy.convolution import Tophat2DKernel
from scipy.ndimage import binary_dilation, median_filter

from jwst.datamodels import dqflags


# Snowball masking algorithm

def snowball_tier_mask(dqarray,detection_kernel,growing_kernel):
    # Find big groups of jumps
    jumps = dqarray == dqflags.pixel['JUMP_DET']
    filtered_jumps = median_filter(jumps, footprint=detection_kernel)
    # Incorporated saturated pixels that are encompassed within these big groups
    saturated = dqarray == dqflags.pixel['SATURATED']
    surrounded_saturated = saturated & filtered_jumps
    new_jumps = jumps | surrounded_saturated
    # Find the big groups again, having incorporated these saturated pixels
    filtered_jumps = median_filter(new_jumps, footprint=detection_kernel)
    # Only need to grow the mask if there is something in it...
    if (filtered_jumps.sum() > 0):
        grown_snowball = ndimage.binary_dilation(filtered_jumps,growing_kernel)
    else:
        #print("No snowballs")
        grown_snowball = filtered_jumps
    snowball_mask = grown_snowball * dqflags.pixel['JUMP_DET']
    return snowball_mask

def create_kernels(detection_radii, grow_radii):
    detection_kernels = []
    growing_kernels = []
    for dr,gr in zip(detection_radii,grow_radii):
        detection_kernels += [Tophat2DKernel(dr).array]
        growing_kernels   += [Tophat2DKernel(gr).array]
    return detection_kernels, growing_kernels

def snowball_mask_groupdq(groupdq, detection_radii=[7,15], grow_radii=[7,35]):
    # Set up for looping through the groupdq planes
    n_integrations = groupdq.shape[0]
    n_groups = groupdq.shape[1]
    new_groupdq = groupdq.copy()
    snowball_mask_4d = np.zeros(groupdq.shape,groupdq.dtype)
    detection_kernels, growing_kernels = create_kernels(detection_radii,grow_radii)
    # Step through the integrations and groups
    for integration in range(n_integrations):
        for group in range(n_groups):
            dqarray = groupdq[integration,group,:,:]
            plane_mask = snowball_mask_single_plane(dqarray,detection_kernels,growing_kernels)
            snowball_mask_4d[integration,group,:,:] = plane_mask
            new_groupdq[integration,group,:,:] = np.bitwise_or(dqarray,plane_mask)
    return snowball_mask_4d, new_groupdq

def snowball_mask_single_plane(dqarray,detection_kernels,growing_kernels):
    #print("starting mask_single")
    snowball_mask = np.zeros(dqarray.shape,bool)
    for dk,gk in zip(detection_kernels, growing_kernels):
        snowball_mask = snowball_mask | snowball_tier_mask(dqarray,dk,gk)
    return snowball_mask

