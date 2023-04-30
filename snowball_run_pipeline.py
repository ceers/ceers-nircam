#!/usr/bin/env python
# Run the calwebb detector1 pipeline once to get the ramp.fits file.
# Use that to mask snowballs.
# Pick up the pipeline in the ramp-fitting step, using this new ramp file.

__author__ = "Henry C. Ferguson, STScI"
__version__ = "0.2.3"
__license__ = "BSD3"

# Version history
# 0.2.1 -- fixed bug where I was appending '_uncal.fits' twice 
# 0.2.2 -- 20221220 M Bagley edit to pass maximum_cores to ramp fitting step
# 0.2.3 -- 20230324 M Bagley edit to use run method in 
#          run_pipelin_to_create_ramp() for consistency. Now saving the 
#          rampints output only for exposures with multiple integrations. 
#          Optional ramp fitting output not saved by default. Switching 
#          numbered extensions in HDU files to keywords, i.e., hdu['GROUPDQ'] 
#          instead of hdu[3]

import numpy as np
from astropy.io import fits
import os

from snowball_mask import snowball_mask_groupdq

# Pipeline imports
from jwst import datamodels
from jwst.pipeline import Detector1Pipeline
from jwst.ramp_fitting import RampFitStep


def run_pipeline_to_create_ramp(dataset,output_dir="."):
    # Instantiate the pipeline step 1
    d1p = Detector1Pipeline()

    # Configure to skip the last two steps
    d1p.ramp_fit.skip = True
    d1p.gain_scale.skip = True
    d1p.save_calibrated_ramp = True
    d1p.save_results = True
    d1p.output_dir = output_dir

    result = d1p.run(f"{dataset}_uncal.fits")

    return result


def make_snowball_mask(rampfile):
    # Read in the ramp file
    hdu = fits.open(rampfile)
    groupdq = hdu['GROUPDQ'].data
    # Make the snowball mask
    snowball_mask, new_groupdq = snowball_mask_groupdq(groupdq)
    # Substitute this new groupdq file and write the ramp file back out
    hdu['GROUPDQ'].data = new_groupdq
    hdu.writeto(rampfile,overwrite=True)

    hdu.close()


def run_ramp_fitting(datadir,dataset,rampfile,maxcores,output_dir="."):
    # Read in the ramp as a datamodel
    jump = datamodels.RampModel(rampfile)

    # Process using the run() method
    ramp_fit_step = RampFitStep()
    ramp_fit_step.output_dir = output_dir
    ramp_fit_step.save_results = True
    ramp_fit_step.maximum_cores = maxcores

    # Set to True to save optional ramp fit info, good for visualization
    ramp_fit_step.save_opt = False

    # Run the step
    ramp_fit = ramp_fit_step.run(jump)

    return ramp_fit


def detector1_with_snowball_correction(dataset,input_dir=".",output_dir=".",maxcores='none'):
    # Run the pipeline the first time
    input_file = os.path.join(input_dir,f"{dataset}")
    result = run_pipeline_to_create_ramp(input_file,output_dir=output_dir)

    # Revise the ramp file to put in the snowball mask
    rampfile = os.path.join(output_dir,f"{dataset}_ramp.fits")
    make_snowball_mask(rampfile)

    # Re run the ramp fitting step with the new ramp file
    ramp = run_ramp_fitting(input_dir,dataset,rampfile,maxcores,
                            output_dir=output_dir)

