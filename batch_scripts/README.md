<a name='batch'></a>
## Batch Processing

Here we provide scripts that will run each of the reduction steps on all 
imaging from the Epoch 1 observations. We used the Texas Advanced Computing 
Center (TACC) at the University of Texas at Austin for our image reduction, 
and so these scripts were designed to be executed in parallel using a 
distributed computing system. Almost all of these scripts therefore provide 
one command per line.

The scripts assume input directory called `uncals` (for stage 1 processing
only) and all processed outputs are saved in a directory called `calibrated`.
We describe each reduction step in the READMEs for [Stage1](../stage1.md), 
[Stage2](../stage2.md) and [Stage3](../stage3.md), and also explain how to 
change the input/output options in each reduction script.

We provide the output WCS calculated by our runs of the tweakreg step for 
each pointing, along with the Python code that will apply this tweaked WCS
to each image header. Therefore, you do not need to run the TweakReg step 
yourself, and can use the `updatewcs.py` script instead to use the CEERS
astrometry from DR0.5.

We provide scripts to reduce the CEERS NIRCam pointings individually in the
directories `nircam1`, `nircam2`, etc., as well as scripts to reduce all 
690 images at once in the `all` directory. In each case, the batch reduction 
scripts are executed in the following order:

* `detector1_and_snowballs.run` - Runs `snowball_wrapper.py` on each uncal 
  image, which runs the Detector1 Pipeline up to the ramp fitting step, 
  identifies and flags snowballs, and then reruns the ramp fitting step with 
  the updated DQ arrays
* `wisps_subtraction.run` - Runs `wispsub.py` on each rate image, using 
  the scaling factors we determined for DR0.5
* `remstriping.run` - Runs `remstriping.py` on each rate image, using the 
  masking thresholds we determined for DR0.5
* `image2.run` - Runs the Image2 Pipeline on each rate image
* `update_wcs.run` - Updates the WCS in each image header using the output 
  saved in `tweakreg_wcs` from our run of TweakReg
* `image3_part1.run` - Runs the outlier detection step of the Image3 Pipeline
* `skysub_varscale.run` - Runs `skywcsvar.py` on each image, which measures
  and removes a sky pedestal, scales the readnoise variance maps to include 
  an estimate of the sky RMS and fills holes in the variance maps
* `make_mosaics.run` - Creates a mosaic for each pointing and filter 
* `mosaic_background.run` - Runs the 2D background subtraction on the mosaic 
  images

Note that the Stage 3 scripts in the `all` directory (the `*final.json` 
association files and `image3_full.asdf` parameter file) are set up to create 
combined mosaics of all four pointings. These mosaics can require a 
very significant amount of memory. 
