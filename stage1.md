# Stage 1: Detector-level Corrections

Reduction scripts related to Stage 1 and our custom routines


## Table of contents

* [Reduction Steps](#summary)
* [JWST Calibration Detector1 Pipeline](#detector1)
* [Snowball Correction](#snowballs)
* [Wisp Subtraction](#wisps)
* [1/f Noise Subtraction](#1fnoise)


<a name='summary'></a>
## Reduction Steps

As a summary, running a single file (`jw01345001001_10201_00001_nrca3_uncal.fits`) through all of the detector-level corrections (from the Pipeline and our custom scripts) would look like:
```
python snowball_wrapper.py jw01345001001_10201_00001_nrca3_uncal.fits
python wispsub.py jw01345001001_10201_00001_nrca3_rate.fits --fit_scaling
python remstriping.py jw01345001001_10201_00001_nrca3_rate.fits --save_patterns
```
Each of these scripts is described in detail below. 


<a name='detector1'></a>
## JWST Calibration Detector1 Pipeline

Stage 1 of the JWST Calibration Pipeline performs detector-level corrections,
many of which are common to all instruments and observing modes. The reduction
steps involve initializing the data quality (DQ) arrays for flagging pixels,
identifying saturated pixels, subtracting the superbias, using reference
pixels to correct for readout noise, correcting pixels for nonlinearity,
subtracting the dark current, identifying cosmic rays as jumps in each pixel's
up-the-ramp signal, and calculating a linear fit to the unflagged ramp data
to determine the average count rate per pixel. The end product of Stage 1 is
a count-rate image in units of counts/s. See the [Pipeline Documentation](https://jwst-pipeline.readthedocs.io/en/latest/jwst/pipeline/calwebb_detector1.html#calwebb-detector1) for more information on Stage 1.

In our reduction, we adopt the default parameter values for these steps. We
run the reduction steps of Stage 1 together with the
[Snowball Correction](#snowballs) (described below), using the Python script
`snowball_wrapper.py`. For reference we also provide the pipeline parameter
file `detector1_1.7.2.asdf` with the equivalent set of default parameter
values.

<a name='snowballs'></a>
## Snowball Correction

The snowball correction is described in Section 3.1.1 of [Bagley et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJ...946L..12B/abstract).

We perform this correction with our Python script `snowball_wrapper.py`, which 
runs Stage 1 of the Pipeline, saves the ramps from the first run of the 
ramp-fitting step, identifies snowballs and increases their footprint, flags 
them as cosmic rays, and then runs the ramp-fitting step again to create a 
count-rate map that excludes the flagged portions of the ramps. 

Running the step on a single image:
```
python snowball_wrapper.py jw01345001001_10201_00001_nrca3_uncal.fits
```

**Customization Options:**

At the top of `snowball_wrapper.py`: 
* Input/output: provide the relative (or absolute) paths to the directory 
  containing the raw `uncal` files and the directory for the Detector1 outputs.
  Defaults are `uncals` and `calibrated`. 
* Multi-processing: Set `MAXCORES` to use multi-processing during the Ramp
  Fitting step. Options are `'none'` (for no multi-processing), `'quarter'`
  (to use 1/4 of available cores), `'half'`, or `'all'`
* Deleting interim files: There are 3 interim files that are produced during
  the Stage 1 Pipeline. Set `DELINTERIM` to delete these files to save disk 
  space. 


<a name='wisps'></a>
## Wisp Subtraction

Wisp subtraction is described in Section 3.1.2 of [Bagley et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJ...946L..12B/abstract).

We subtract the "wisp" features from detectors A3, A4, B3 and B4 of the F150W 
and F200W images. We use the wisp templates provided by the NIRCam team that 
were released on 2022 August 26. (Wisps are described in the [JDox](https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-features-and-caveats/nircam-claws-and-wisps), and the templates are [hosted at STScI](https://stsci.app.box.com/s/1bymvf1lkrqbdn9rnkluzqk30e8o2bne) and called wisps_2022_08_26.tgz.).
The wisp features can have a variable brightness from exposure to exposure, 
and so our subtraction routine scales the templates to match the brightness 
of the feature in each individual image. 

Running the step on a single image:

```
python wispsub.py jw01345001001_10201_00001_nrca3_rate.fits --fit_scaling 
```

The `--fit_scaling` option will attempt to determine the best scale factor
that matches the strength of the feature in the given image, by minimizing 
the variance of (image - scale * template).

However, we find that for some images, it is necessary to override the 
automatically-determined scaling factor, especially when bright sources, 
artifacts, etc., cause the scaling function to over or under fit the wisp
feature. 

To run the step on a single image when supplying a manual scale factor:
```
python wispsub.py jw01345001001_10201_00001_nrca3_rate.fits --scale 0.7
```

We inspected all of the CEERS wisp-subtracted rate files by eye and made a 
few manual adjustments. We provide the scaling factor for each individual 
image in the `batch_scripts` directory for those who aim to reproduce our
reduction.

**Note:** This version of `wispsub.py` uses Photutils to identify sources 
in the images before fitting the wisp feature. It has been tested successfully
with Photutils v1.5.0, but may break for newer versions of Photutils, in 
which the `make_source_mask` function was moved to a different part of the 
module.


**Customization Options:** 

At the top of `wispsub.py`:
* Input/output: provide the relative (or absolute) paths to the directory
  containing the input `rate` files and the directory for the wisp-subtracted 
  outputs. Default is `calibrated` for both
* Also, provide the relative (or absolute) path to where you are storing
  the wisp templates.


<a name='1fnoise'></a>
## 1/f Noise Subtraction

1/f noise subtraction is described in Section 3.1.3 of [Bagley et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJ...946L..12B/abstract).

We measure and subtract the horizontal and vertical striping patterns in 
the images introduced during detector readout with the Python script 
`remstriping.py`, which applies the flat field, subtracts a background 
pedestal, masks source flux, and then measures a sigma-clipped median along 
first rows and then columns. For the row correction (horizontal striping), we 
measure and remove the pattern individually for each of the four amplifiers. 
In some images, especially in the SW filters, the difference from amplifier 
to amplifier can vary by ~3%-5%, and so this amplifier-dependent correction 
works better than using the median from the full row. However, in some cases,
especially around bright or extended sources, a large number of pixels in a 
given amp-row are masked, leaving too few to calculate a robust median. In 
these cases, we use the median of the entire row for the affected amp-row. 

The threshold used to determine the cutoff in number of masked pixels per
amp-row varies from image to image to allow for the best pattern subtraction. 
The script `remstriping.py` uses a threshold of `0.8` for all images, meaning
if more than 80% of the pixels in an amp-row are masked, the full-row median
is used for the amp-row. The 80% value is set at the top of the script via 
the global variable `MASKTHRESH`.

To run the 1/f noise subtraction step on a single image:
```
python remstriping.py jw01345001001_10201_00001_nrca3_rate.fits
```

This script creates several outputs:
* `*rate_1fmask.fits` - The segmentation map used to mask sources. This mask
  will be used again during Stage 3 processing
* `*rate_horiz.fits`, `*rate_vert.fits` - images saving the horizontal and 
  vertical striping patterns that were subtracted from the image. These two
  are not saved by default. 

To save the striping patterns for later reference:
```
python remstriping.py jw01345001001_10201_00001_nrca3_rate.fits --save_patterns
```

For DR0.5, we visually inspected each 1/f noise-subtracted image and adjusted
the masking threshold for the images with over- or under-subtracted pattern 
noise. For example, images with large, bright sources usually benefit from 
a threshold closer to 0.1 than 0.8. The threshold can either be changed 
using the global variable `MASKTHRESH`, or by setting a value manually for 
each image:

```
python remstriping.py jw01345001001_10201_00001_nrca3_rate.fits --thresh 0.15
```

We provide the threshold used for each individual image in the `batch_scripts`
directory for those who aim to reproduce our reduction.

**Customization Options:**

At the top of `remstriping.py`:
* Input/output: provide the relative (or absolute) paths to the directory
  containing the input `rate` files and the directory for the 1/f-subtracted
  outputs. Default is `calibrated` for both
* Also, change the masking threshold (fraction of masked pixels in an amp-row)
  for switching from using the amp-row median to a full-row median with the 
  global variable `MASKTHRESH`. 

