# Stage 3: Ensemble Processing

Reduction scripts related to Stage 1 and our custom routines


## Table of contents

* [Reduction Steps](#summary)
* [JWST Calibration Image3 Pipeline](#image3)
* [Astrometric Alignment](#tweakreg)
* [Outlier Detection](#outliers)
* [Sky Subtraction, Variance Map Updates](#skysubvar)
* [Mosaic Creation](#mosaics)
* [Background Subtraction](#bkgsub)


<a name='summary'></a>
## Reduction Steps

The Stage 3, ensemble processing is run on a set of images, for example all
images from a single filter.

As a summary, running the ensemble processing steps on all F150W images in 
CEERS NIRCam pointing 1 would look like:
```
python run_tweakreg.py jw01345001001 f150w 0
strun image3_part1.asdf nircam1_f150w.json
python skysub_wcs_varscale.py --all_images
strun image3_nircam1.asdf nircam1_f150w_final.json
```

Each of these scripts is described in detail below. 


<a name='image3'></a>
## JWST Calibration Image3 Pipeline

Stage 3 of the JWST Calibration Pipeline performs ensemble reduction steps, 
where the output product is a single mosaic per filter combining all images 
and dithers. These steps include astrometric alignment, background matching, 
outlier detection, and resampling the images onto a common output grid. We 
break this stage up into individual steps, each described below. 
See the [Pipeline Documentation](https://jwst-pipeline.readthedocs.io/en/latest/jwst/pipeline/calwebb_image3.html) for more information on Stage 3.

The Stage 3 steps are run on a set of images using association files. We have
created association files for every CEERS NIRCam pointing and filter. They 
are available in the `batch_scripts` directory. We also provide the script
used to create these files (`prep_stage3.py`), though there is no need to 
run this script if you use the already prepared files in `batch_scripts`.


<a name='tweakreg'></a>
## Astrometric Alignment

The Astrometric alignment is described in Section 3.3.1 of [Bagley et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJ...946L..12B/abstract).

We perform an astrometric calibration using a modified version of the 
Pipeline TweakReg routine, with the modifications primarily aimed at 
exposing more fitting parameters and being able to input our own catalogs. 
The main addition to the default TweakReg routine is the ability to 
specify a user-provided catalog to use as the absolute astrometric reference
as well as user-provided individual source catalogs for each image.
These modifications are required for Pipeline versions <1.8, and so must 
be used to fully recreate our DR0.5 reduction. These modifications were
incorporated into the Pipeline starting with version 1.8.

To use our modified Tweakreg, you first replace the installed Pipeline
routine with our version of `tweakreg_step.py`. First, identify where the 
Pipeline is installed:
```
from jwst.tweakreg import tweakreg_step 
print(tweakreg_step.__file__)
```
This will print the absolute path to the Pipeline version of 
`tweakreg_step.py`. Replace the Pipeline script with the `tweakreg_step.py` 
from this GitHub repo.

Next, we provide a wrapper to call the TweakReg step, `run_tweakreg.py`. 
The wrapper runs Source Extractor on each individual input image to 
detect the source positions using the windowed centroid coordinates, and 
our modified `tweakreg_step.py` matches these input positions with the 
sources in the absolute reference catalog.

TweakReg should be run on each filter and visit separately. To run it on 
F115W images in the first visit of observation 1:
```
python run_tweakreg.py jw01345001001 f115w 
```
The wrapper will look for all F115W images with the prefix `jw01345001001` 
that are in the current working directory. It will group them into an 
association file, run Source Extractor, prepare the input catalogs, invoke
TweakReg, and parse the TweakReg output into an output file called 
`tweakreg_results_all.jw01345001001.txt`. This output file will summarize
the X,Y shifts, rotations and scalings of the relative and absolute 
astrometric fits as well as the RMS of the alignments. 

The parameters for the fits are set using configuration files, for example
`f115w.cfg`. 

```
python run_tweakreg.py jw01345001001 f115w --save_fit
```



<a name='outliers'></a>
## Outlier Detection

The outlier detection step of the Pipeline attempts to identify bad pixels or 
cosmic rays that were not detected during the jump step from Stage 1. This 
step identifies outliers from a median stack of all available images and 
flags them in the DQ array of each input image. We adopt all of the default 
parameter values for the Pipeline v1.7.2, noting that in future reduction 
versions we will improve on this outlier detection. 

We perform this step with the `strun` command and the Pipeline parameter 
file `image3_part1.asdf`. This parameter file runs the SkyMatch and 
OutlierDetection steps of Stage 3 with the v1.7.2 pipeline defaults. We 
include the SkyMatch step so that the headers are appropriately updated with 
all the necessary keywords, though the sky background is not subtracted at 
this time (see [Sky Subtraction](#skysubvar) below).

All Stage 3 processing works on sets of images, using association files 
that list all images for processing

Running the step on all F150W images from CEERS NIRCam pointing 1:

```
strun image3_part1.asdf nircam1_f150w.json
```

This will output files called `*crf.fits` ('cosmic ray flagged') images in 
the `calibrated` directory. 


**Customization Options:**

Within `image3_part1.asdf`:
* Input/output: The input and output directory locations can be changed via
  the `input_dir` and `output_dir` parameters. The association files we
  provide in `batch_scripts` assume input images can be found in a directory
  called `calibrated`.


<a name='skysubvar'></a>
## Sky Subtraction, Variance Map Updates

We perform three additional corrections to the individual exposures, which are 
described in Section 3.3.2 of [Bagley et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJ...946L..12B/abstract).

First, rescale the variance maps to include the sky RMS fluctuations. we run 
a background subtraction routine that creates a 2D model of the background 
in an image. This saves a 2D background-subtracted version of the cal file 
but does not alter the original cal file. We only want to subtract a pedestal 
from the images before stacking, the 2D background subtraction routine can 
oversubtract in regions around extended/low surface brightness sources. 

Next we subtract a pedestal value in MJy/sr from each exposure, necessary
because the SkyMatch routine of jwst does not successfully match the background
across all detectors. The tiered source masks created during the 
[1/f noise step](stage1.md) are used to mask source flux. We then fit a 
Gaussian to the distribution of unmasked, sigma-clipped pixel fluxes and 
subtract the peak of the Gaussian from the image. 

Finally, we fill 'holes' in the variance maps. In v1.7.2 and previous versions 
of the jwst Pipeline, we found that known bad pixels had values of exactly 
zero in the variance arrays. In these areas, the input error array with the 
missing data or bad pixel did not contribute to the rms of the affected output
pixel during drizzling. We set these pixels to infinity in the individual 
variance maps to correctly down-weight bad pixels during drizzling.

All three of these steps are performed with the script `skywcsvar.py`. If
necessary, this script will also update the WCS from saved Tweakred `asdf` 
files. In this way, you can create multiple reduction versions without needing
to rerun Tweakreg. If the Tweakreg step is not identified as 'complete' in the
image header, `skywcsvar.py` will look for the correspdonding `*_tweakreg.asdf`
and update the datamodel accordinly.

To run `skywcsvar.py` on a single image output by `image3_part1.asdf`:
```
python skywcsvar.py --image jw01345001001_10201_00001_nrca1_a3001_crf.fits
```

Alternatively, you can run it on all images with a given filename suffix by 
setting the global variable FILE_SUFFIX and running it as:
```
python skywcsvar.py --all_images
```

**Customization Options:**

At the top of `skywcsvar.py`:
* Input/output: provide the relative (or absolute) paths to the directory
  containing the input calibrated files and the directory for the output,
  corrected images. Default is `calibrated` for both.
* Provide the directory containing the saved Tweakreg asdf files if the 
  WCS in the input images has not already been updated/tweaked
* Provide the file suffixes for the input images (default = 'crf'), the
  output 2D background subtracted images ('bkgsub1'), and the output 
  corrected images ('match')


<a name='mosaics'></a>
## Mosaic Creation

Mosaic creation is described in Section 3.3.2 of [Bagley et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJ...946L..12B/abstract).

We create pixel-aligned mosaics for each CEERS NIRCam pointing and filter
using the Resample step of Stage 3. We drizzle all images onto a 30mas pixel
scale using a pixfrac of 1. The images are pixel-aligned by specifying the 
same output CRPIX, CRVAL, rotation, pixel scale, shape and size when drizzling.
These values are set in the Pipeline parameter file 
`image3_nircam[pointing].asdf`, where `[pointing]` refers to the CEERS NIRCam
pointing ID. These are all calculated to have the same tangent point as the
HST mosaics available on our [website](http://ceers.github.io/hdr1.html).
While this is technically a parameter file for all of Stage 3, all other 
steps are skipped. You may wish to run the Source Catalog step as well to 
play with source detection parameters, though note that the values currently
in the parameter file are unoptimized defaults.

To create a F150W mosaic in CEERS NIRCam pointing 1:
```
strun image3_nircam1.asdf nircam1_f150w_final.json
```

The `nircam[pointing]_[filter]_final.json` association files list the same 
images as `nircam[pointing]_[filter].json`, but with the updated filename
suffixes that include 

We note that in Bagley et al. (2023), we mentioned using two different 
parameter files for creating SW and LW mosaics, where the only difference was
the input-to-output pixel size ratio. This was a holdover from an earlier
version of the pipeline where that distinction was necessary. The pixel 
size ratio parameter is now ignored if the pixel scale is set, and so we 
have only a single parameter file for all filters.

**Customization Options:**

Within `image3_nircam[pointing].asdf`:
* Input/output: The input and output directory locations can be changed via
  the `input_dir` and `output_dir` parameters. The association files we
  provide in `batch_scripts` assume input images can be found in a directory
  called `calibrate`.


<a name='bkgsub'></a>
## Background Subtraction

Mosaic creation is described in Section 3.3.3 of [Bagley et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJ...946L..12B/abstract).

Finally, we estimate and subtract any remaining background in the mosaics 
using a custom Python script that efficiently masks source flux before fitting 
the unmasked pixels with a two-dimensional model. The source masking is 
performed in several tiers to detect and mask progressively smaller sources.
Once a source mask is created for each individual mosaic, the masks from all
available filters in the pointing are combined to create one final, merged
mask. The background in each filter is then estimated using the unmasked 
pixels. The result of this step is a background-subtracted image, a 
background map, and a "tiermask," where anything set to 0 is used for 
estimating background and anything non-zero is not used.

**Steps in the process:**

 - Mask sources in each individual image using a tiered approach
 - Merge the masks from each filter for the pointing
 - In some cases "hand edit" the masks to remove spurious sources near the 
   borders or patches of scattered light that were identified as astronomical 
   sources to be masked.
 - After this editing, merge all of the masks from HST and JWST (with the 
   possible excception of bands where the mask was not very useful -- e.g. 
   F105W for some fields where there is not much area overlap).
 - Re-measure the background for each image using this merged mask

We perform this background subtraction with `mosaic_background.py`, which 
will do all of the above steps on all imaging for a single CEERS NIRCam
pointing. To run `mosaic_background.py` on the mosaics from NIRCam pointing 1 
(`nircam1`):
```
python mosaic_background.py nircam1
```

The parameters we used to generate the tiermasks for each filter are 
provided in `mosaic_background.cfg`. The `mosaic_background.py` routine 
uses this config file when running the initial background subtraction and 
source masking on each individual filter, and so updating these parameters
in the config file will change them for `mosaic_background.py`.

To include HST imaging:
```
python mosaic_background.py nircam1 --add_hst
```
This will include the source masks from the HST images in the merged mask. 
The HST images will also be background subtracted, though those images already
have backgrounds close to zero. The HST imaging is available for each CEERS 
NIRCam pointing on the CEERS website: 
[ceers.github.io/dr05.html](https://ceers.github.io/dr05.html). 
These cutouts are taken from the full HST mosaics available at 
[ceers.github.io/hdr1.html](https://ceers.github.io/hdr1.html).


**Customization Options:**

At the top of `mosaic_background.py`:
* Input/output: provide the relative (or absolute) path to the directory
  containing the mosaics. The defauls is `calibrated`, and the 
  background-subtracted mosaics and merged mask will be saved to the same
  directory.
* Provide the directory containing the HST imaging. This is necessary if 
  the HST imaging will be included in the merged mask.
* Provide the filename for the merged source mask.
* Provide the file suffixes for the output 2D background subtracted images
  based on individual source masks (default = 'bkgsub1'), and the output 2D 
  background subtracted images based on the merged mask ('mbkgsub1')

