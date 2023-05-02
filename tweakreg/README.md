# TweakReg

Reduction scripts related to running TweakReg for astrometric alignment

<a name='top'></a>
## Table of contents

* [Astrometric Alignment](#tweakreg)
  * [Updating the WCS in Headers](#updatewcs)
  * [Running TweakReg in Pipeline v1.8+](#newtweakreg)
  * [Running TweakReg in Pipeline v<1.8](#oldtweakreg)


<a name='tweakreg'></a>
## Astrometric Alignment

The Astrometric alignment is described in Section 3.3.1 of [Bagley et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJ...946L..12B/abstract).

We perform an astrometric calibration using a modified version of the 
Pipeline TweakReg routine. Our main modification to the default TweakReg
routine is the ability to specify user-provided catalogs. In v1.7, the 
Pipeline had a limited number of absolute reference catalogs available,
none of which worked for CEERS as they had too few sources. Additionally, 
we found that the internal source detection routine was detecting too many 
spurious sources in the input images, and the source detection parameters 
were as tunable as needed. Our modified version therefore accepts user-provided
source catalogs for each input image as well as an absolute reference catalog.

Our modifications are required for Pipeline versions <1.8, and so must 
be used to fully recreate our DR0.5 reduction from scratch. To use our 
modified TweakReg, you must replace the installed Pipeline routine with our 
version of `tweakreg_step.py` that is available in this directory. However, 
these modifications were incorporated into the Pipeline starting with 
version 1.8, and so it may not be worth setting up our modified version. 
We provide an alternative option in the [following subsection](#updatewcs).

<a name='updatewcs'></a>
### Updating the WCS in Headers
 
Instead of replacing the Pipeline routine with our modified version, we 
recommend using the tweaked WCS models that we have saved from our run of 
TweakReg to update the WCS in the image headers. This will allow you to 
recreate our DR0.5 reduction without needing to alter your Pipeline 
installation. We provide the tweaked WCS models for all 690 images
in the `batch_scripts/tweakreg_wcs` folder.
To use our saved WCS model to update the header of a single cal
file (`jw01345001001_02201_00001_nrca1_cal.fits`):
```
python updatewcs.py --image jw01345001001_02201_00001_nrca1_cal.fits
```
To update the headers of all images in the specified `INPUTDIR`:
```
python updatewcs.py --all_images
```

Whether you choose to apply our tweaked WCS models for each image using
`updatewcs.py` or run the routine yourself [as described in the following
subsections](#newtweakreg), you must align the images before running the 
outlier detection routine described in the [next section](#outliers).


**Customization Options:**

At the top of `updatewcs.py`:
* Input/output: provide the relative (or absolute) paths to the directory
  containing the input calibrated files and the directory for the output,
  corrected images. Default is `calibrated` for both.
* Provide the path to the directory containing the saved Tweakreg asdf files.
  Default is `tweakreg_wcs`, which is packaged in the `batch_scripts` directory 
  of this repo
* Provide the file suffixes for the input images (default = 'cal') and the
  the output images with updated WCS headers (default = 'tweakreg')


<a name='newtweakreg'></a>
### Running TweakReg in Pipeline v1.8+

We provide a wrapper for TweakReg (`run_tweakreg.py`) that works with 
**Pipeline versions 1.8+**. You may find it easier to run the wrapper in a 
separate directory, as it creates many output files in both the current working 
directory and the input directory containing the images to be tweaked. 
We suggest moving or copying the `*cal.fits` files output by Stage 2 to the 
`tweakreg` directory or a subdirectory separate from the `calibrated` directory
that contains the majority of the processed products.

TweakReg should be run on each filter and visit separately. To run the 
wrapper on F115W images in the first visit of observation 1:
```
python run_tweakreg.py jw01345001001 f115w
```
The wrapper will look for all F115W images with the prefix `jw01345001001`
that are in the specified `INPUTDIR` directory. It will group them into
associations for each detector, run Source Extractor to detect source 
positions using the windowed centroid coordinates, prepare the input
catalogs, invoke TweakReg, and parse the TweakReg output into an output file
called `tweakreg_results_all.jw01345001001.txt`. 

TweakReg will first perform a relative alignment (e.g., all A1 detector 
images aligned with each other) and will then perform an absolute aligment 
with the provided reference catalog. For CEERS DR0.5, we used a catalog of 
sources detected in the EGS HST F160W mosaic (the mosaic is available at
[CEERS HDR1](https://ceers.github.io/hdr1.html)) as the absolute astrometric
reference: `CEERS_EGS_HST_v1.9_cat_radecmag.ecsv`.

The wrapper creates several files. For each input image:
* `*_sci.fits` - the science extension of the input cal image saved as a 
  separate file for Source Extractor
* `*_rms.fits` - the ERR map of the input cal image saved as a separate 
  file for Source Extractor
* `*_cat.idxy/ecsv` - the Source Extractor catalog reformatted for TweakReg

For each set of detector images (for example all detector A1 images from F115W):
* `jw01345003001_nrca1.json` - association files grouping all input images 
  from the specified observation, visit, and filter
* `tweakreg_jw01345003001_nrca1.log` - output log file from the relative and 
  absolute astrometric fits for this group of images

Overall output files:
* `tweakreg_results_all.jw01345003001.txt` - the overall output file 
  summarizing the relative and absolute astrometric fits for all groups of 
  images
* `tweakreg_out.log` - the TweakReg log file, that is overwritten each time 
  the wrapper calls TweakReg. This log file is copied to the detector-specific
  log files (`tweakreg_jw01345003001_nrca1.log`) that are saved for reference.
  The info for `tweakreg_results_all.jw013450*.txt` is extracted from this log 
  file. 


The output file (`tweakreg_results_all.jw013450*.txt`) will summarize the X,Y 
shifts, rotations and scalings of the relative and absolute astrometric fits 
as well as the RMS of the alignments. To quickly view the results:
```
grep Relative tweakreg_results_all.jw01345001001.txt | sort
grep Absolute tweakreg_results_all.jw01345001001.txt | sort
```

The wrapper requires a few files:

* `CEERS_EGS_HST_v1.9_cat_radecmag.ecsv` - the absolute astrometric reference
  catalog (or specify your own catalog with the `ABS_REFCAT` global variable)
* Source Extractor files - `se.config`, `se.conv`,
  `se.nnw` and `se.outputs`. These are required to run 
  the source detection on each individual input catalog, and should be 
  in the working directory.
* `tweakreg_log.cfg` - a Pipeline configuration file that sets the name of 
  the TweakReg log file. The wrapper extracts a summary of the fit results
  from the TweakReg log file. This config file should be present in the 
  working directory.
* `[filt].cfg` - a configuration file for each filter that specifies the 
  TweakReg parameters for the fit. These config files (described more below)
  should be in the input directory with the images.

The parameters for the fits are set using configuration files, for example
`f115w.cfg`. The wrapper will look for files called `filt.cfg`, where 
`filt` is specified in the call to `run_tweakreg.py`.
These config files should have one section for each visit during which 
the filter was observed. For example, F115W was observed in visit `01` for 
CEERS NIRCam1, so `f115w.cfg` would have:
```
[jw01345001001]
# source finding parameters
detect_thresh  = 4
detect_minarea = 16
deblend        = 4
# Relative alignment: 
minobj         = 20
searchrad      = 1.0
separation     = 0.1
tolerance      = 0.05     
nclip          = 1
sigma          = 0.7
# Absolute alignment:
abs_minobj     = 20
abs_searchrad  = 1.0
abs_separation = 0.2
abs_tolerance  = 0.15
abs_nclip      = 1
abs_sigma      = 0.6
abs_fitgeom    = rshift
```
If F115W had also been observed in visit `02`, there would need to be a second 
section with the ID `[jw01345001002]` so that the visits could be aligned
separately. We provide the config files for each filter and pointing used for
DR0.5 in `batch_scripts/tweakreg_cfgs`.

We recommend using a fit geometry for the absolute alignment that
allows for only shifts and rotations (`rshift`) for short wavelength filters 
and shifts, rotations and scalings (`rscale`) for the long wavelength filters. 
Our wrapper is hard-coded to use only shifts for the relative alignment.
See the documentation for the [TweakReg routine](https://jwst-pipeline.readthedocs.io/en/latest/jwst/tweakreg/README.html) for more information on each 
parameter.

You may find the need to run TweakReg several times for a set of images in 
order to determine the optimal parameters. 
When you are happy with the astrometric fit, rerun the wrapper with the
`save_results` kwarg:
```
python run_tweakreg.py jw01345001001 f115w --save_results
```
This will write the output files, each image saved with the updated WCS and
a file suffix `*tweakreg.fits`.
Additionally, it will save the WCS models (`*asdf` files) for later use so 
you can update the WCS of future versions of the images without rerunning
TweakReg (see [Updating the WCS in Headers](updatewcs) above).


**Customization Options:**

At the top of `run_tweakreg.py`:
* Input/output: provide the relative (or absolute) path to the directory
  containing the input calibrated files. Default is `calibrated`.
* Provide the path and filename of the absolute reference catalog.
* Provide the file suffixes for the input images (default = 'cal') 
* Provide the directory for saving the Tweakreg asdf files. Default 
  is `tweakreg_wcs`


<a name='oldtweakreg'></a>
### Running TweakReg in Pipeline v<1.8

To reproduce our DR0.5 reduction, we recommend using `updatewcs.py` with 
the WCS models we provide in `batch_scripts/tweakreg_wcs`. However, if you
want to run TweakReg from scratch using a Pipeline version <1.8, we provide a
wrapper `run_tweakreg_1.7.py` and modified version of the TweakReg routine 
`tweakreg_step.py` to work with it. 

To use our modified TweakReg, you must replace the installed Pipeline 
routine with our version of `tweakreg_step.py`. First, identify where the 
Pipeline is installed:
```
from jwst.tweakreg import tweakreg_step 
print(tweakreg_step.__file__)
```
This will print the absolute path to the Pipeline version of 
`tweakreg_step.py`. Replace the Pipeline script with the `tweakreg_step.py` 
from this GitHub repo. This step is necessary as the Pipeline version does 
not allow for user-specified input and absolute reference catalogs.

Next, you can run the wrapper `run_tweakreg_1.7.py` as described in the 
above section ([Running TweakReg in Pipeline v1.8+](newtweakreg)), using
the same method of config files for setting the parameters.

[Return to Top](#top)

