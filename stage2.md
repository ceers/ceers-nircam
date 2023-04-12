# Stage 2: Individual Image Calibrations

Reduction scripts related to Stage 2 


## Table of contents

* [Reduction Steps](#summary)
* [JWST Calibration Image2 Pipeline](#image2)


<a name='summary'></a>
## Reduction Steps

As a summary, running a single rate file (`jw01345001001_02201_00001_nrca1_rate.fits`) through Stage 2 would look like:
```
strun image2_1.7.2.asdf jw01345001001_02201_00001_nrca1_rate.fits
```


<a name='image2'></a>
## JWST Calibration Image2 Pipeline

Stage 2 of the JWST Calibration Pipeline involves individual image calibrations
such as flat-fielding and flux calibrating the data. The input to Stage 2 is 
a wisp- and 1/f noise-subtracted `*rate.fits` image that is in units of 
counts/s. The end product of Stage 2 is a flux calibrated image in units of 
MJy/sr. See the [Pipeline Documentation](https://jwst-pipeline.readthedocs.io/en/latest/jwst/pipeline/calwebb_image2.html#calwebb-image2) for more 
information on Stage 2.

In our reduction, we adopt the default parameter values for the Stage 2 steps. 
We use a Pipeline parameter file (`image2_1.7.2.asdf`) to run Stage 2 using 
the command line `strun` method. In `image2_1.7.2.asdf`, we have turned off 
the Resample step (which produces a rectified, distortion-corrected image 
for quick-look use) to save time. 


