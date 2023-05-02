# nircam3

Reduction scripts and information related to CEERS NIRCam pointing 3

## Pointing Information

There are 120 input images in CEERS NIRCam3. 

The WCS paramters for the drizzled mosaic are:

* CRVAL: (214.825,52.825)  
* CRPIX: (6589.5,-813.5)  (Note that these are the CRPIX listed in the mosaic
  headers. They are `+1` of the CRPIX in `image3_nircam3.asdf` to account for 
  Python's indexing)  
* Rotation: -49.7 degrees 
* Pixel Scale: 0.03"/pixel
* Pixfrac: 1.0
* Array Shape: (10500,4800)


## Reduction Steps

To fully reduce the 120 images in CEERS NIRCam3, run the following scripts
in this order:

```
detector1_and_snowballs.run
wisp_subtraction.run
remstriping.run
image2.run
update_wcs.run
image3_part1.run
skysub_varscale.run
make_mosaics.run
python mosaic_background.py nircam3 --add_hst
```

## File List

In addition to the 8 reduction scripts listed above, the following files 
are needed to reduce CEERS NIRCam3:

* `image3_nircam3.asdf` - Stage3 parameter file set up for the Resample step
  with parameters for creating the NIRCam3 mosaic
* `nircam3_[filter].json` - Association files of all images for the given 
  filter, listing images with the `*tweakreg.fits` suffix for use with 
  `image3_part1.asdf`
* `nircam3_[filter]_final.json` - Association files of all images for the 
  given filter, listing images with the `*a3001_match.fits` suffix for use
  with `image3_nircam3.asdf`


## Note on TweakReg

We encountered issues aligning NIRCam3 to the HST F160W mosaic, specifically 
for a ~1' region in the quadrant of the mosaic mostly covered by the B2 
detector where astrometric offsets from NIRCam-to-HST are large (~0.05"). 
One possible explanation is that the HST F160W exposures for the particular 
visit in that part of the mosaic were affected by a low-level guide star 
tracking issue. 

We therefore followed a different process for aligning this NIRCam pointing. 
We first aligned F277W to F160W, excluding all sources in the upper left 
quadrant of detector BLONG from the fit. We then used F277W rather than F160W
as the reference catalog for all other filters in this pointing. See
Section 3.3.1 of [Bagley et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJ...946L..12B/abstract) for more information. 
The `*tweakreg.asdf` files in `tweakreg_wcs/` account for this change in our 
TweakReg method for NIRCam3, and we again recommend using these saved WCS
models instead of running TweakReg from scratch. 

