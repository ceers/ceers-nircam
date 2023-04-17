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
mosaic_background.run
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

