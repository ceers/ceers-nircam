# nircam6

Reduction scripts and information related to CEERS NIRCam pointing 6

## Pointing Information

There are 120 input images in CEERS NIRCam6. 

The WCS paramters for the drizzled mosaic are:

* CRVAL: (214.825,52.825)  
* CRPIX: (7352.5,3501.5)  (Note that these are the CRPIX listed in the mosaic
  headers. They are `+1` of the CRPIX in `image3_nircam6.asdf` to account for 
  Python's indexing)  
* Rotation: -49.7 degrees 
* Pixel Scale: 0.03"/pixel
* Pixfrac: 1.0
* Array Shape: (10500,4800)


## Reduction Steps

To fully reduce the 120 images in CEERS NIRCam6, run the following scripts
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
python mosaic_background.py nircam6 --add_hst
```

## File List

In addition to the 8 reduction scripts listed above, the following files 
are needed to reduce CEERS NIRCam6:

* `image3_nircam6.asdf` - Stage3 parameter file set up for the Resample step
  with parameters for creating the NIRCam6 mosaic
* `nircam6_[filter].json` - Association files of all images for the given 
  filter, listing images with the `*tweakreg.fits` suffix for use with 
  `image3_part1.asdf`
* `nircam6_[filter]_final.json` - Association files of all images for the 
  given filter, listing images with the `*a3001_match.fits` suffix for use
  with `image3_nircam6.asdf`

