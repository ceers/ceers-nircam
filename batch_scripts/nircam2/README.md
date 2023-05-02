# nircam2

Reduction scripts and information related to CEERS NIRCam pointing 2

## Pointing Information

There are 240 input images in CEERS NIRCam2. 

The WCS paramters for the drizzled mosaic are:

* CRVAL: (214.825,52.825)
* CRPIX: (16681.5,923.5)  (Note that these are the CRPIX listed in the mosaic
  headers. They are `+1` of the CRPIX in `image3_nircam2.asdf` to account for
  Python's indexing)
* Rotation: -49.7 degrees 
* Pixel Scale: 0.03"/pixel
* Pixfrac: 1.0
* Array Shape: (11000,6450)


## Reduction Steps

To fully reduce the 240 images in CEERS NIRCam2, run the following scripts
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
python mosaic_background.py nircam2 --add_hst
```

Note that CEERS NIRCam2 has additional imaging in F200W and F444W 
(observation 52) that was obtained one week later and was heavily affected by 
persistence. This additional imaging is discussed in Section 3.1.4 of 
[Bagley et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJ...946L..12B/abstract). The reduction scripts for NIRCam2 include this extra imaging. For 
the mosaicking step, we provide the option to create mosaics with and 
without the extra imaging.

## File List

In addition to the 8 reduction scripts listed above, the following files 
are needed to reduce CEERS NIRCam2:

* `image3_nircam2.asdf` - Stage3 parameter file set up for the Resample step
  with parameters for creating the NIRCam2 mosaic
* `nircam2_[filter].json` - Association files of all images for the given 
  filter, listing images with the `*tweakreg.fits` suffix for use with 
  `image3_part1.asdf`
* `nircam2_[filter]_final.json` - Association files of all images for the 
  given filter, listing images with the `*a3001_match.fits` suffix for use
  with `image3_nircam2.asdf`
* `nircam2b_[filter]_final.json` - Association files of all the extra 
  imaging for the F200W and F444W filters, for use with `image3_nircam2.asdf`
  to make a mosaic of only the additional imaging
* `nircam2all_[filter]_final.json` - Association files of all F200W and F444W
  images, including the extra imaging for, for use with `image3_nircam2.asdf`
  to make a mosaic of all available imaging in these two filters


