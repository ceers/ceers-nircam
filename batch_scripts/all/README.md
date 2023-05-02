# all

Reduction scripts and information related to CEERS Epoch 1 NIRCam imaging
(pointings 1, 2, 3 and 6)

## Pointing Information

There are 690 input NIRCam images in CEERS Epoch 1. 

The WCS paramters for the drizzled mosaic are:

* CRVAL: (214.825,52.825)  
* CRPIX: (31585.5,9044.5)  (Note that these are the CRPIX listed in the mosaic
  headers. They are `+1` of the CRPIX in `image3_fullceers.asdf` to account for 
  Python's indexing)  
* Rotation: -49.7 degrees 
* Pixel Scale: 0.03"/pixel
* Pixfrac: 1.0
* Array Shape: (46200,15400)


## Reduction Steps

To fully reduce the 690 images in CEERS Epoch 1, run the following scripts
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
python mosaic_background.py epoch1 --add_hst
```

The scripts in this directory are set up to create combined mosaics 
of all four pointings. These mosaics can require a very significant amount of 
memory. We did not release these Epoch1 mosaics in DR0.5, but provide the 
scripts here in case you wish to create the full mosaics. Note that if you 
want to include HST images in the merged source mask for the mosaic background
subtraction step (`add_hst`), then you will need to create cutouts of the 
HST mosaics that are available at [CEERS HDR1](https://ceers.github.io/hdr1.html).

## File List

In addition to the 9 reduction scripts listed above, the following files 
are needed to reduce CEERS Epoch 1:

* `image3_fullceers.asdf` - Stage3 parameter file setup for the Resample step
  with parameters for creating the Epoch 1 mosaic on the full CEERS footprint
* `nircam[pointing]_[filter].json` - Association files of all images for the 
  given filter, listing images with the `*tweakreg.fits` suffix for use with 
  `image3_part1.asdf`. We suggest running this step pointing by pointing
  rather than for the entire Epoch 1 at once to save memory, especially in 
  the outlier detection step, and have therefore provided the association 
  files separately for each pointing. We note that there are some small 
  overlaps with neighboring pointings at the edges of images, and the 
  neighboring pointings will not be included in the outlier detection. 
* `epoch1_[filter]_final.json` - Association files of all images for the 
  given filter, listing images with the `*a3001_match.fits` suffix for use
  with `image3_fullceers.asdf`

