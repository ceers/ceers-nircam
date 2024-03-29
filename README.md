# ceers-nircam 

Reduction scripts for CEERS NIRCam imaging public data releases

## Table of contents
* [CEERS Introduction](#intro)
  * [CEERS NIRCam Imaging](#intro-imaging)
* [This Repository](#repo-purpose)
  * [Reduction Code Versions](#versions)
* [Installation and Set-up](#install)
  * [Additional packages and software](#requiredpackages)
  * [Downloading Data](#download)
* [Stage 1: Detector-level Corrections](#stage1)
* [Stage 2: Individual Image Calibrations](#stage2)
* [Stage 3: Ensemble Processing](#stage3)
* [Batch Processing](#batch)


<a name='intro'></a>
## CEERS Introduction

The Cosmic Evolution Early Release Science Survey (CEERS) has covered ~100 sq. 
arcmin of the EGS field with JWST imaging and spectroscopy using NIRCam, MIRI, 
and NIRSpec. CEERS is demonstrating, testing, and validating efficient 
extragalactic surveys with coordinated, overlapping parallel observations in a 
field supported by a rich set of HST/CANDELS multi-wavelength data. In Cycle 1,
CEERS obtained imaging and spectroscopy of the EGS field with six parallel 
JWST instrument modes: 
* Four pointings with MIRI imaging (primary) + NIRCam imaging (parallel) 
* Six pointings with NIRSpec MSA R\~1000 and four additionally with R\~100 
  (primary) + NIRCam imaging (parallel)
* Four pointings with NIRCam grism (primary) + MIRI imaging (parallel)

See our team website for more information on the program: [ceers.github.io](https://ceers.github.io)


<a name='intro-imaging'></a>
### CEERS NIRCam Imaging 

All told there are ten CEERS NIRCam pointings obtained in parallel to MIRI and 
NIRSpec observations. The tightly-constrained, overlapping CEERS observations 
were schedulable in June or December with a 180 degree position angle rotation.
Due to time constraints in the 2022 June observing window following telescope 
and instrument commissioning, the CEERS observations were split into two 
epochs. CEERS pointings 1, 2, 3 and 6 (MIRI primary + NIRCam parallel) were 
observed in Epoch 1 in June 2022. The remaining CEERS observations (six 
NIRSpec MSA+NIRCam imaging pointings and four NIRCam WFSS+MIRI imaging 
pointings) were observed in Epoch 2 in December 2022. 


<a name='repo-purpose'></a>
## This Repository

Here we provide the scripts used to reduce the NIRCam imaging that is 
available in our [public data release 0.5 (DR0.5)](https://ceers.github.io/dr05.html). 
The details of this reduction are presented in [Bagley et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJ...946L..12B/abstract). While we continue to 
update and improve our processing, here we provide the versions of the 
reduction scripts that will reproduce our public data release. Where possible 
in this repository, we also provide information on our updates and 
improvements that are not included in this version of these scripts.

We have organized this as a set of README files for each stage of the JWST 
Calibration Pipeline and any other steps that may be included. These READMEs
are called [stage1.md](stage1.md), [stage2.md](stage2.md), 
and [stage3.md](stage3.md). Each file describes the scripts and routines 
needed to perform the corresponding level of image reductions. We give examples
on how to run the scripts and how to customize them with paths, file 
suffixes, optional parameters, etc. 

The Epoch 1 observations include 690 individual images, and performing each 
reduction step on these images in series takes quite a while. We therefore 
prefer to reduce CEERS NIRCam imaging in parallel via batch scripts (see 
[Batch Processing](#batch) below), and so we have not created Jupyter 
notebooks to run these reduction steps. Please see the CEERS JWebbinar 
([JWebbinar 13](https://www.stsci.edu/jwst/science-execution/jwebbinars)) and 
associated notebook for examples on how to run each stage of the pipeline. We 
note that these notebooks were created for simulated data and use a much older 
version of the JWST Calibration Pipeline. However, they provide more detailed
information on each step of the Pipeline as well as demonstrating various 
methods for running the Pipeline. We also recommend JWebbinar3 and the
[STScI Pipeline Documentation](https://jwst-pipeline.readthedocs.io/en/latest/) 
for a more detailed explanation of each Pipeline step.

Our CEERS NIRCam reduction makes use of a combination of the Pipeline `strun` 
command and Python scripts/wrappers. The somewhat inconsistent nature of our 
coding practices is a testament to the large number of collaborators who 
contributed quickly to developing our procedures. 


<a name='versions'></a>
### Reduction Code Versions

This repository currently contains the code version we used to produce 
[DR0.5](https://ceers.github.io/dr05.html).
As time goes on, we may provide updates to our scripts here. In that case, we
will tag this version of the repository and update the table below. In this 
way, if you wish to download the version we used to produce DR0.5, you will 
be able to clone the corresponding tag. 
| ceers-nircam Tag | Corresponding CEERS Data Release | Released   | jwst Version | CRDS pmap |                                                          Notes |
|------------------|----------------------------------|------------|--------------|-----------|-------|
| --               | DR0.5                            | 2023-04-05 | v1.7.2       | 0989      | Initial code release, reduction intended to be run in parallel |


<a name='install'></a>
## Installation and Set-up

To fully reproduce our reduction, please install version 1.7.2 of the JWST 
Calibration Pipeline. We recommend doing so in a conda environment:

```
conda create -n <env_name> python=3.9
conda activate <env_name>
pip install jwst==1.7.2
```

When running a Pipeline or Pipeline step, the Pipeline will automatically
look for any required reference files in a pre-defined local directory. 
If the required reference files are not present, they will automatically be
downloaded from the [Calibration Reference Data System (CRDS)](https://jwst-crds.stsci.edu) at STScI. 

You will need to specify a local directory to store the reference files
and provide the server to use to download files from CRDS. To do this,
define two environment variables:

For example, in Bash: 
```
export CRDS_PATH="$HOME/crds_cache"
export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
```
This assumes you will store reference files in a directory `crds_cache` in 
your home directory.

Next, to reproduce our reduction, please specify the same CRDS pipeline
mapping (pmap) that we used in creating the DR0.5 reduction:
```
export CRDS_CONTEXT=jwst_0989.pmap
```

<a name='requiredpackages'></a>
### Additional packages and software
Installing `jwst` will install any necessary dependencies, including 
other STScI packages. If you have previously installed a higher version of 
`jwst`, you may need to downgrade the `stcal` and `stdatamodels` packages to 
work with `jwst` v1.7.2. If you encounter errors in Stage 1 processing related
to missing modules or inconsistent datamodel structures, you may try:
```
pip install stcal==1.1.0
pip install stdatamodels==0.4.3
```
These are the versons of the `stcal` and `stdatamodels` packages that were 
available when v1.7.2 was the current release of `jwst`.

Our source masking procedures (used for the [1/f noise subtracion](stage1.md) 
and [sky background subtraction](stage3.md) require a `Photutils` function 
called `make_source_masks`, which was moved to a different part of the module 
in version 1.6.0. If you have previously installed a higher version of 
`Photutils`, you may need to downgrade. If you encounter errors when running
`remstriping.py`, you may try:
```
pip install photutils==1.5.0
```

Finally, our wrapper for [running TweakReg](tweakreg/README.md#newtweakreg) 
will require Source Extractor, which can be installed via conda:
```
conda install -c conda-forge astromatic-source-extractor
```
Source Extractor is not needed if you choose to skip running TweakReg from 
scratch and instead update the image headers with our saved WCS models 
[as described here](stage3.md#updatewcs).


<a name='download'></a>
### Downloading Data

You can download data from the MAST web portal or using Astropy astroquery.
We provide a script to download all CEERS NIRCam imaging raw files, sorted
by pointing ID. First, download astroquery if it isn't already a part of your
environment:
```
conda install -c conda-forge astroquery
```

To download all raw files from CEERS NIRCam pointing 1:
```
python ceersdownload.py 1
```

To download all 1936 raw CEERS NIRCam files:
```
python ceersdownload.py all
```

All products are downloaded into individual subdirectories within a parent
directory `DOWNLOAD_DIR/mastDownload/JWST`. To move all downloaded files from
their subdirectories (and delete the subdirectories) to `DOWNLOAD_DIR`:
```
python ceersdownload.py 1 --cleanup
```

The download script can optionally also sort all raw files by CEERS NIRCam 
pointing. To move all downloaded files into directories called `nircam1`, 
`nircam2`, etc.: 
```
python ceersdownload.py 1 --sort_files
```
If the `--sort_files` option is set, the subdirectory structure is cleaned
up (`--cleanup`) by default.

**Customization Options:**

At the top of `ceersdownload.py`:
* `DOWNLOAD_DIR`: provide the relative (or absolute) path to the directory
  into which to download the MAST data products. 
* `DATATYPE`: provide the type of data to download (`UNCAL`, `RATE`, `CAL`,etc.) 


<a name='stage1'></a>
## Stage 1: Detector-level corrections 

Stage 1 of the JWST Calibration Pipeline performs detector-level corrections,
many of which are common to all instruments and observing modes.

We have created our own procedures for snowball removal, wisp subtraction 
and 1/f noise subtraction. See [stage1.md](stage1.md) for more information.


<a name='stage2'></a>
## Stage 2: Individual Image Calibrations

Stage 2 of the JWST Calibration Pipeline involves individual image calibrations
such as flat-fielding and flux calibrating the data.

We adopt the default options for this stage. See [stage2.md](stage2.md) for 
more information.


<a name='stage3'></a>
## Stage 3: Ensemble Processing

Stage 3 of the JWST Calibration Pipeline performs ensemble reduction steps,
where the output product is a single mosaic per filter combining all images
and dithers.

We have created our own procedures for astrometric alignment, sky subtraction,
and rescaling the variance maps. See [stage3.md](stage3.md) for more 
information.


<a name='batch'></a>
## Batch Processing

In the directory `batch_scripts`, we provide scripts used to run each of the 
reduction steps on all imaging from the Epoch 1 observations. We used the 
Texas Advanced Computing Center (TACC) at the University of Texas at Austin 
for our image reduction, and so these scripts were designed to be executed in 
parallel using a distributed computing system. These scripts provide one 
command per line. 

The scripts assume an input directory called `uncals` (for stage 1 processing 
only) and all processed outputs are saved in a directory called `calibrated`. 
The READMEs for each Pipeline stage describe how the I/O paths and filenames
can be changed in each reduction script.


