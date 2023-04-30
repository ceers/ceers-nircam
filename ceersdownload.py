import os
import shutil
from glob import glob
import argparse
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack, unique
from astroquery.mast import Observations


#################################
### SET I/O AND OTHER SETUP HERE 
# Destination directory for downloaded images
DOWNLOAD_DIR = 'uncals'
# Datatype to download, options are UNCAL, RATE, CAL, I2D
DATATYPE = 'UNCAL'
#################################


def download_ceers(pointing, cleanup=False, sort_files=False):
    """ """
    # if files will be sorted into pointing directories after downloading, 
    # then the cleanup option must also be True
    if sort_files:
        cleanup = True

    # Download CEERS
    propid = '1345'
    instrument = 'NIRCAM/IMAGE'
    filts = ['F115W','F150W','F200W','F277W','F356W','F410M','F444W']

    # get list of CEERS observations
    ceers = Table.read('ceersobs.dat', format='ascii')

    if pointing == 'all':
        obs_table = Observations.query_criteria(obs_collection='JWST',
                        filters=filts, proposal_id=[propid],
                        instrument_name=instrument, 
                        dataproduct_type='IMAGE')
        
    else:
        w = np.where(ceers['pointing'] == int(pointing))
        obs = np.unique(ceers['observation'][w])
        for o in obs:
            obsid = '*o%03d*'%o
            obs_t = Observations.query_criteria(obs_collection='JWST',
                        filters=filts, proposal_id=[propid],
                        instrument_name=instrument, obs_id=obsid,
                        dataproduct_type='IMAGE')

            try:
                obs_table = vstack([obs_table, obs_t])
            except NameError:
                obs_table = obs_t.copy()

    products_list = Observations.get_product_list(obs_table)
    products = Observations.filter_products(products_list,
                                          productType=['SCIENCE', 'INFO'],
                                          productSubGroupDescription=DATATYPE)

    files = unique(products, keys='productFilename')
    files.pprint()

    try:
        os.makedirs(DOWNLOAD_DIR)
    except FileExistsError:
        # directory already exists
        pass
        

    # the development version of astroquery has an option to download files
    # with a flattened structure (flat=True), but this doesn't exist in 
    # current, 0.4.6 version of astroquery, so do this manually after download
    manifest = Observations.download_products(files, #flat=True, 
                                               download_dir=DOWNLOAD_DIR)

    if cleanup:
        # get list of all downloaded files
        downloads = glob(os.path.join(DOWNLOAD_DIR, 
                        'mastDownload/JWST/jw01345*/*%s.fits'%DATATYPE.lower()))
        for downloaded_file in downloads:
            if sort_files:
                # determine pointing by checking obs number
                obs = int(fits.getheader(downloaded_file)['OBSERVTN'])
                pointing = ceers['pointing'][ceers['observation'] == obs][0]
                pointing_dir = os.path.join(DOWNLOAD_DIR, 'nircam%i'%pointing)
                try:
                    os.makedirs(pointing_dir)
                except FileExistsError:
                    pass
                shutil.move(downloaded_file, pointing_dir)

            else:
                shutil.move(downloaded_file, DOWNLOAD_DIR)

            # remove the corresponding MAST directory if it's empty
            download_dir = os.path.dirname(downloaded_file)
            if len(os.listdir(download_dir)) == 0:
                os.rmdir(download_dir)
        
   

def main():
    parser = argparse.ArgumentParser(description='Download CEERS NIRCam images from MAST')
    parser.add_argument('pointing', type=str,
                        help='CEERS NIRCam pointing ID to download: 1, 2, ... or all to get all 1936 CEERS images.')
    parser.add_argument('--cleanup', action='store_true',
                        help='Clean up downloaded files by moving them from mastDownload/JWST/* subdirectories to a single directory specified by DOWNLOAD_DIR.')
    parser.add_argument('--sort_files', action='store_true',
                        help='Sort downloaded files into directories by CEERS NIRCam pointing (nircam1, nircam2, ...). Sets --cleanup by default')

    args = parser.parse_args()

    download_ceers(args.pointing, cleanup=args.cleanup, 
                   sort_files=args.sort_files)


if __name__ == '__main__':
    main()

