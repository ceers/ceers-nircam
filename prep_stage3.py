import os
from glob import glob
import numpy as np
from astropy.table import Table, Column
from astropy.io import fits

#################################
### SET I/O AND OTHER SETUP HERE 
INPUTDIR = 'calibrated'
OUTPUTDIR = 'calibrated'
#################################


def make_pool(output, filetype='cal'):
    """Create table for use as asn pool file for stage 3

    Args:
        output (str): output filename for table
        filetype (Optional [str]): file name/type to use when searching 
            for the files to include in table. Default is 'cal'
    """
    cals = glob(os.path.join(INPUTDIR, '*%s.fits'%filetype))
    cals.sort()

    base = [x.split('_%s'%filetype)[0] for x in cals]
    t = Table([cals, ['science']*len(cals), base], names=('expname', 
               'exptype','name'))
    t.add_column(Column(data=np.zeros(len(cals), dtype='S10'), name='filter'))
    for i,row in enumerate(t):
        t['filter'][i] = fits.getheader(row['expname'])['FILTER']
    t.write(output, format='ascii', overwrite=True)


def make_asn(asnpool, output, product=None):
    """Make association file using `asn_from_list` and DMS_Level3_Base rule

    Args:
        asnpool (str): name of asn pool file that lists images, base
            file names, and filters
        output (str): base name for output association files, the filter name
            will be added output as in [output]_[filter].json
        product (Optional [str]): base name for Stage 3 output image, the 
            filter name will be added to product as [product]_[filter]. 
            This base name is then used for Stage 3 output, as in 
            [product]_[filter]_i2d.fits. If None, product will be the same
            as output. Default is None.
    """
    if product is None:
        product = output

    t = Table.read(asnpool, format='ascii')
    filts = np.unique(t['filter'])
    print('Found filters: \n', filts)
    for filt in filts:
        w = np.where(t['filter'] == filt)
        calfiles = list(t['expname'][w])
        files = ' '.join(calfiles)
        jsonfile = '%s_%s.json'%(output, filt.lower())
        cmd = 'asn_from_list -o %s -r DMS_Level3_Base --product-name %s_%s %s'%(jsonfile,product,filt.lower(),files)
        os.system(cmd)


def main():
    field = 'nircam1'

    # create a pool to use in generating an asn file for stage 3
    make_pool('%s_pool.dat'%field, filetype='cal')

    # create asn files from pool file
    make_asn('%s_pool.dat'%field, field)


if __name__ == '__main__':
    main()


