#!/usr/bin/env python

import pyfits,sys

f='/projects/quiet/external_data/wmap_lcdm_sz_lens_wmap7_cl_v4.fits'

hdulist=pyfits.open(f)
data=hdulist[1].data
table_header=hdulist[1].header
hdulist.close()
#print header

for row_number,row in enumerate(data):
    print row_number,row[0]

sys.exit(0)
