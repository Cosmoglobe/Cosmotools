#!/usr/bin/python

import sys

def usage():
    sys.stderr.write( \
"""
Usage: ./makeDifferenceMap.py map1.fits map2.fits output_difference.fits (divide=[True|False])
Output = map1 - map2
Will only work on full healpix maps.\n
""")
    sys.exit(1)

# -------------------------------------------------------------

if __name__ == '__main__':

    # *Divide* first map by second map; otherwise
    # *subtract* second from first
    divide=False
    # divide=True

    # Read the file names from the arguments to the program
    # Optionally read the divide keyword
    try:
        mapname1,mapname2,outname = sys.argv[1:4]
    except:
        try:
            mapname1,mapname2,outname,divide = sys.argv[1:5]
	except:
	    usage()

    import pyfits
    from numpy import where,logical_and

    # The healpix "unseen pixel" value. We do not want to difference any
    # pixels with this value in either field
    from healpy import UNSEEN

    # The set of fields we want to difference.  We want to exclude any
    # fields like "pixel number" or "signal to noise"
    fields = ('TEMPERATURE',
	      'Q-POLARISATION','U-POLARISATION',
	      'Q-POLARIZATION','U-POLARIZATION',
	      'Q_POLARISATION','U_POLARISATION',
	      'Q_POLARIZATION','U_POLARIZATION',
	      'Q_MAP','U_MAP')

    # Use pyfits to load the list of HDUs in the map
    map1 = pyfits.open(mapname1)
    map2 = pyfits.open(mapname2)

    # Fix tod2map-type header keywords
    map1[1].header['TTYPE1']='PIX_NUM'
    map1[1].header['TTYPE2']='TEMPERATURE'
    map1[1].header['TTYPE3']='Q-POLARIZATION'
    map1[1].header['TTYPE4']='U-POLARIZATION'
    map2[1].header['TTYPE1']='PIX_NUM'
    map2[1].header['TTYPE2']='TEMPERATURE'
    map2[1].header['TTYPE3']='Q-POLARIZATION'
    map2[1].header['TTYPE4']='U-POLARIZATION'

    # Get the data in the first extension of each map
    # This is where healpix stores its data
    data1=map1[1].data
    data2=map2[1].data

    # Loop through the fields in the map differencing each one
    for field in data1.names:
        # Ignore miscellaneous fields
	if field not in fields:
	    print "Skipping field %s" % field
	    continue
        print "Differencing field %s" % field
	
	# Extract the data in this column
	col1=data1.field(field)
	col2=data2.field(field)

	# Find the pixels that are seen in both maps
	good_pixels = where(logical_and(col1!=UNSEEN,col2!=UNSEEN))
	print "%d pixels remain in field" % len(good_pixels[0])
	
	# Get the difference of the full columns
	if divide:
	    # This will probably cause problems if there are
	    # any zeros in the second map
	    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.divide.html
	    col_diff = col1/col2
	else:
	    col_diff = col1-col2
	
	# First re-set column1 to unseen
	col1[:]=UNSEEN
	# then set the twice-seen pixels to the right difference.
	# this modifies the first map but does not save the difference
	# back to original file
	col1[good_pixels] = col_diff[good_pixels]

    # save the first map back to the new file name
    map1.writeto(outname)

    sys.exit(0)

