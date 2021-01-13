#!/usr/bin/python

"""
For a healpix map, calculate statistics for each polarization
(min, max, mean, rms)
"""

#-------------------------------------------------------------------

if __name__ == '__main__':

    import os,sys,pyfits
    from math import sqrt
    import numpy
    from healpy import UNSEEN,nside2npix

    f=sys.argv[-1]

    fits=pyfits.open(f)

    fields = ['TEMPERATURE',
              'Q-POLARISATION','U-POLARISATION',
              'Q-POLARIZATION','U-POLARIZATION',
              'Q_POLARISATION','U_POLARISATION',
              'Q_POLARIZATION','U_POLARIZATION',
              'Q_MAP','U_MAP']

    maps=fits['FULL SKY MAP']
    nside=maps.header['NSIDE']

    print 'For %s (nside %i),' % (os.path.basename(f), nside)
    for field in maps.data.names:
        # Ignore miscellaneous fields
        if field not in fields:
            print 'Skipping field %s' % field
            continue
        if field == 'TEMPERATURE': continue # Skip temperature data
        #print 'Working on field %s' % field
        map=maps.data.field(field)
        
        vector=map.flatten()
        assert(len(vector)==nside2npix(nside)), '** npix discrepant!!'

        mask=numpy.where(vector != UNSEEN)

        j=vector[mask]

        n=len(j)
        assert(n != 0), '** zero elements in masked array!'
        mean=j.mean()
        std=j.std()
        rms=std/sqrt(float(n))
        min=j.min()
        max=j.max()

        print '%s: %i elements, %f < (%f +/- %f) < %f uK' \
              % (field[0],n,min,mean,rms,max)

    sys.exit(0)
