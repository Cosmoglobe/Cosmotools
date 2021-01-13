#!/usr/bin/python

"""
Utilities to convert pixels.dat list into FITS using a template FITS
and (optionally) to weight that template FITS by the pixels.dat
and to output the weighted map as FITS
and (optionally) to calculate and output a Chi^2(map|null) map
"""

# ----------------------------------------------------------

def plot_chisq(array, label, fig, chiofmap):
    from matplotlib import pyplot as plt
    import pylab
    import string
    params={'axes.labelsize':8,'ytick.labelsize':6,'xtick.labelsize':6,'legend.font size':6}
    pylab.rcParams.update(params)

    if label!='TEMPERATURE':
       lower = 0
       upper = 0.05*max(array)
       chilabel=str(chiofmap)
       plt.hist(array,bins=1000, histtype='bar',range=(lower,upper))
       plt.ylabel('counts')
       plt.xlabel(label)
       #plt.title(string.join(['Chi: ',chilabel],''))
       plt.savefig(fig)
       plt.clf()
       plt.close('all')

    return

#-----------------------------------------------------------

def get_weights(noisemap):

    import numpy
    from numpy import where,logical_and
    import pyfits, healpy
    import quiet, differenceNoiseMap
    #from quiet import read_pixfile

    #The healpix "unseen pixel" value.  We do not want to difference
    #any pixels with this value in either field
    from healpy import UNSEEN

    # get noise pixels
    #nside,pix,vals=quiet.read_pixfile(noisemap)
    nside,pix,vals=differenceNoiseMap.read_map(noisemap)
    full_map=numpy.repeat(UNSEEN,healpy.nside2npix(nside))
    full_map[pix]=vals

    return full_map,nside

#-----------------------------------------------------------

def pixels2map(pixels,vals,nside):

    """
    Convert a list of pixels to a 1-D map, initialized with UNSEENs
    """

    from numpy import repeat
    from healpy import UNSEEN, nside2npix

    npix=nside2npix(nside)
    mapp=repeat(UNSEEN,npix)
    mapp[pixels]=vals

    return mapp

#-----------------------------------------------------------

#def 1dmapto2dmap(1dmap,shape):
#    """
#    """
#    2dmap=1dmap
#    return 2dmap

#-----------------------------------------------------------

def goodpix(map1,map2):

    """
    Check two Healpix arrays for non-UNSEEN pixels and return mask of
    good pixels, i.e. those seen in both maps
    """

    from numpy import where,logical_and
    from healpy import UNSEEN

    return where(logical_and(map1!=UNSEEN,map2!=UNSEEN))

#-----------------------------------------------------------

def pixels2fits(pixelsf=None,template=None,chi2=False,write=True,divide=True):

    """
    Convert a healpix pixels file to FITS format, using a FITS map as
    a template
    """

    import os, pyfits, numpy, healpy, string
    import differenceNoiseMap
    from healpy import UNSEEN
    from secateur import StripExtension, FixPerms

    rootdir='/usit/titan/u1/jonathaz/quiet_data/nulltests/pcl/pcl-2009123-001/maps/'
    if pixelsf is None: pixelsf=os.path.join(rootdir,'pcl-2009123-001-patch_6a_tandem.diff_pixels.txt')
    if template is None: template=os.path.join(rootdir,'pcl-2009123-001-patch_6a_tandem.diff_map.fits')

    # Handle the input files
    if not os.path.isabs(pixelsf): pixelsf=os.path.abspath(pixelsf)
    if not os.path.isabs(template): template=os.path.abspath(template)

    # Read in the pixels file
    nside,pixels,vals=differenceNoiseMap.read_map(pixelsf)
    #nside,pix,vals=quiet.read_pixfile(noisemap)
    #rms,nside=get_weights(pixels)

    # Turn the pixels list into a numpy 1-D map
    rmsmap=pixels2map(pixels,vals,nside)

    # Read in the template FITS
    template_file=pyfits.open(template)

    # Check the pixels map and template FITS map have the same nside
    location='FULL SKY MAP'
    templatemap=template_file[location]
    nside_template=templatemap.header['NSIDE']
    assert(nside_template==nside,'Maps have different nsides')

    # The list of fields we want to difference. We want to exclude any
    # fields like "pixel number" or "signal to noise"
    fields = ['TEMPERATURE',
              'Q-POLARISATION','U-POLARISATION',
              'Q-POLARIZATION','U-POLARIZATION',
              'Q_POLARISATION','U_POLARISATION',
              'Q_POLARIZATION','U_POLARIZATION',
              'Q_MAP','U_MAP']

    ff=[]
    chi2f=[]
    rmsdone=False
    for field in templatemap.data.names:
        # Ignore miscellaneous fields
        if field not in fields:
            print "Skipping field %s" % field
            continue
        print "Working on field %s" % field
        fieldmap=templatemap.data.field(field)
        # Reshape the 1-D rms map into a 2-D map (once only)
        if not rmsdone:
            rms2d=numpy.reshape(rmsmap,numpy.shape(fieldmap))
            rmsdone=True
        fieldgood=goodpix(rms2d,fieldmap)
        wfield=numpy.tile(UNSEEN,fieldmap.shape)
        wfield[fieldgood]=fieldmap[fieldgood]/(rms2d[fieldgood]**2)
        fieldflat=wfield.flatten()
        ff.append(fieldflat)
        if chi2:
            chi2field=numpy.tile(UNSEEN,fieldmap.shape)
            chi2field[fieldgood]=fieldmap[fieldgood]**2/(rms2d[fieldgood]**2)
            chi2flat=chi2field.flatten()
            chi2f.append(chi2flat)
            #dohist=False
            dohist=True
            if dohist:
                # (n,bins) = numpy.histogram(chi2flat,bins=nbins,normed=False)
                print 'Max', max(chi2flat)
                seen=numpy.where(chi2flat!=UNSEEN)
                c2=chi2flat[seen].sum()
                n=len(chi2flat[seen])
                if n != 0:
                    c2red=c2/float(n)
                else:
                    c2red = -99.0
                print 'Map Chi^2 %6.1f / %i = %4.2f' % (c2,n,c2red) # Laura look here
                # I think Laura's 'chiofmap' should = Jon's c2red
                chiofmap=sum(chi2flat)/len(chi2flat)
                plot_chisq(chi2flat, field, string.join([pixelsf,'_',field,'_','hist.png'],''), chiofmap)

    if write:
        # The x temp file is pretty horrible!
        # Needed because of healpy's dodgy output headers
        outfile=StripExtension(os.path.basename(template))[0]+'.weighted.fits'
        location = 'XTENSION'
        if os.path.exists(outfile):
            os.rename(outfile,outfile+'.old')
        healpy.write_map(outfile,ff)
        # Fix the healpy header..
        rein=pyfits.open(outfile)
        rein[location].header=templatemap.header
        if os.path.exists('x'+outfile): os.remove('x'+outfile)
        rein.writeto('x'+outfile)
        os.remove(outfile)
        os.rename('x'+outfile,outfile)
        FixPerms(outfile)
        # ? NEED ? rein['PRIMARY'].header.update('BITPIX',32,'array data type')
        #rein[location].header.update('TTYPE1','TEMPERATURE','')
        #rein[location].header.update('TTYPE2','Q-POLARISATION','')
        #rein[location].header.update('TTYPE3','U-POLARISATION','')
        #rein[location].header.update('EXTNAME','FULL SKY MAP','')

        rmsfile=StripExtension(os.path.basename(pixelsf))[0]+'.rms.fits'
        if os.path.exists(rmsfile):
            os.rename(rmsfile,rmsfile+'.old')
        healpy.write_map(rmsfile,[rmsmap,rmsmap,rmsmap])
        # Fix the healpy header..
        rein=pyfits.open(rmsfile)
        rein[location].header=templatemap.header
        if os.path.exists('x'+rmsfile): os.remove('x'+rmsfile)
        rein.writeto('x'+rmsfile)
        os.remove(rmsfile)
        os.rename('x'+rmsfile,rmsfile)
        FixPerms(rmsfile)

        if chi2:
            chi2file=StripExtension(os.path.basename(pixelsf))[0]+'.chi2.fits'
            if os.path.exists(chi2file):
                os.rename(chi2file,chi2file+'.old')
            healpy.write_map(chi2file,chi2f)
            # Fix the healpy header..
            rein=pyfits.open(chi2file)
            rein[location].header=templatemap.header
            if os.path.exists('x'+chi2file): os.remove('x'+chi2file)
            rein.writeto('x'+chi2file)
            os.remove(chi2file)
            os.rename('x'+chi2file,chi2file)
            FixPerms(chi2file)
            return ff, rmsmap, rmsfile, outfile, chi2file

    return ff, rmsmap, rmsfile, outfile

#-----------------------------------------------------------

if __name__ == '__main__':

    import sys

    narg=2
    if len(sys.argv)!=narg+1:
        print >> sys.stderr, '%s <pixels.dat> <template.fits>' % sys.argv[0]
        #print ''
        sys.exit(2)

    import os,re

    pixels,template=sys.argv[1:narg+1]

    # Give absolute paths to the input files
    if not os.path.isabs(pixels): pixels=os.path.abspath(pixels)
    if not os.path.isabs(template): template=os.path.abspath(template)

    chi2=True
    #chi2=False
    outputs=pixels2fits(pixels,template,chi2,write=True)
    if chi2:
        (weightedmaps,rmsmap,rmsfile,weightf,chi2f)=outputs
    else:
        (weightedmaps,rmsmap,rmsfile,weightf)=outputs

    print '\n%s, %s -> %s, %s\n' % \
          (os.path.basename(pixels),os.path.basename(template),\
           os.path.basename(rmsfile),os.path.basename(weightf))

    if chi2: print 'chi2 map is %s' % os.path.basename(chi2f)

    secondhalf=re.split('patch',rmsfile)[1]
    patch=re.split('_',secondhalf)[1]

    print 'now ~ do patch2img patch_%s %s -cel -P -min -100 -max 100' \
          % (patch,rmsfile)
    print 'and ~ do: patch2img patch_%s %s -cel -P -min -1e-4 -max 1e-4' \
          % (patch,weightf)\

    if chi2:
        print 'and ~ do: patch2img patch_%s %s -cel -P -min -XXX -max XXX' \
              % (patch,chi2f)

    sys.exit(0)

# -----------------------------------------------------------
