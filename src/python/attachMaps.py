#!/usr/bin/env python

"""
Read a TOD gif and a map png
Determine their dimensions
Attach them and output a joined png
"""

#-------------------------------------------------------------------

def readImage(f):

    """
    """

    import Image

    image = Image.open(f)

    return image

#-------------------------------------------------------------------

def joinImages(i1,i2,style):

    """
    """

    import Image,sys

    (s1a,s1b)=i1.size
    (s2a,s2b)=i2.size

    crop=False
    if crop:
        i2=i2.crop((0,0,s2a/2,s2b)).copy()
        i2.size=[2*el for el in i2.size]
        (s2a,s2b)=i2.size
        print (s1a,s1b)
        i1.resize((int(0.5*s1a),int(0.5*s1b)),Image.ANTIALIAS)
        (s1a,s1b)=i1.size
        print i1.size,i2.size


    # this needs to be much more sophisticated
    if style == 'oandu':
        h=max(s1a,s2a)
        v=s1b+s2b
    elif style == 'sbys':
        h=s1a+s2a
        v=max(s1b,s2b)
    else:
        print 'style not set'
        sys.exit(0)

    j=(h,v)
    #print j

    format='RGB'
    merged_image = Image.new(format,j)

    # as does this..
    o1=(0,0)
    if style == 'oandu':
        o2=(0,s1b)
    elif style == 'sbys':
        o2=(s1a,0)

    #print merged_image.size

    merged_image.paste(i1,o1)
    merged_image.paste(i2,o2)

    print '%s + %s -> %s' % (i1.size,i2.size,merged_image.size)

    return merged_image

#-------------------------------------------------------------------

def writeImage(image,f):

    """
    """

    image.save(f)

    print 'Joint image written to %s' % f

    return

#-------------------------------------------------------------------

def parseTODName(f):

    """
    """

    # To avoid full path issues
    t=f.split('/')[-1]

    s=t.split('_')
    tces  = int(s[0][2:])
    tmod  = int(s[1][2:])
    tdi   = int(s[3])
    tstat = int(s[4].split('.')[0][2:])

    return (tces,tmod,tdi,tstat)

#-------------------------------------------------------------------

def unparseTODName(tces,tmod,tdi,tstat,joint=True):

    """
    """

    if not joint:
        print 'Assumed wrongly that you want to protect input files'
        print 'Exiting...'
        sys.exit(0)

    fout = 'ML%i_RW%02i_joint_%i_st%i.png' % (tces,tmod,tdi,tstat)

    return fout

#-------------------------------------------------------------------

def parseMapName(f):

    """
    """

#    f='patch_6a_set002_ces3436_P_bin.fits'
#    f='patch_7b_set001_ces1670_ass080_P_bin.png'

    j=f.split('/')[-1]
    s=j.split('_')

    mass=None
    # I think this test could be shortened (i.e. for loop removed) - works OK
    for el in s:
        if 'ass' not in el:
            mpatch = s[1]
            mces   = int(s[3][3:])
            mpol   = s[-2]
        else:
            mpatch = s[1]
            mces   = int(s[3][3:])
            mass   = int(s[4][3:])
            mpol   = s[-2]

    return (mpatch,mces,mass,mpol)

#-------------------------------------------------------------------

def attachMap(f1,f2,style):

    """
    Attach a single TOD image to a single map image
    """

    (tces,tmod,tdi,tstat)=parseTODName(f1)
    (mpatch,mces,mass,mpol)=parseMapName(f2)

    print 'TOD ML%i RW%02i di%i st%i' % (tces,tmod,tdi,tstat)
    print 'MAP ML%i RW%02i patch_%s poln=%s' % (mces,mass-1,mpatch,mpol)
    #print tces,tmod,tdi,tstat
    #print mces,mpatch,mpol


    mces=tces ### TEMPORARY FOR TESTING

    assert(mces==tces), '***CES mismatch!! %i /= %i' % (mces,tces)
    assert(mass==tmod+1), '***CES mismatch!! %i /= %i' % (mass,tmod)

    i1=readImage(f1)
    i2=readImage(f2)

    i=joinImages(i1,i2,style)

    fout = unparseTODName(tces,tmod,tdi,tstat,joint=True)
    writeImage(i,fout)

    return fout

#-------------------------------------------------------------------

if __name__ == '__main__':

    from optparse import OptionParser

    print
    # Set up options
    usage = "./attachMaps.py -h\n"
    parser = OptionParser(usage=usage)

    #parser.add_option("-v", "--verbose", action="store_false",\
    #                  dest="verbose",help="make lots of noise")

    parser.add_option("-b", "--batch", type="string",\
                      dest="batch",help="batch mode (give dir)")

    parser.add_option("-s", "--style", type="string",\
                      dest="style",help="Set join style",
                      metavar="(oandu|sbys)")

    parser.add_option("-i", "--infi1", dest="infi1",\
                      help="specify TOD", type="string",\
                      metavar='')

    parser.add_option("-j", "--infi2", dest="infi2",\
                      help="specify map", type="string",\
                      metavar='')

    parser.add_option("-d", "--indi1", dest="indi1",\
                      help="specify TODs dir", type="string",\
                      metavar='')

    parser.add_option("-e", "--indi2", dest="indi2",\
                      help="specify maps dir", type="string",\
                      metavar='')


    # Assert defaults
    parser.set_defaults(infi1=None,infi2=None,dir1=None,\
                        dir2=None,style=None,batch=None)

    # Parse command line
    (options, args) = parser.parse_args()

    f1=options.infi1
    f2=options.infi2
    dir1=options.indi1
    dir2=options.indi2
    style=options.style
    batch=options.batch

    import os,sys,glob

    #dir1='/projects/quiet/jonathaz/wband/cuts/tods/test-dir2'
    #dir2='/projects/quiet/jonathaz/wband/early-maps-1103/7b/patch_7b/test'
    #f1='patch_7b_set001_ces3571_P_bin.gif'
    #f1='ML861_RW03_all_0_st0.png'
    #f2='patch_7b_set003_ces6223_P_bin.gif'

    if batch is None:
        attachMap(os.path.join(dir1,f1),os.path.join(dir2,f2),style)
    else:
        dir1=os.path.join(batch,'tods')
        dir2=os.path.join(batch,'ass-pngs')

        tods=glob.glob(os.path.join(dir1,'*_?_*.png'))
        maps=glob.glob(os.path.join(dir2,'*'))

        for map in maps:
            (mpatch,mces,mass,mpol)=parseMapName(map)
            for tod in tods:
                (tces,tmod,tdi,tstat)=parseTODName(tod)
                if mces==tces and mass==tmod+1:
                    print mces,tces,mass,tmod+1
                    attachMap(tod,map,style)

    print '\nFinished!\n'
    sys.exit(0)

