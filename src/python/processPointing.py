#!/usr/bin/env python

"""
Read an ASCII file containing pointing fits
Tabulate them
Output to file
"""

#-------------------------------------------------------------------

def djd2jd(djd):
    return djd + 2415020.0
def jd2mjd(jd):
    return jd - 2400000.5
def djd2mjd(djd):
    return jd2mjd(djd2jd(djd))

#-------------------------------------------------------------------

def utc2mjd(utc):

    """
    Convert '%d-%b-%Y %H:%M:%S' string to MJD float
       e.g. 09-FEB-2010 01:37:08 -> 55236.0674537
    """

    import time,ephem

    time_struct=time.strptime(utc,'%d-%b-%Y %H:%M:%S')
    ephem_object=ephem.Date(tuple(time_struct)[:6])

    mjd = djd2mjd(float(ephem_object))

    return mjd

#-------------------------------------------------------------------

def convertFromOmegaTheta((omega,theta)):

    """
    Does what it says on the tin
    Arguments and outputs are in **degrees**
    """

    from math import sin,cos,asin,degrees,radians

    (omega,theta) = [radians(el) for el in (omega,theta)]

    sin_theta_lat =  cos(omega) * sin(theta)
    sin_theta_HA  = -sin(omega) * sin(theta)

    theta_lat = asin(sin_theta_lat)
    theta_HA  = asin(sin_theta_HA)

    (theta_lat,theta_HA) = [degrees(el) for el in (theta_lat,theta_HA)]

    return (theta_lat,theta_HA)

#-------------------------------------------------------------------

def convertToOmegaTheta((theta_lat,theta_HA)):

    """
    Does what it says on the tin
    Arguments and outputs are in **degrees**
    """

    from math import sin,asin,atan2,sqrt,radians,degrees

    (theta_lat,theta_HA) = [radians(el) for el in (theta_lat,theta_HA)]

    omega = atan2(-sin(theta_HA),sin(theta_lat))

    arg_sin = sin(theta_lat)**2 + sin(theta_HA)**2
    sin_theta = sqrt(arg_sin)
    theta = asin(sin_theta)

    (omega,theta) = [degrees(el) for el in (omega,theta)]

    return (omega,theta)

#-------------------------------------------------------------------

if __name__ == '__main__':

    import sys,os,string
    from optparse import OptionParser
    import scanf

    print
    # Set up options
    usage = "./processPointing.py -h.py\n"
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose", action="store_false",\
                      dest="verbose",help="make lots of noise")

    parser.add_option("-i", "--infile", dest="infile",\
           help="specify input file", type="string",\
           metavar='fits.txt')

    parser.add_option("-o", "--outfile", dest="outfile",\
           help="specify outfile", type="string",\
           metavar='optical.dat')

    parser.add_option("-p", "--postfile", dest="postfile",\
           help="specify postfile", type="string",\
           metavar='post.txt')

    # Assert defaults
    parser.set_defaults(infile='fits.txt',object='optical.dat',\
                        postfile='post.txt',verbose=False)

    # Parse command line
    (options, args) = parser.parse_args()

    verbose=options.verbose
    infile=options.infile
    outfile=options.outfile
    postfile=options.postfile

    assert(infile is not None), 'ces_stats required!'

    indata=open(infile,'r')

    print 'Reading %s....' % infile

    #lflex=' [0]   Flexure parameter       = %f deg = %f arcsec'

    form={\
        0:' [0]   Flexure parameter      = %f deg = %f arcsec',\
        1:' [1]   AZ axis tilt in HA     = %f deg = %f arcsec',\
        2:' [2]   AZ axis tilt in Lat    = %f deg = %f arcsec',\
        3:' [3]   EL axis tilt           = %f deg = %f arcsec',\
        4:' [4]   Collimation error      = %f deg = %f arcsec',\
        5:' [5]   PA offset of coll. err =     %f deg',\
        6:' [6]   AZ encoder offset      = %f deg = %f arcsec',\
        7:' [7]   EL encoder offset      = %f deg = %f arcsec',\
        8:' [A]   x tilt zero            = %f deg = %f arcsec',\
        9:' [B]   y tilt zero            = %f deg = %f arcsec',\
        10:'# UT range %s to %s',\
        11:'# Residuals: rms AZ=%f, AZcos(EL)=%f, EL=%f, sky=%f; max=%f arcsec'}

    key={0:' [0]',\
         1:' [1]',\
         2:' [2]',\
         3:' [3]',\
         4:' [4]',\
         5:' [5]',\
         6:' [6]',\
         7:' [7]',\
         8:' [A]',\
         9:' [B]',\
         10:'# UT range',\
         11:'# Residuals:'}

    if verbose:
        print '# Run runid flex_deg flex_asec aztilt_HA_deg aztilt_HA_asec',
        print 'aztilt_lat_deg aztilt_lat_asec eltilt_deg eltilt_asec',
        print 'coll_deg coll_asec PA_deg azenc_deg azenc_asec elenc_deg',
        print 'elenc_asec x_deg x_asec y_deg y_asec start_date start_time',
        print 'stop_date stop_time MJD_start MJD_stop niter'
#        print 'resid_az resid_coaz resid_el resid_sky resid_max'

    if os.path.exists(outfile): os.rename(outfile,'%s.old'%outfile)
    outdata=open(outfile,'w')
    outdata.write('# runid flex_deg flex_asec aztilt_HA_deg aztilt_HA_asec ')
    outdata.write('aztilt_lat_deg aztilt_lat_asec eltilt_deg eltilt_asec ')
    outdata.write('coll_deg coll_asec PA_deg azenc_deg azenc_asec elenc_deg ')
    outdata.write('elenc_asec x_deg x_asec y_deg y_asec start_date start_time ')
    outdata.write('stop_date stop_time MJD_start MJD_stop niter\n')
#    outdata.write('resid_az resid_coaz resid_el resid_sky resid_max\n')

    #print utc MJD

    #data={}
    master=[]
    print '..Massaging formats..'
    for line in indata:
        try:
            if line.startswith('Did'):
                data[13]=[int(line.split()[1])]
                # FAIL data[13] = [scanf.sscanf(line,'Did %i iterations of gridls')]
            if line.startswith('Model'):
                data={}

            for param in sorted(key.keys())[0:10]:
                if line.startswith(key[param]):
                    data[param]=(scanf.sscanf(line,form[param]))
                    #print data[param]

            # Overwrite the time collected
            if line.startswith('# UT range'):
                data[10]=\
                  (' '.join(line.split()[3:5]),' '.join(line.split()[6:8]))
                data[12]=[utc2mjd(el) for el in data[10]]
                #print data[10],data[12]

            # Wrap up
            if line.startswith('collimate'):
                master.append(data)
                #print 'xx',data[10],data[12]

        except:
            #print 'Skipping %s' % line
            pass

    indata.close()

    #print master[0][11]
    #sys.exit(0)
    #for run in range(len(master)):
    #    print run,len(master[run].keys()),master[run][10]
    #print len(master[5].keys())

    for run in range(len(master)):
        if verbose: print 'Run %i' % run,
        outdata.write('%i ' % run)
        for key in sorted(master[run].keys()):
            s=' '.join([str(e) for e in master[run][key]])
            if verbose: print s,
            outdata.write('%s '% s)
            #print '%i %s' %(key, master[0][key]),
        if verbose: print ''
        outdata.write('\n')

    outdata.close()
    print '....Wrote %s' % outfile

    postprocess=True

    if postprocess:
        print '..Massaging formats again..'
        if os.path.exists(postfile): os.rename(postfile,'%s.old'%postfile)
        finaldata=open(postfile,'w')

        finaldata.write('# MJD     -- Starting date/time for parameters from this row to take effect\n')
        finaldata.write('# DECK    -- Deck offset, in degrees\n')
        finaldata.write('# FLEX    -- CBI mount flexure parameter\n')
        finaldata.write('# OMEGA   -- Orientation of the az axis tilt, in degrees E of N\n')
        finaldata.write('# THETA   -- Magnitude of the az axis tilt\n')
        finaldata.write('# THETA_E -- El axis non-perpendicularity\n')
        finaldata.write('# THETA_C -- Radio collimation angle\n')
        finaldata.write('# PSI_C   -- Radio collimation orientation\n')
        finaldata.write('# THETA_O -- Optical collimation angle\n')
        finaldata.write('# PSI_O   -- Optical collimation orientation\n')

        units='deg'
        #units='arcsec'
        if units == 'deg':
            unit=0
        elif units == 'arcsec':
            unit=1

        #verbose=True

        # Sort data by (starting) MJD and remove online run
        sorted_master = sorted(master,key=lambda master:master[12][0])

        for entry in sorted_master[1:]:
            mjd=(entry[12][0]+entry[12][1])/2.0
            mj='%13.7f ' % mjd
            if verbose: print mj,
            finaldata.write(mj)

            dk_offset = -99.0
            d_o='%3.1f ' % dk_offset
            if verbose: print d_o,
            finaldata.write(d_o)

            flex=entry[0][unit]
            fl='%f ' % flex
            if verbose: print fl,
            finaldata.write(fl)

            theta_HA = entry[1][unit]
            theta_lat = entry[2][unit]
            (omega,theta)=convertToOmegaTheta((theta_lat,theta_HA))
            #print 'x',theta_lat,theta_HA
            #print 'x',convertFromOmegaTheta((omega,theta))
            ot='%f %f ' % (omega,theta)
            if verbose: print ot,
            finaldata.write(ot)

            el_axis_tilt=entry[3][unit]
            theta_e = el_axis_tilt
            t_e='%f ' % theta_e
            if verbose: print t_e,
            finaldata.write(t_e)

            theta_c = -99.0
            psi_c   = -99.0
            theta_o = -99.0
            psi_o   = -99.0

            r='%f %f %f %f' % (theta_c,psi_c,theta_o,psi_o)
            if verbose: print r
            finaldata.write(r)
            finaldata.write('\n')

#        0:' [0]   Flexure parameter      = %f deg = %f arcsec',\
#        1:' [1]   AZ axis tilt in HA     = %f deg = %f arcsec',\
#        2:' [2]   AZ axis tilt in Lat    = %f deg = %f arcsec',\
#        3:' [3]   EL axis tilt           = %f deg = %f arcsec',\
#        4:' [4]   Collimation error      = %f deg = %f arcsec',\
#        5:' [5]   PA offset of coll. err =     %f deg',\
#        6:' [6]   AZ encoder offset      = %f deg = %f arcsec',\
#        7:' [7]   EL encoder offset      = %f deg = %f arcsec',\
#        8:' [A]   x tilt zero            = %f deg = %f arcsec',\
#        9:' [B]   y tilt zero            = %f deg = %f arcsec',\

        finaldata.close()
        print '....Wrote %s' % postfile

    sys.exit(0)

