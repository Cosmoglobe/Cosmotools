#!/usr/bin/python

##!/usr/bin/python2.4

"""
Collection of utilities for processing cuts on QUIET L2 data \n
$Author: roberdum $
$Date: 2010-10-19 06:45:29 +0200 (Tue, 19 Oct 2010) $
$Id: secateur.py 1739 2010-10-19 04:45:29Z roberdum $
History: Little worth mentioning
Dependencies: numpy, APEX_weather.py, quiet.py
Usage:
      developed on python 2.4-2.5
      ./secateur.py cuttype [cut-params] runlist.txt outlist.cut \n
      ./secateur.py -f cuts.cmd \n
      ./secateur -h  <--- Type this for help! \n
Footnotes: \n
0. There may be a pause of a few seconds the first time you run, while
     Python compiles the code .py -> .pyc; fast subsequently
1. Tests for equality to zero should be replaced by machine epsilon test
"""

import sys,platform,os,time,re
from numpy import *
from APEX_weather import *
from optparse import OptionParser
#import numarray
#import numpy

# Set up the default L2 runlist, plus backups
default_rundir  = '/data4/quiet/level2_new' # Not used (maybe externally)?
default_runlis  = '/data4/quiet/runlist_l2_new.txt'
#ORIG_RUNLIST    = '/data4/quiet/runlist_l2.txt' # Not used
BACKUP_RUNLIST  = '/data4/quiet/runlist_l2_new.txt.100113'
#BACKUP_RUNLIST  = '/data4/quiet/runlist_l2.txt.090922'
BACKUP2_RUNLIST = '/path/to/backup2_runlist_l2.txt'
#BACKUP2_RUNLIST = '/usit/titan/u1/jonathaz/quiet_svn/oslo/src/python/all.run'

#if platform.system() == 'Linux' :
#    try:
#        import slalib
#    except ImportError:
#        print "You need slalib"
#        sys.exit(1)
#else:
#    __name__='muesli'
#    help(__name__)
#    sys.exit(0)

#try:
#    import pyfits
#except ImportError:
#    print "You need modified pyfits and to be on Linux"
#    print "Get pyfits from /home/jtlz2/cobra/lib/"
#    sys.exit(1)


#def precision():

#    """Return machine precisions"""

#    from scipy import limits

#    limits.double_epsilon
#    limits.float_epsilon
#    limits.double_tiny
#    limits.float_tiny

#    return deps,feps

#-----------------------------------------------------------
# from http://code.activestate.com/recipes/66062
# use sys._getframe() -- it returns a frame object, whose attribute
# f_code is a code object, whose attribute co_name is the name:
#import sys
#this_function_name = sys._getframe().f_code.co_name
# the frame and code objects also offer other useful information:
#this_line_number = sys._getframe().f_lineno
#this_filename = sys._getframe().f_code.co_filename
# also, by calling sys._getframe(1), you can get this information
# for the *caller* of the current function.  So you can package
# this functionality up into your own handy functions:
def whoami():
    "Usage: me=whoami()"
    import sys
    return sys._getframe(1).f_code.co_name
# this uses argument 1, because the call to whoami is now frame 0.
# and similarly:
def callersname():
    "Usage: him = callersname()"
    import sys
    return sys._getframe(2).f_code.co_name

#-----------------------------------------------------------

def grep(pattern, *files):

    """
    Grep function from http://bytes.com/topic/python/answers/169012-grep
    """

    import re,glob

    result=[]
    search = re.compile(pattern).search
    for file in files:
        for index, line in enumerate(open(file)):
            if search(line):
                result.append((file, str(index+1), line[:-1]))
                #print ":".join((file, str(index+1), line[:-1]))

    return result

#-----------------------------------------------------------

def imported(mod):

    """
    Returns true (false) if a module is (not) already imported
    """

    boolean = False
    if mod in dir():
        boolean = True

    return boolean

#-----------------------------------------------------------

# From http://www.daniweb.com/code/snippet806.html
def find_key(dic, val):

    """return the key of dictionary dic given the value"""

    return [k for k, v in dic.iteritems() if v == val][0]

#-----------------------------------------------------------

def print_dict(dict):

    """
    print_dict is a function to print the contents of a dictionary.
    It takes a dictionary as its argument. For each key in the
    dictionary it prints the key followed by the value.
    """

    for key in dict:
        print key, dict[key]

    return

#-----------------------------------------------------------

def StripExtension(filename) :
    
    """Split a filename into a basename and an extension"""

    for i in range(2, 0, -1) :
        fparts = filename.rsplit('.', i)
        if len(fparts) == 3 :
            extension = fparts[1] + '.' + fparts[2]
        elif len(fparts) == 2 :
            extension = fparts[1]
        elif len(fparts) > 3:
            print 'More than two dots on which to split - exiting'
            sys.exit(1)

    basename=fparts[0]

    return basename,extension

#-----------------------------------------------------------

def ReadList(filename):
    """Reads an ASCII list of files into a list"""
    import string
    clist=[]
    file=open(filename)
    clist=file.readlines()
    file.close()
    # Remove trailing carriage return
    clist = [line[:-1] for line in clist]
    # Remove blank last line, if it's there
    if clist[-1]==' ':
        clist.remove(' ')
    # Remove initial comment lines, if there
    if string.count(clist[0],'#') > 0:
        del clist[0]
    if string.count(clist[0],'#') > 0:
        del clist[0]

    file.close()

    return clist


#-----------------------------------------------------------

def ReadList2(filename):
    """Reads an ASCII list of files into a list"""

    clist=[]
    file=open(filename)

    for line in file:
        line = line.split()
        size = len(line)
        for i in range (0,size):
            clist.append(line[i])
    file.seek(0,0)

    # Remove blank last line, if it's there
    if clist[-1]==' ':
        clist.remove(' ')

    file.close()

    return clist

#-----------------------------------------------------------

def ReadCutList(filename):

    """
    Read a cut list from file
    """

    import quiet

    return quiet.CutList(filename)

#-----------------------------------------------------------

def ReadAcceptlist(filename):

    """
    Reads an ASCII acceptlist into a quiet.py acceptlist
    """

    import quiet

    return quiet.AcceptedScanList(filename)

#-----------------------------------------------------------

def modulo(initial,target):

    """
    Subtract multiples of target from initial until initial is
    within +/- target of zero
    """
    #sign = lambda x:+(x > 0) or -(x < 0)
    #while True:
    #    if abs(initial) < target : return initial
    #    initial-=sign(initial)*target
    #return

    return initial % target

#-----------------------------------------------------------

def ParamDict(hin):

    "Return the dictionary of the random parameters"

    """
    This is required to connect field names like "c1", "c2", etc, to parameter names like uu, vv etc.
    """

    ## Look for header items starting with PTYPE -- they tell us what
    ## the random parameters are.
    pre=re.compile(r"PTYPE(?P<i>\d+)")
    
    res={}
    for k,v in hin.header.items():
        m=pre.match(k)
        if m :
            cname = "C%i" % int(m.group("i"))
            res[ cname ] = v.upper()
    return res

#-----------------------------------------------------------

def sepn2(r1,d1,r2,d2):

    """
    Calculate the separation between 2 sources, RA and Dec must be
    given in radians. Returns the separation in radians [TWS]
    """

    # NB slalib sla_dsep does this
    # www.starlink.rl.ac.uk/star/docs/sun67.htx/node72.html
    # pyephem also does this? (but what accuracy?)

    cos_sepn=sin(d1)*sin(d2) + cos(d1)*cos(d2)*cos(r1-r2)
    sepn = arccos(cos_sepn)
    
    return sepn

#-----------------------------------------------------------

def ReadCat(file,version='v2'):

    """
    Read source catalogue into dictionary by source ID number
    """

    try:
        if not imported('slalib'):
            import slalib
    except ImportError:
        print "You need slalib (on Linux systems)"
        sys.exit(1)

    catalogue={}
    
    lines=ReadList(file)
    #print 'LL',lines[0]

    for line in lines:
        words=re.split('\s+',line)
        del words[0]
        #print words

        id   = int(words[0])    #-- id
        rah  = float(words[1])  #}
        ram  = float(words[2])  #}- RA
        ras  = float(words[3])  #}
        decd = float(words[4])  #}
        decm = float(words[5])  #}- Dec
        decs = float(words[6])  #}
        if (version == 'v2'):
            flux = float(words[7])  #-- Flux
            alpha = float(words[8]) #-- Spectral index (cam)
        elif (version == 'v1'):
            flux = float(words[10]) #-- Flux
            alpha = 0.5             #-- Spectral index (cam)
        #date = float(words[8]) 

        (rarad,status)  = slalib.sla_dtf2r(rah,ram,ras) # hms2rad
        #rarad=rarad*15.0 # rad(h) -> rad(deg) => No matches!!
        (decrad,status) = slalib.sla_daf2r(decd,decm,decs)   # dms2rad

        if status != 0:
            print 'slalib status', status
            print 'Exiting...'
            sys.exit(status)

        catalogue[id]=[rarad,decrad,flux,alpha]

    #print_dict(catalogue)

    return catalogue

#-----------------------------------------------------------

def MachineEps():

    import numpy

    eps = numpy.finfo(numpy.double).eps

    return eps

#-----------------------------------------------------------

def QUIETParams():

    """
    Define QUIET parameters
    """

    return


#-----------------------------------------------------------

def ReadSidelobe(file):

    """
    Read a sidelobe file, dbi(u,v), and return as a cube
    e.g. /usit/titan/u1/jonathaz/sidelobes/fake/sidelobes_fake.txt
    """

    import os,scanf,numpy

    data=open(file,'r')

    # I cheated by hardcoding these..
    # Could avoid a second pass by putting at the start of the file
    umin=-1.0
    vmin=-1.0
    umax=1.0
    vmax=1.0
    du=0.02
    dv=0.02
    nu=int((umax-umin)/du)+1
    nv=int((vmax-vmin)/dv)+1
    print nu,nv
    arr=numpy.zeros((nu,nv))
    uu=numpy.zeros(nu)
    vv=numpy.zeros(nv)

    linef='%9f %9f %9f'
    for line in data:
        if line.startswith('#'): continue
        try:
            u,v,dbi=scanf.sscanf(line,linef)
        except scanf.IncompleteCaptureError:
            print 'Skipping %s' % line

        i=int(round(((u-umin)/du),0))
        j=int(round(((v-vmin)/dv),0))
        uu[i]=u
        vv[j]=v
        arr[i,j]=dbi
#        print i,j,u,v,uu[i],vv[i],dbi,arr[i,j]

    data.close()
    return arr,uu,vv

#-----------------------------------------------------------

def QuerySidelobe(arr,uarr,varr,u,v):

    """
    This is a very pedestrian bilinear interpolation routine
    Follows F code based on Numerical Recipes pp 116-7 (apparently)
    In this form, supply:
    arr = 2-D array of points
    uarr } values of indices of arr
    varr }
    u,v are the points you seek
    dbi is the interpolated value
    """

    import numpy

    umin=uarr[0]; umax=uarr[-1]
    vmin=varr[0]; vmax=varr[-1]
    nu=numpy.size(uarr)
    nv=numpy.size(varr)

    du=(umax-umin)/(nu-1) # -1 is a fix!
    dv=(vmax-vmin)/(nv-1) # -1 is a fix!

    # u = umin + i*du
    # i = u - umin / du
    i=int((u-umin)/du) # Deliberately rounds down
    j=int((v-vmin)/dv) #

#    print i,j
#    print u,v
#    print uarr[i],varr[j]
#    print uarr,varr

    dif_u=u-uarr[i] # Distance in u from point before
    dif_v=v-varr[j]

    t1 = dif_u/du  # Fractional distance from point before
    u1 = dif_v/dv

#    print du,dv,nu,nv,umin,umax,vmin,vmax
#    print i,j,uarr[i],varr[j],dif_u,dif_v
#    print t1,u1,uarr[i],uarr[i+1],varr[j],varr[j+1]

    assert(u >= umin), 'u out of range %f < %f' % (u,umin)
    assert(u <= umax), 'u out of range %f > %f' % (u,umax)
    assert(v >= vmin), 'v out of range %f < %f' % (v,vmin)
    assert(v <= vmax), 'v out of range %f > %f' % (v,vmax)

    dbi=-1.0e4
#    print arr[i,j],arr[i+1,j],arr[i,j+1],arr[i+1,j+1]
    if (i+1 < nu and j+1 < nv) and (i >= 0 and j >= 0):
        dbi = (1-t1)*(1-u1)*arr[i,j] \
              + t1*(1-u1)*arr[i+1,j] \
              + t1*u1*arr[i+1,j+1] \
              + (1-t1)*u1*arr[i,j+1]
    
    return dbi

#-----------------------------------------------------------

def QuerySidelobeCube(cube,du,dv,umin,umax,vmin,vmax,u,v):

    """
    Query gridded sidelobe data
    Return bilinearly-interpolated value
    """

    dbi=-99.0
        
    #print du,dv

    uu=int(u)
    vv=int(v)
 
    dif_u=u-uu
    dif_v=v-vv

    t1 = 1.0 - dif_u
    u1 = 1.0 - dif_v

    if (uu <= umax and vv <= vmax) and (uu >= umin and vv >= vmin):
        dbi=t1*u1*cube[uu][vv]\
             +dif_u*u1*cube[uu+du][vv]\
             +dif_u*dif_v*cube[uu+du][vv+dv]\
             +t1*dif_v*cube[uu][vv+dv]
        return dbi
    else:
        return None

    return

#-----------------------------------------------------------

def ReadParamDatabase(paramf,types=None):

    """
    Read parameter database file
    """

    import scanf

    if types is None:
        types = ['sigma0']

    sigma={}
    database=open(paramf,'r')

    linef = 'sigma0 %f %f %i %i %i %i %f'
    for line in database:
        #if not line.startswith(types[0]): pass
        try:
            mjdstart,mjdstop,scan,seg,mod,di,\
                     sigma[scan,seg,mod,di]=scanf.sscanf(line,linef)
        except scanf.IncompleteCaptureError:
            pass

    database.close()

    return sigma

#-----------------------------------------------------------

def ReadTypeBSigmas(indir=None,ces=None):

    """
    Read Type-B sigma values from file to memory
    Assumes CES translation has already been done on these files
    Optionally for a given CES (scan,seg)
    """

    import os,sys,scanf

    if indir is None:
        indir='/usit/titan/u1/jonathaz/quiet_data/calib_data/jonathaz/typeb_data/v2_translated'

    if not indir.endswith('translated'):
        print 'Are you sure these files are in Oslo format??'
        return -1

    if ces is None:
        infiles=os.listdir(indir)
        # Ignore files starting with '.' using list comprehension
        infiles=[filename for filename in infiles if \
             filename[0] != '.' and os.path.isdir(os.path.join(indir,filename))]
        inscans=[tuple([int(s) for s in infile.split('.')]) for infile in infiles]
        infiles.sort()
        inscans.sort()
    else:
        inscans=ces


    sigmas={}; preamp={}; mod=''
    linef='%2i %f %f %f %f'
    for (scan,seg) in inscans:
        mods=set() # Reset the module counter for each CES
        typebf = '%s/%i.%i/anomaly.txt' % (indir,scan,seg)
        if not os.path.exists(typebf):
            print 'typebf missing %s' % typebf 
        typeb=open(typebf,'r')
        for line in typeb:
            #if line.startswith('#') or line.startswith('__'): continue
            try:
                if mod not in mods:
                    mod,\
                          preamp[scan,seg,mod,0],\
                          preamp[scan,seg,mod,1],\
                          preamp[scan,seg,mod,2],\
                          preamp[scan,seg,mod,3],\
                          =scanf.sscanf(line,linef)
                    mods.add(mod)
                else:
                    mod,\
                          sigmas[scan,seg,mod,0],\
                          sigmas[scan,seg,mod,1],\
                          sigmas[scan,seg,mod,2],\
                          sigmas[scan,seg,mod,3],\
                          =scanf.sscanf(line,linef)
            except scanf.IncompleteCaptureError:
                pass

        typeb.close()

    noise={}
    factor=6880.0
    sampling_rate=100.0

    for (scan,seg) in inscans:
        noise[scan,seg,mod,di]=sqrt(factor)*sampling_rate*\
                                preamp[scan,seg,mod,di]*\
                                sigmas[scan,seg,mod,d]

    return preamp,sigmas,noise

#-----------------------------------------------------------

def QUIETBookends(band,format='UTC'):

    """
    Return QUIET observing periods (Q, W, [Q2, W2, ...])
    Optionally, return as MJD rather than UTC
    """

    if (band == 'Q'):
        UTCstart = '2008/10/22 20:30:00.0' # Started double-demod
        UTCstop  = '2009/06/13 13:52.00.0' # Laura suspended schedules to begin installing
                                           # the signal generator on top of the ground screen
    elif (band == 'W'):
        UTCstart = '2009/06/13 13:52.00.0' # Placeholder
        UTCstop  = '2010/12/31 23:59:00.0' # ???

    if format == 'MJD':
        MJDstart = UTC2MJD(UTCstart)
        MJDstop  = UTC2MJD(UTCstop)
        return MJDstart,MJDstop
    elif format == 'UTC':
        return UTCstart,UTCstop

    return

#-----------------------------------------------------------

def UTC2MJD(time_string):

    """
    Convert UTC YYYY/MM/DD HH:MM:SS.SS to MJD
    """

    import re,string,ephem

    #time_tuple = tuple(time.strptime(time_string, '%Y-%m-%d %I:%M:%S'))

    #if string.fnd('/') > n:

    if re.search('-', time_string) != None: # "if '-' in string:"
        time_string = string.replace(time_string,'-','/')

    date = ephem.Date(time_string)

    MJD = ephem.julian_date(date) - 2400000.5

    return MJD

#-----------------------------------------------------------

def mjd2djd(mjd):
    "Convert modified Julian date to Dublin date"
    return jd2djd(mjd,modified=True)

#-----------------------------------------------------------

def jd2mjd(jd): return jd - 2400000.5
def mjd2jd(mjd): return mjd + 2400000.5

#-----------------------------------------------------------

def jd2djd(JD,modified=False):

    """
    Convert Julian date to Dublin Julian date
    """

    if modified:
        JD += 2400000.5

    return JD - 2415020.0

#-----------------------------------------------------------

def getLengthOfMonth(key,doleap):
    """
    Return length of given month in that (leap) year
    """

    if doleap:
        month_map = {'jan': 31, 'feb': 29, 'mar': 31, 'apr': 30, 'may': 31, 'jun': 30, 'jul': 31, 'aug': 31, 'sep': 30, 'oct': 31, 'nov': 30, 'dec': 31}
    else:
        month_map = {'jan': 31, 'feb': 28, 'mar': 31, 'apr': 30, 'may': 31, 'jun': 30, 'jul': 31, 'aug': 31, 'sep': 30, 'oct': 31, 'nov': 30, 'dec': 31}
    key = key.lower()
    if month_map.has_key(key):
        return month_map[key]

    return 0

#-----------------------------------------------------------

def MJD2UTC(MJD):
    """
    """

    JD = MJD + 2400000.5
    UTC = (JD-2440587.5)*86400.0

    return UTC


#-----------------------------------------------------------

def OpenFile(filename,operation='w',permissions='urwgrwor'):

    """
    Open file for write
    Set permissions
    """
    
    import os, stat

    if os.path.exists(filename) and operation == 'w':
        os.rename(filename,filename+'.old')

    output=open(filename,operation)

    FixPerms(filename,permissions)

    return output

#-----------------------------------------------------------

def FixPerms(filename,permissions='urwgrwor'):

    """
    Set permissions; handle directories automatically
    """
    
    import os, stat

    dir=os.path.isdir(filename)

    if os.path.exists(filename):
        if permissions=='urwgrwor' and not dir:
            perms = stat.S_IREAD | stat.S_IWRITE | \
                    stat.S_IRGRP | stat.S_IWGRP  | \
                    stat.S_IROTH
        elif os.path.isdir(filename):
            perms = stat.S_IREAD | stat.S_IWRITE | \
                    stat.S_IRGRP | stat.S_IWGRP  | \
                    stat.S_IROTH | stat.S_IXUSR  | \
                    stat.S_IXGRP | stat.S_IXOTH  | \
                    stat.S_ISGID # Set group SUID bit
            # www.dartmouth.edu/~rc/help/faq/permissions.html
        else:
            print 'not sure how to handle permissions for this file'
            return 2

        os.chmod(filename,perms)
        return 0
    else:
        print filename, 'does not exist'
        return 1

#-----------------------------------------------------------


def WriteHeader(outfile,extras,style='cut',handle=None):

    """
    Write header for a .cut or jackknife file
    """

    import os, re, glob
    from datetime import date,datetime

    user = os.getenv('USER')
    host = os.getenv('HOSTNAME').split('.')[0]
    host = re.sub('login','titan',host) # Switch login -> titan

    #date = str(date.today())
    date = datetime.utcnow().strftime("%Y-%m-%d %A %I:%M%p UT")

    # Attempt to locate this source code
    thiscode=os.path.realpath(__file__)
    repos=os.path.dirname(thiscode)#+'/'
    scode=os.path.basename(thiscode)

    if os.path.exists(os.path.join(repos,scode)):
        sversion = grep('\$Id', *glob.glob(os.path.join(repos,scode)))[0][2]
    else:
        sversion = 'Where is your SVN? I tried '+os.path.join(repos,scode)
    qcode='quiet.py'
    if os.path.exists(os.path.join(repos,qcode)):
        qversion = grep('\$Id', *glob.glob(os.path.join(repos,qcode)))[0][2]
    else:
        qversion = 'Where is your SVN? I tried '+repos+qcode

    if not os.path.isabs(outfile):
        outfile=os.path.abspath(outfile)

    header = '# header_start\n'+ \
             '# This is ' +outfile+ '\n'+ \
             '# date: ' +date+ '\n'+ \
             '# run by ' +user+ '@' +host+ '\n'+ \
             '# SVN ' +sversion+ '\n' + \
             '# SVN ' +qversion+ '\n' + \
             '# patch: patch_6a\n'+ \
             '# info: weather, static, dead\n'+ \
             '# header_end\n'
       #      '# patch : ' +extras+ '\n' \

    if style=='cut': handle.write(header)

    return header

#-----------------------------------------------------------

def WriteCutLine(handle,mjdstart,mjdstop,ces,patch,modulelist,\
                 diodelist,comment):

    """
    Write a single cut line to file
    """

    cmd = 'cut'

    #line = #cmd+\
    line = '%s %14.8f %14.8f %s %s %s %s # %s' \
           %(cmd,mjdstart,mjdstop,str(ces),patch,\
             modulelist,diodelist,comment)
           #str(mjdstart) +' '+ \
           #+'%12.5f ' %(mjdstop)\
           #str(mjdstop) +' '+ \
           #+str(ces) +' '+ \
           #patch +' '+ \
           #modulelist +' '+ \
           #diodelist +' '+ \
           #'# '+comment

    handle.write(line+'\n')

    return

#-----------------------------------------------------------

def GenerateModuleList(band):

    """
    Populate the list of modules given the band
    """

    modulelist=[]
    if (band=='Q'):
        modulelist = range(19)
    elif (band=='W'):
        modulelist = range(91)

    return modulelist

#-----------------------------------------------------------

def SnipModule(cutmodule,modulelist):

    """
    Cut a module out of the active module list
    """

    modulelist.remove(cutmodule)

    return modulelist

#-----------------------------------------------------------

def ModuleStringToList(modulestring):

    """
    Convert concise module string from a .cut entry into a list
    """

    import re,string

    modulelist=[]
    if modulestring == 'all':
        modulelist=range(19)
    elif modulestring == 'none':
        modulelist=[]
    else:
        # Parse the string
        blocks=modulestring.split(',')
        print blocks
        for block in blocks:
            if re.search('-',block) == None:
                modulelist.append(int(block))
            else:
                ends=block.split('-')
                print ends
                for end in range(int(ends[0]),int(ends[1])+1):
                    modulelist.append(int(end))

    modulelist.sort()

    return modulelist

#-----------------------------------------------------------

def ModuleListToString(modulelist=range(19)):

    """
    Convert module list into concise format for a .cut entry
    Examples: modulelist=[0,1,2,3,5,6,8,14,15,18]
              modulelist=[0,2,4,6,8,14,15,17]
              [0,1,3,5,6,7,11,12,13,14,15,18]
              [0,1,3,5,6,7,11,12,   14,15,18]
    """

    print modulelist

    modulelist.sort()   # } Sort elements in situ...
    import sets         # } ... and remove duplicates
    modulelist=list(sets.Set(modulelist))

    #print modulelist

    if not modulelist: return 'none'

    # First element will always appear
    string=str(modulelist[0])

    # For just one element in list
    if len(modulelist) < 2: return string

    # For just two elements i list, return them
    if len(modulelist) == 2:
        if modulelist[0] == modulelist[1] - 1:
            string+='-'
        else:
            string+=','
        string+=str(modulelist[-1])
        return string

    # For all other cases (> 2 elements)
    prev = modulelist[0]
    if modulelist[1] == modulelist[0]+1:
        string+='-'
    else:
        string+=','

    for position,module in enumerate(modulelist[:-1]):
        prev = modulelist[position-1]
        next = modulelist[position+1]
        if position > 0: # Skip first element, but need it to hand
            #print position,prev,module,next
            if module == prev+1 and module == next-1:
                # Mid-run
                pass
            if module == prev+1 and module != next-1:
                # End of run; print value and terminate run with comma
                string += str(module)+','
            if module != prev+1 and module == next-1:
                # Start of new run, so print value with dash
                string += str(module)+'-'
            if module != prev+1 and module != next-1:
                # Lone value - print it and terminate with comma
                string += str(module)+','

    string+=str(modulelist[-1]) # Always append last element

    return string

#-----------------------------------------------------------

# Define CMB patches to be used for null tests
CMBPatches=['patch_2a','patch_4a','patch_6a','patch_7b']
#CMBPatches=['patch_2a','patch_4a','patch_6a','patch_7a','patch_7b']
AllPatches=['patch_21','patch_2a','patch_3','patch_4','patch_4a',\
            'patch_6','patch_6a','patch_7','patch_7a','patch_7b',\
            'patch_a','patch_b','patch_ga','patch_gb','patch_gc',\
            'patch_ggc']

#-----------------------------------------------------------

# Define ordering of diodes
Diodes=['Q1','U1','U2','Q2']
DiodesLookup={'Q1':0,'U1':1,'U2':2,'Q2':3}
Q1='Q1'; U1='U1'; U2='U2'; Q2='Q2'

#-----------------------------------------------------------

def ReadL2Runlist(level=2,runfile=None):

    """
    Read a [L1/]L2 runlist
    """

    import os
    import quiet

    if runfile is None: runfile=default_runlis

    # Check L2 runlist exists (disk quota problems etc)
    # JTLZ This logic tree is horrible.. need a good way to simplify
    if os.path.exists(runfile):
        print 'reading runlist', runfile
        if level == 1:
            runlist = quiet.L1Runlist(runfile)
        elif level == 2:
            runlist = quiet.L2Runlist(runfile)
        # Assume you want all the CMB patches:
        targetsfetched=runlist.targetTypes['cmb'].targets.keys()
        if set(CMBPatches).issubset(set(targetsfetched)):
            return runlist
        else:
            print 'problem with runlist', runfile
            # Try backup L2 runlist
            runfile = BACKUP_RUNLIST
            if os.path.exists(runfile):
                print 'reading runlist', runfile
                if level == 1:
                    runlist = quiet.L1Runlist(runfile)
                elif level == 2:
                    runlist = quiet.L2Runlist(runfile)
                # Assume you want all the CMB patches:
                targetsfetched=runlist.targetTypes['cmb'].targets.keys()
                if set(CMBPatches).issubset(set(targetsfetched)):
                    return runlist
                else:
                    print 'problem also with runlist', runfile
                    print runlist.targetTypes['cmb'].targets.keys()
                    print "You have other options, but I'm giving up for now"
                    print 'Exiting...'
                    sys.exit(2)
    else:
    # [Should really now try default L2 list, if it's not that already]
    # Try backup L2 runlist
        print "Can't find %s" % runfile
        runfile = BACKUP_RUNLIST
        if os.path.exists(runfile):
            print 'reading runlist', runfile
            if level == 1:
                runlist = quiet.L1Runlist(runfile)
            elif level == 2:
                runlist = quiet.L2Runlist(runfile)
            return runlist
        else:
    # Try cut-down (CMB only) L2 runlist
            runfile = BACKUP2_RUNLIST
            if os.path.exists(runfile):
                print 'reading runlist', runfile
                if level == 1:
                    runlist = quiet.L1Runlist(runfile)
                elif level == 2:
                    runlist = quiet.L2Runlist(runfile)
                return runlist
            else:
                print 'Could not find a decent runlist! Exiting...'
                sys.exit(1)

    return None

#-----------------------------------------------------------

def CleanL2Runlist(runlist,*keep):

    """
    Delete targets from L2 runlist, except for those specified
    
    JAZ: Fixed a bug that meant the default "keep" would not work
         and also edited so that you can now call in three ways:
    CleanL2Runlist(runlist)
    CleanL2Runlist(runlist,"patch_gc","patch_2a")
    myListOfPatches=["patch_gc","patch_2a"]
    CleanL2Runlist(runlist,myListOfPatches)
    """

    import quiet

    if not keep: keep=CMBPatches  #If no arguments are passed use the default list  #JAZ
    elif len(keep)==1 and type(keep[0])==list: keep=keep[0]   #if only one arg is passed and it's a list, use that #JAZ
    #otherwise use any passed #JAZ

    # Remove non-CMB target types
    for target in runlist.targetTypes.keys(): # .keys() needed!
        if target != "cmb":
            runlist.removeTargetType(target)

    # Remove non-main patches
    for target in runlist.targetTypes["cmb"].targets.keys(): # .keys() needed!
        if target not in keep:   #JAZ edited from [keep]
            runlist.targetTypes["cmb"].removeTarget(target)

    #print runlist.targets()
    
    return runlist

#-----------------------------------------------------------

def Test():

    """
    Dumping ground for useful commands etc
    """

    # Get list of run ids
    #print runlist.targetTypes['cmb'].targets['patch_2a'].scans.keys()

    # Get at segments for run 1041
    #runlist.targetTypes['cmb'].targets['patch_2a'].scans[1041].segments

    return

#-----------------------------------------------------------


def WriteL2Runlist(runlist,patch,outfile='out.run'):

    """
    Write L2 runlist to file for specified patch(es)
    """

    import quiet

    #if patch == 'all':
    #    patch = CMBPatches

    # Should now be done already, but repeat anyway..
    runlist = CleanL2Runlist(runlist,patch)

    runlist.toFile(outfile)

    return

#-----------------------------------------------------------

def passover(acceptlist,bookofsaints):

    """
    Remove every CES from an acceptlist, except those in the
    bookofsaints
    """

    for scan,seg in acceptlist.scans.keys():
        cesname=str(scan)+'.'+str(seg)
        if (scan,seg not in bookofsaints.scans.keys()):
            zap(acceptlist,cesname)

    return acceptlist

#-----------------------------------------------------------

def multizap(acceptlist,removelist):

    """
    Remove one removelist of CESs from an acceptlist
    """

    for scan,seg in removelist.scans.keys():
        cesname=str(scan)+'.'+str(seg)
        zap(acceptlist,cesname)

    return acceptlist

#-----------------------------------------------------------

def zap(acceptlist,ces,convention='oslo'):

    """
    Delete one CES from acceptlist, IFF it's in the acceptlist already
    """

    import quiet, re

    scan,seg = re.split('\.',ces)
    scan,seg = int(scan),int(seg)

    # Check CES is in acceptlist:
    if acceptlist.scans.has_key((scan,seg)):
        acceptlist.removeScan(scan,seg)

    return acceptlist

#-----------------------------------------------------------

def acopy(acceptlist):

    """
    Make a deep copy of an acceptlist
    """

    import copy

    return copy.deepcopy(acceptlist)

#-----------------------------------------------------------

def acount(acceptlist):

    """
    Count number of CESs in an acceptlist or dictionary of acceptlists
    """

    count={}
    if type(acceptlist)==type(dict()):
        for patch in acceptlist:
            count[patch] = len(acceptlist[patch].scans.keys())
        return count
    else:
        return acceptlist.acount()

#-----------------------------------------------------------

def GenerateCMBRunlists(patch=None,silent=True,runfile=None):

    """
    Write (and return) L2 runlists for the CMB patches (specified)
    silent argument switches off file write
    Runlists and acceptlists are returned in dictionaries keyed by patch
    """

    import quiet
    import copy

    rundict={}
    acceptdict={}

    level=2

    if patch is None:
        patches=CMBPatches[:]
    elif type(patch) == list:
        patches=patch[:]
    else:
        patches=[]
        patches.append(patch)

    runlist = ReadL2Runlist(level,runfile=runfile)
    for patch in patches:
        runlistcopy = copy.deepcopy(runlist)
        runlistcopy = CleanL2Runlist(runlistcopy,patch)

        if not silent:
            outfile = patch + '.run'
            if os.path.exists(outfile):
                os.rename(outfile,outfile+'.old')
            # ** NB CleanL2Runlist is repeated within WriteL2Runlist
            WriteL2Runlist(runlistcopy,patch,outfile)

        rundict[patch]    = runlistcopy
        acceptdict[patch] = quiet.acceptedListFromRunList(runlistcopy)

    return rundict, acceptdict

#-----------------------------------------------------------

def MatchLevels(L2_rundict,L1_runlist):

    """
    Take L2 runlist and match the associated L1 files from L1 runlist
    """

    import quiet

    mapping={}
    # For each patch
    for patch in L2_rundict.keys():
        # Pull the runlist out from the runlist dictionary
        runlist=L2_rundict[patch].targetTypes['cmb'].targets[patch]
        level1run=L1_runlist.targets[patch]

        # For each scan
        for scan in runlist.scans.keys():
            #print scan, runlist.scans[scan]
            # For each segment in the scan
            for nseg in range(0,runlist.scans[scan].nSegment):
                seg=runlist.scans[scan].segments[nseg].id
                # Construct CES
                ces='.'.join([str(s) for s in [scan,seg]])
                # Identify L1 file for this CES
                # NOT SURE this is quite right??
                mapping[ces]=level1run.scans[scan].files[nseg]
                #  print runlist.targetTypes[type].targets[patch].scans[scan]

    return mapping

#-----------------------------------------------------------

def MJDBisect(runlist,patch,type='cmb'):

    "Find the MJD midpoint of the patch CESs in a runlist"

    import quiet

    # Need to find patch start and end MJDs (i.e. real patch bookends)
    starttime = 0.0
    stoptime = 1.0

    #scanlengths={}
    scanlengths = 0.0
    integrator = 0.0
    for scan,seg in sorted(runlist[patch].scans.keys()):
        t0 = runlist.targetTypes[type].targets[patch].scans[scan].startDate
        t1 = runlist.targetTypes[type].targets[patch].scans[scan].endDate
        #scanlength = t1-t0
        integrator += t1-t0

    halftime = integrator/2.0

    # Find MJD of patch bisector
    bisectorMJD = starttime + halftime

    # Need to find bisectorID

    # Summarize:
    print 'Patch', patch
    print '  start MJD =', starttime
    print '  stop MJD =', stoptime
    print '  total scanning time/hours =', integrator/24.0

    print '  midpoint is at MJD =', bisectorMJD
    print '  which falls in scan', bisectorID

    return bisectorID, bisectorMJD, scanlengths

#-----------------------------------------------------------


def write_cmd(commands,filename):
    """Write an ASCII cmd file from a list"""

    import os

    if os.path.exists(filename):
        os.rename(filename,filename+'.old')

    output=open(filename,'w')
    for line in commands:
        output.write(line+'\n')

    output.close()

    return

#-----------------------------------------------------------

def TuneDiodeThreshold(acceptlist=None,postprocess=True,threshold=None):

    """
    For an acceptlist, allow tuning of the min_diode_accept so that
    nmin-nmax CESs are left
    Report the corresponding min_diode_accepts
    """

    import quiet

    # Apply the tod2map cuts to the Columbia list
    if acceptlist is None: # Default
        cutsdir='/data4/quiet/nulltests/cuts'
        colf='mlaccept-100117a.txt'
        acol=quiet.AcceptedScanList(os.path.join(cutsdir,colf))
        dir='/data4/quiet/jonathaz4/maps_columbia_100117a/cmb/patch_6a'
        tod2mapf='accepted_6a_columbia_100117a.sort'
        aout=quiet.AcceptedScanList(os.path.join(dir,tod2mapf))
        if postprocess:
            removed,acol=PostProcessAcceptlists(acceptin=acol,acceptout=aout)

        acol.printSummary()

        nbefore=len(acol.scans.keys())

        if threshold is None: threshold=0.6
        acol=ApplyMinDiodeAccept(threshold=threshold,acceptlist=acol,\
                            verbose=True,patch='patch_6a')

        nafter=len(acol.scans.keys())
        nlost=nbefore-nafter

        print 'After threshold %1.1f, %i are lost' % (threshold,nlost)
        acol.printSummary()

    return

#-----------------------------------------------------------

def DifferenceAcceptlistsToMatrix(a1,a2,outf='diff.png',master=None):

    """
    For the union of scans in two acceptlists (or for a master),
    plot those two acceptlists to a matrix
    If a1 and a2 aren't acceptlists, assume they are files thereof
    """

    import quiet
    import os
    import pylab,numpy

    #assert(master is not None), 'master is required for %s!' % whoami()

    if type(a1) != quiet.AcceptedScanList:
        a1=quiet.AcceptedScanList(a1)
    if type(a2) != quiet.AcceptedScanList:
        a2=quiet.AcceptedScanList(a2)

    #nmodules=quiet.N_MODULE_MAX-2
    ndiodes=quiet.N_MODULE_MAX*quiet.NDIODE_MODULE

    s1=set(a1.scans.keys())
    s2=set(a2.scans.keys())

    scanlist=[]
    if master is None:
        either  = s1 | s2  # Union
        scanlist=sorted(list(either))
        print "Can't yet go beyond here..."
        return
    else:
        scanlist=sorted(master.scans.keys())

    nces=len(scanlist)

    Z=numpy.zeros((ndiodes,nces))

    # Set values (bearing in mind colourmap)
    #map=pylab.cm.jet
    #map=pylab.cm.gray
    #map=pylab.cm.spectral
    map=pylab.cm.hot
    if map==pylab.cm.jet:
        inclval=0.5
        exclval=0.0
        absentval=1.0
    elif map==pylab.cm.hot:
        inclval=1.0
        exclval=0.0
        absentval=0.0

    agreepositive=1.0
    agreenegative=-1.0
    #inneither=0.0
    a1cutonly=-0.5 # red
    a2cutonly=0.5  # yellow

    # Initialize counters
    counters={'a1cutonly':0,'a2cutonly':0,'agreekeep':0,'agreecut':0}

    mods=range(quiet.N_MODULE_MAX) # This is designed to speed up the for loops
    diodes=range(quiet.NDIODE_MODULE)

    # Run through each CES in the scanlist
    for ices,(scan,seg) in enumerate(scanlist):
        idiode=0
        for module in mods:
            for diode in diodes:
                # If CES is in the 
                if (scan,seg) in s1 and (scan,seg) in s2:
                    # If the diode is included in a1
                    ina1=a1.isIncluded(scan,seg,module,diode)
                    ina2=a2.isIncluded(scan,seg,module,diode)
                    if ina1 and ina2:
                        Z[idiode,ices]=agreepositive
                        counters['agreekeep']+=1
                    elif not ina1 and not ina2:
                        Z[idiode,ices]=agreenegative
                        counters['agreecut']+=1
                    elif ina1 and not ina2:
                        Z[idiode,ices]=a2cutonly
                        counters['a2cutonly']+=1
                    elif not ina1 and ina2:
                        Z[idiode,ices]=a1cutonly
                        counters['a1cutonly']+=1
                elif (scan,seg) in s1 and (scan,seg) not in s2:
                    ina1=a1.isIncluded(scan,seg,module,diode)
                    if ina1:
                        Z[idiode,ices]=a2cutonly
                        counters['a2cutonly']+=1
                    elif not ina1:
                        Z[idiode,ices]=agreenegative
                        counters['agreecut']+=1
                elif (scan,seg) not in s1 and (scan,seg) in s2:
                    ina2=a2.isIncluded(scan,seg,module,diode)
                    if ina2:
                        Z[idiode,ices]=a1cutonly
                        counters['a1cutonly']+=1
                    elif not ina2:
                        Z[idiode,ices]=agreenegative
                        counters['agreecut']+=1
                elif (scan,seg) not in s1 and (scan,seg) not in s2:
                    Z[idiode,ices]=agreenegative
                    counters['agreecut']+=1
                else:
                    print "This shouldn't happen!"
                    return 1
                idiode+=1

    scannames=[]
    for (scan,seg) in scanlist:
        scannames.append('%s.%s' % (str(scan),str(seg)))

    #counters={'a1cutonly':(Z==getattr(whoami(),'a1cutonly').sum())}
    counters2={'a1cutonly':(Z==a1cutonly).sum(),\
               'a2cutonly':(Z==a2cutonly).sum(),\
               'agreekeep':(Z==agreepositive).sum(),\
               'agreecut':(Z==agreenegative).sum()}
    assert(counters==counters2), '%s %s' % (counters,counters2)
    print 'counters: %s (sum %i)' % (str(counters), sum(counters.values()))

    # Calculate 1-D stats, by diode and by CES
    stats={}
    for axis in range(Z.ndim):
        stats[axis]={'a1cutonly':(Z==a1cutonly).sum(axis=axis),\
                     'a2cutonly':(Z==a2cutonly).sum(axis=axis),\
                     'agreekeep':(Z==agreepositive).sum(axis=axis),\
                     'agreecut':(Z==agreenegative).sum(axis=axis)}

    outroot='.'.join(outf.split('.')[:-1])
    # stats by diode
    statsfd='%s.bydiode'%outroot
    if os.path.exists(statsfd): os.rename(statsfd,'%s.old'%statsfd)
    statsd=open(statsfd,'w')
    statsd.write('# mod di keep cut a1cut a2cut\n')
    statsd.write('# a1=%s\n' % a1.title)
    statsd.write('# a2=%s\n' % a2.title)
    ax=1
    for mod in mods:
        for di in diodes:
            statsd.write('%2i %1i %3i %3i %3i %3i\n' %(mod,di,\
                     stats[ax]['agreekeep'][4*mod+di],\
                     stats[ax]['agreecut'][4*mod+di],\
                     stats[ax]['a1cutonly'][4*mod+di],\
                     stats[ax]['a2cutonly'][4*mod+di]))
    statsd.close()
    print '\nLook in %s for stats by diode' % statsfd


    # stats by CES
    statsfc='%s.byces'%outroot
    if os.path.exists(statsfc): os.rename(statsfc,'%s.old'%statsfc)
    statsc=open(statsfc,'w')
    statsc.write('# run seg keep cut a1cut a2cut\n')
    statsc.write('# a1=%s\n' % a1.title)
    statsc.write('# a2=%s\n' % a2.title)
    ax=0
    for n,(scan,seg) in enumerate(scanlist):
        statsc.write('%4i %1i %3i %3i %3i %3i\n' %(scan,seg,\
                     stats[ax]['agreekeep'][n],\
                     stats[ax]['agreecut'][n],\
                     stats[ax]['a1cutonly'][n],\
                     stats[ax]['a2cutonly'][n]))
    statsc.close()
    print 'Look in %s for stats by CES' % statsfc

    caption='%s (%i red) against %s (%i yellow) [%i white, %i black]' % \
             (a1.title,counters['a1cutonly'],a2.title,counters['a2cutonly'],\
              counters['agreekeep'],counters['agreecut'])


    #pylab.pcolor(Z,cmap=pylab.cm.gray_r)
    # See www.scipy.org/Cookbook/Matplotlib/Show_colormaps
    pylab.pcolor(Z,cmap=map)
    cmin=-1.0
    cmax=1.0
    pylab.clim(cmin, cmax)
    #pylab.pcolor(Z)
    pylab.colorbar()
    plotwidth=200.0
    plotheight=8.0
    pylab.gcf().set_size_inches(plotwidth,plotheight)
    pylab.xlabel('CES')
    pylab.xlim([0.0,nces])
    xlabelskip=int(200/int(plotwidth)) # i.e. 5 at 40in, 1 at 200in
    #xticks=range(0,nces+xlabelskip,xlabelskip)
    xticks=[x+0.5 for x in range(0,nces+xlabelskip,xlabelskip)]
    labels=scannames[0:len(scannames):xlabelskip]
    pylab.xticks(xticks,labels,rotation=80)
    pylab.ylabel('diode number')
    pylab.ylim([0.0,ndiodes])
    yticks=range(0,ndiodes+quiet.NDIODE_MODULE,quiet.NDIODE_MODULE)
    pylab.yticks(yticks)
    pylab.title(caption)
    #pylab.set_size_inches(6.0,4.0)
    pylab.savefig(outf,format=outf.split('.')[-1])
    pylab.close()

    print 'acceptlist plotted to %s' % outf
    print 'white  - agree included (%2.1f for %i diodes)' % \
           (agreepositive,counters['agreekeep'])
    print 'black  - agree excluded (%2.1f for %i diodes)' % \
           (agreenegative,counters['agreecut'])
    print 'red    - %s-cut only (%2.1f for %i diodes)' % \
          (os.path.basename(a1.title),a1cutonly,counters['a1cutonly'])
    print 'yellow - %s-cut only (%2.1f for %i diodes)' % \
          (os.path.basename(a2.title),a2cutonly,counters['a2cutonly'])

    return Z

#-----------------------------------------------------------

def AcceptlistToMatrix(acceptlist,outf='test.png',master=None):

    """
    Take an acceptlist, and plot it to a matrix
    Optionally do this for scans in a master acceptlist rather than
      for the acceptlist itself
    """

    import quiet

    import pylab,numpy

    #nmodules=quiet.N_MODULE_MAX-2
    ndiodes=quiet.N_MODULE_MAX*quiet.NDIODE_MODULE

    scanlist=[]
    if master is None:
        scanlist=sorted(acceptlist.scans.keys())
    else:
        scanlist=sorted(master.scans.keys())
    nces=len(scanlist)

    Z=numpy.zeros((ndiodes,nces))

    sa=set(acceptlist.scans.keys())

    # Set values (bearing in mind colourmap)
    #map=pylab.cm.jet
    #map=pylab.cm.gray
    #map=pylab.cm.spectral
    map=pylab.cm.hot
    if map==pylab.cm.jet:
        inclval=0.5
        exclval=0.0
        absentval=1.0
    elif map==pylab.cm.hot:
        inclval=1.0
        exclval=0.0
        absentval=0.0

    # Run through each CES in the scanlist
    for ices,(scan,seg) in enumerate(scanlist):
        idiode=0
        for module in range(quiet.N_MODULE_MAX):
            for diode in range(quiet.NDIODE_MODULE):
                # If CES is in the master
                if (scan,seg) in sa:
                    # If the diode is included
                    if acceptlist.isIncluded(scan,seg,module,diode):
                        Z[idiode,ices]=inclval
                    # Diode not included
                    else:
                        Z[idiode,ices]=exclval
                # CES not in master
                else:
                    Z[idiode,ices]=absentval
                idiode+=1

    # Colour map:
    # 1.0 
    # 0.5 
    # 0.0 

    scannames=[]
    for (scan,seg) in scanlist:
        scannames.append('%s.%s' % (str(scan),str(seg)))

    #pylab.pcolor(Z,cmap=pylab.cm.gray_r)
    # See www.scipy.org/Cookbook/Matplotlib/Show_colormaps
    pylab.pcolor(Z,cmap=map)
    #pylab.pcolor(Z)
    #pylab.colorbar()
    plotwidth=40.0
    plotheight=8.0
    pylab.gcf().set_size_inches(plotwidth,plotheight)
    pylab.xlabel('CES')
    pylab.xlim([0.0,nces])
    xlabelskip=5
    xticks=range(0,nces+xlabelskip,xlabelskip)
    labels=scannames[0:len(scannames):xlabelskip]
    pylab.xticks(xticks,labels,rotation=80)
    pylab.ylabel('diode number')
    pylab.ylim([0.0,ndiodes])
    yticks=range(0,ndiodes+quiet.NDIODE_MODULE,quiet.NDIODE_MODULE)
    pylab.yticks(yticks)
    pylab.title(acceptlist.title)
    #pylab.set_size_inches(6.0,4.0)
    pylab.savefig(outf,format=outf.split('.')[-1])
    pylab.close()
    print 'acceptlist plotted to %s' % outf
    print 'black/red   - diode is included'
    print 'grey/green  - diode is excluded'
    if master is not None:
        print 'white/blue  - diode not in master list (whole scan cut)'

    return

#-----------------------------------------------------------

def ExtractPatchFromAcceptlist(acceptlist,master=None,patch=None):

    """
    For an acceptlist (e.g. Chicago acceptlist),
    strip out all CESs not from the given patch,
    by querying a master acceptlist dictionary
    """

    import copy

    assert(patch is not None), 'Patch is required!'

    # Build a dictionary of patches keyed by (scan,seg)
    if master is None:
        run,acceptdict=GenerateCMBRunlists()
    else:
        acceptdict=master

    acceptlistforpatch=copy.deepcopy(acceptlist)
    backup=copy.deepcopy(acceptlist)

    patchset=set(acceptdict[patch].scans.keys())
    for (scan,seg) in sorted(acceptlist.scans.keys()):
        if (scan,seg) not in patchset:
            acceptlistforpatch.removeScan(scan,seg)

    assert(acceptlist==backup), 'acceptlist has been modified!'

    print 'Extracted CESs for patch %s' % patch[-2:]

    return acceptlistforpatch

#-----------------------------------------------------------

def WrapperSplitAcceptlistIntoPatches(acceptf='accept_chicago_v8.txt',\
        outdir='/data4/quiet/nulltests/cuts/chicago8/'):

    """
    Take a Chicago (but Oslo-numbered) acceptlist with several patches
    split them by CMB patch
    and write them to file, one for each patch
    """

    import os,quiet

    #outdir='/data4/quiet/nulltests/cuts/chicago8/'
    #acceptf='accept_chicago_v8.txt'
#    outdir='/usit/titan/u1/roberdum/get_cuts/'
#    acceptf='accept_chicago_v15i.txt'

    chiclist=quiet.AcceptedScanList(os.path.join(outdir,acceptf))

    chicdict=SplitAcceptlistIntoPatches(chiclist)

    for patch in sorted(chicdict.keys()):
        patchf='accept_chic8_%s.txt' % patch[-2:]
        chicdict[patch].toFile(os.path.join(outdir,patchf))
        print 'wrote %s' % patchf

    return

#-----------------------------------------------------------

def SplitAcceptlistIntoPatches(acceptlist):

    """
    Take an acceptlist (already in Oslo numbering) for several patches
    and split it into a dictionary of acceptlists keyed by CMB patch
    """

    import copy

    run,masterdict=GenerateCMBRunlists()

    backup=copy.deepcopy(acceptlist)
    acceptdict={}
    for patch in sorted(masterdict.keys()):
        acceptdict[patch]=\
           ExtractPatchFromAcceptlist(acceptlist,master=masterdict,patch=patch)

    assert(acceptlist==backup), 'acceptlist has been modified!'

    return acceptdict

#-----------------------------------------------------------

def ConvertChicagoFilesToAcceptlists(files=None,master=None,patches=None):

    """
    Take a list of Chicago SQL files (preprocessed for Oslo by processSQL.py)
    and write them out to acceptlists, for each patch specified
    """

    import os,quiet

    if files is None or os.path.isdir(files):
        if os.path.isdir(files):
            chicagodir=files
        else:
            chicagodir='/data4/quiet/nulltests/cuts/chicago/'
        files=os.listdir(chicagodir)
        for f in files[:]:
            if not f.endswith('_oslo.txt'):
                files.remove(f)

    if patches is None:
        patches=CMBPatches

    run,acc=GenerateCMBRunlists()

    achic={}
    for patch in patches:
        for f in files:
            achic[(f,patch)]=ConvertChicagoListToAcceptlist(\
                f,master=acc[patch],patch=patch)
            fout='%s_%s.acc' % (StripExtension(\
                os.path.basename(f))[0],patch.split('_')[-1])
            achic[(f,patch)].toFile(fout)
            print '%s - %s --> %s' % (patch,f,fout)

    return

#-----------------------------------------------------------

def ConvertChicagoListToAcceptlist(chicf,master=None,patch=None):

    """
    Read a Chicago list of accepted diodes, and
    return it as an acceptlist for a given patch
    """

    import copy
    import quiet

    chicdata=open(chicf,'r')

    if patch is None:
        if master is None:
            print 'Need a master or patch to generate an acceptlist!'
        else:
            master=master
    elif patch is not None:
        if master is None:
            # Regenerate the acceptlist for this patch (slowish)
            rundict,acceptdict=GenerateCMBRunlists()
            master=acceptdict[patch]
        else:
            master=master

    chicagoaccepted=set()
    for line in chicdata:
        if line.startswith('#'): continue
        scan,seg,module,diode=line.split()
        chicagoaccepted.add((int(scan),int(seg),int(module),int(diode)))

    acceptlist=copy.deepcopy(master)

    acceptlist.printSummary()

    for (scan,seg) in sorted(master.scans.keys()):
        for module in range(quiet.N_MODULE_MAX):
            for diode in range(quiet.NDIODE_MODULE):
                # Test for (fast!) set membership
                if (scan,seg,module,diode) not in chicagoaccepted:
                    #print 'killing %i %i %i %i' % (scan,seg,module,diode)
                    acceptlist.excludeSegmentDiode(scan,seg,module,diode)

    acceptlist.removeEmptyScans()

    acceptlist.title=chicf

    acceptlist.printSummary()

    return acceptlist

#-----------------------------------------------------------

def CompareAcceptlistFiles(chicf=None,colf=None,
                           postprocess=False,t2mf=None,patch=None):

    """
    Compare two acceptlist files and report findings
    """

    import os
    import quiet

    if patch is None:
        patch='patch_6a'

    if chicf is None:
        cutsdir='/data4/quiet/nulltests/cuts'
        chicf='acc_chic_q_v6_pol_translat.txt'
        achic=quiet.AcceptedScanList(os.path.join(cutsdir,chicf))
    else:
        achic=quiet.AcceptedScanList(chicf)

    if colf is None:
        cutsdir='/data4/quiet/nulltests/cuts'
        colf='mlaccept-100117a.txt'
        acol=quiet.AcceptedScanList(os.path.join(cutsdir,colf))
    else:
        acol=quiet.AcceptedScanList(colf)

    # Apply the tod2map cuts to the Columbia list
    if postprocess:
        if t2mf is None:
            dir='/data4/quiet/jonathaz4/maps_columbia_100117a/cmb/patch_6a'
            t2mf='accepted_6a_columbia_100117a.sort'
            aout=quiet.AcceptedScanList(os.path.join(dir,t2mf))
        else:
            aout=quiet.AcceptedScanList(t2mf)

        removed,acol=PostProcessAcceptlists(acceptin=acol,acceptout=aout)

    common,neither,chiconly,colonly,chic,col,patchset,\
         =CompareAcceptlists(achic,acol,patch)

    print '%s: %i common, %i (chicago only) %i (columbia only), %i neither' \
          % (patch,len(common),len(chiconly),len(colonly),len(neither))
    print '%s total %i, %i (chicago) %i (columbia)' \
          % (patch,len(patchset),len(chic),len(col))

    outf='chiconly.txt'
    f=open(outf,'w')
    for (scan,seg,module,diode) in chiconly:
        line = '%i %i %i %i\n' %(scan,seg,module,diode)
        f.write(line)
    f.close()

    return removed,achic,acol

#-----------------------------------------------------------

def CompareAcceptlists(a1,a2,patch='patch_6a',acceptlist=None):

    """
    Compare two acceptlists as robustly as possible,
    given a patch
    Reference to a third superset acceptlist (assumed to be all included)
    """

    import copy
    import secateur,quiet

    if acceptlist is None:
        rundict,acceptdict=GenerateCMBRunlists()
        acceptlist=acceptdict[patch]

    masterlist=[]
    masterexlist=[]
    # Turn acceptlist into a set
    for (scan,seg) in acceptlist.scans.keys():
        for module in range(quiet.N_MODULE_MAX):
            for diode in range(quiet.NDIODE_MODULE):
                # Add the CES-diode to a list
                if acceptlist.isIncluded(scan,seg,module,diode):
                    masterlist.append((scan,seg,module,diode))
                else:
                    masterexlist.append((scan,seg,module,diode))
    # And turn the list into a set
    masterset=set(masterlist)
    masterexset=set(masterexlist)

    list1=[]
    exlist1=[]
    # Run through first list
    for (scan,seg) in a1.scans.keys():
        # Only consider scans for that patch
        if (scan,seg) in acceptlist.scans.keys():
            for module in range(quiet.N_MODULE_MAX):
                for diode in range(quiet.NDIODE_MODULE):
                    # Add the CES-diode to a list
                    if a1.isIncluded(scan,seg,module,diode):
                    # Could add here acceptlist.isIncluded()
                    # But then there are three sets to consider really..
                        list1.append((scan,seg,module,diode))
                    else:
                        exlist1.append((scan,seg,module,diode))
    # And turn the list into a set
    set1=set(list1)
    exset1=set(exlist1)

    list2=[]
    exlist2=[]
    # Run through first list
    for (scan,seg) in a2.scans.keys():
        # Only consider scans for that patch
        if (scan,seg) in acceptlist.scans.keys():
            for module in range(quiet.N_MODULE_MAX):
                for diode in range(quiet.NDIODE_MODULE):
                    # Add the CES-diode to a list
                    if a2.isIncluded(scan,seg,module,diode):
                        list2.append((scan,seg,module,diode))
                    else:
                        exlist2.append((scan,seg,module,diode))
    # And turn the list into a set
    set2=set(list2)
    exset2=set(exlist2)

    # Now make the comparisons
    # NB the sets contain INCLUDED diodes
    a1only  = set1 - set2  # a1 only
    a2only  = set2 - set1  # a2 only
    common  = set1 & set2  # Intersection
    either  = set1 | set2  # Union
    neither = masterset - either # Complement

    print len(a1only),len(a2only),len(common),len(either),\
          len(neither),len(masterset),len(set1),len(set2)

    exa1only  = exset1 - exset2  # a1 only
    exa2only  = exset2 - exset1  # a2 only
    excommon  = exset1 & exset2  # Intersection
    exeither  = exset1 | exset2  # Union
    exneither = masterexset - exeither # Complement

    print len(exa1only),len(exa2only),len(excommon),len(exeither),\
          len(exneither),len(masterexset),len(exset1),len(exset2)

    # Write the Columbia-only CES-diodes, for that patch, to file
    outf='colonly.txt'
    colonly=copy.deepcopy(a2)
    for (scan,seg) in sorted(a2.scans.keys()):
        for module in range(quiet.N_MODULE_MAX):
            for diode in range(quiet.NDIODE_MODULE):
                if (scan,seg,module,diode) in a2only:
                    colonly.includeSegmentDiode(scan,seg,module,diode)
                else:
                    colonly.excludeSegmentDiode(scan,seg,module,diode)
    #colonly.printSummary()
    colonly.toFile(outf)
    print 'a2 written to %s' % outf

    return common,neither,a1only,a2only,set1,set2,masterset

#-----------------------------------------------------------

def ApplyMinDiodeAccept(threshold=0.6,acceptlist=None,verbose=False,
                        patch='patch_6a',discount=6):

    """
    Cut any CESs with fewer than threshold diodes
    Return the new acceptlist
    Polarization modules only
    """

    import copy
    import quiet

    assert(float(threshold) >= 0.0 and float(threshold) <= 1.0)

    if acceptlist is None:
        rundict,acceptdict=GenerateCMBRunlists()
        acceptlist=acceptdict[patch]

    acceptlist.printSummary()

    hamlet=copy.deepcopy(acceptlist)

    total=((quiet.N_MODULE_MAX-2)*quiet.NDIODE_MODULE)-discount
    for (scan,seg) in sorted(acceptlist.scans.keys()):
        # For each CES
        ndiodes=0
        for module in range(quiet.N_MODULE_MAX-2):
            for diode in range(quiet.NDIODE_MODULE):
                if acceptlist.isIncluded(scan,seg,module,diode): ndiodes+=1
        if float(ndiodes)/float(total) < threshold:
            hamlet.removeScan(scan,seg)
            if verbose:
                print 'Removed %i.%i with %3.1f per cent' \
                      %(scan,seg,100.0*float(ndiodes)/float(total))

    hamlet.printSummary()

    return hamlet

#-----------------------------------------------------------

def ReadRegulatingFraction(f='/data4/quiet/calib_data/jonathaz/hk/regulating_fraction.txt'):

    """
    Read list of 'regulating fraction' and return it as dictionary keyed by CES
    """

    import scanf

    data=open(f,'r')

    linef='%i %i %f %f %f'

    print 'Reading regulating fraction from %s' % f

    regulating_fraction={}
    for line in data:
        if line.startswith('#'): continue
        scan,seg,mjd_start,mjd_stop,regulating_fraction[(scan,seg)]=\
            scanf.sscanf(line,linef)

    return regulating_fraction

#-----------------------------------------------------------

def ReadAPEXSummary(f='/data4/quiet/calib_data/jonathaz/hk/pwv_humidity.txt'):

    """
    Read APEX summary data and return it as dictionary keyed by CES
    """

    import scanf

    data=open(f,'r')

    linef='%i %i %f %f %f %f %f %f %f'

    print 'Reading APEX summary from %s' % f

    apex_summary={}
    for line in data:
        if line.startswith('#'): continue
        scan,seg,mjd_start,mjd_stop,pwv_avg,pwv_std,hum_avg,hum_std, \
            apex_frac_valid=scanf.sscanf(line,linef)
        apex_summary[(scan,seg)]=[pwv_avg,pwv_std,hum_avg,hum_std,apex_frac_valid]

    return apex_summary

#-----------------------------------------------------------

EXAMPLE_CES_STATS='/data4/quiet/jonathaz4/maps_nocuts_100317_2a/ces_statistics_set001.txt'

def ReadCESStats(f=EXAMPLE_CES_STATS,include_chisq_6hz=False):

    """
    Read a ces_stats file from tod2map into a dictionary keyed by CES mod.di
    """

    import scanf

    data=open(f,'r')

    linef=' %i %i %i %i %i %f %f %f %f %f %f %f %f %f %f %f'
    if include_chisq_6hz: linef += ' %f'

    print 'Reading ces_stats from %s (slow!)' % f

    ces_stats={}
    for line in data:
        if line.startswith(' #'): continue
        if include_chisq_6hz:
            scan,seg,mod,di,accept,ratio,sigma0,alpha,fknee,chisq0,\
                chisq1,chisq2,chisq3,chisq4,chisq_tod,tod_absmax,chisq_6hz=\
                scanf.sscanf(line,linef)
            ces_stats[(scan,seg,mod,di)]=\
                [accept,ratio,sigma0,alpha,fknee,chisq0,\
                 chisq1,chisq2,chisq3,chisq4,chisq_tod,tod_absmax,chisq_6hz]
        else:
            scan,seg,mod,di,accept,ratio,sigma0,alpha,fknee,chisq0,\
                chisq1,chisq2,chisq3,chisq4,chisq_tod,tod_absmax=\
                scanf.sscanf(line,linef)
            ces_stats[(scan,seg,mod,di)]=\
                [accept,ratio,sigma0,alpha,fknee,chisq0,\
                 chisq1,chisq2,chisq3,chisq4,chisq_tod,tod_absmax]

    return ces_stats

#-----------------------------------------------------------

def CutCESStats(patch=None,master=None,verbose=False,cuttype=None):

    """
    Read ces_statistics file and apply cuts to a master acceptlist
    """

    import os,copy,scanf

    assert(patch is not None), 'patch is required for %s!' % whoami()
    assert(master is not None), 'master is required for %s!' % whoami()

    t2mdir='/data4/quiet/jonathaz4/maps_nocuts'
    t2mf='%s_reports/ces_statistics_%s.txt' % (patch[-2:],patch[-2:])
    fftf=os.path.join(t2mdir,t2mf)
    print '** (applying cuts from %s)' % fftf

    linef=' %i %i %i %i %i %f %f %f %f %f %f %f %f %f %f %f'

    f=open(fftf,'r')

    hamlet=copy.deepcopy(master)

    kneethresh=50.0e-3 # Hz
    chi2thresh=4.0
    todthresh=7.0

    i=0
    flushfreq=10000
    for line in f:
        if line.startswith(' #'): continue
        scan,seg,module,diode,accept,ratio,sigma0,alpha,fknee,\
            chi0,chi1,chi2,chi3,chi4,chitod,sigtod=scanf.sscanf(line,linef)
        # Recreate existing cuts (chi2 > 4.0, tod > 7.0)
        if cuttype is None:
            condition = abs(chi1) >= chi2thresh or abs(chi2) >= chi2thresh or \
                        abs(chi3) >= chi2thresh or abs(sigtod) >= todthresh \
                        or fknee >= kneethresh
        elif cuttype == 'low':
            condition = abs(chi1) > chi2thresh or abs(chi2) > chi2thresh
        elif cuttype == 'mid':
            condition = abs(chi3) > chi2thresh
        #elif cuttype == 'chi1only' or cuttype == 'low':
        #    condition = chi1 > chi2thresh
        #elif cuttype == 'chi2only' or cuttype == 'mid':
        #    condition = chi2 > chi2thresh
        #elif cuttype == 'chi3only' or cuttype == 'high':
        #    condition = chi3 > chi2thresh
        #elif cuttype == 'todonly':
        #    condition = sigtod > todthresh
        if condition:
            if verbose:
                print 'Excluding %i %i %i %i' % (scan,seg,module,diode)
                print line
            hamlet.excludeSegmentDiode(scan,seg,module,diode)
            if accept != 0:
                print line
                print accept,chi1,chi2,chi3,sigtod,chi2thresh,todthresh
                print (sigtod==todthresh)
                return
        else:
            if accept == 0:
                print line
                print accept,chi1,chi2,chi3,sigtod,chi2thresh,todthresh
                print (sigtod==todthresh)
                return
        i+=1
        if (i % flushfreq == 0): print 'Reached %i' %i

    hamlet.title = 'columbia-%s-%s' % (patch[-2:],cuttype)

    hamlet.removeEmptyScans()

    return hamlet

#-----------------------------------------------------------

def PostProcessAccept(ain=None,aout=None,verbose=False):

    """
    This is a new, more general version of PostProcessAcceptlists
    """

    import quiet,copy

    mods=range(quiet.N_MODULE_MAX)
    diodes=range(quiet.NDIODE_MODULE)

    hamlet = copy.deepcopy(ain)

    sout=set(aout.scans.keys())
    for (scan,seg) in sorted(ain.scans.keys()):
        assert((scan,seg) in sout), '%i.%i - missing information' % (scan,seg)
        for module in mods:
            for diode in diodes:
                if ain.isExcluded(scan,seg,module,diode) or \
                       aout.isExcluded(scan,seg,module,diode):
                    hamlet.excludeSegmentDiode(scan,seg,module,diode)
                if verbose: print 'Cut %i.%i RQ% %i' % (scan,seg,module,diode)

    return hamlet

#-----------------------------------------------------------

def PostProcessAcceptlists(acceptin=None,acceptout=None,strict=True):

    """
    Compare input and output acceptlists for tod2map
    Assumes output list will have all the CESs from the input list
    IN FACT: assumes output CES list is a subset of the input list
    But that the output diodes will be a further subset to reject
    Print some summaries
    Output an acceptlist containing what was removed
    Also output a merged acceptlist
    """

    import os,copy
    import quiet

    dir='/data4/quiet/jonathaz4/maps_columbia_100117a/cmb/patch_6a'
    if acceptin is None: # Default
        acceptin='mlaccept-100117a.txt'
        print 'Using %s' % acceptin
        acceptin=os.path.join(dir,acceptin)
        input=quiet.AcceptedScanList(acceptin)
    else:
        input=acceptin
    if acceptout is None: # Default
        acceptout='accepted_6a_columbia_100117a.sort'
        print 'Using %s' % acceptout
        acceptout=os.path.join(dir,acceptout)
        output=quiet.AcceptedScanList(acceptout)
    else:
        output=acceptout

    #nout0=output.printSummary()

    output.printSummary()

    print len(input.scans.keys())
    print 'Input:'
    input.printSummary()
    print len(input.scans.keys())

    # If the following is not true, fall over
    # Caused by inappropriate NUM_SUBSETS_PER_TARGET in tod2map par file
    if strict:
        assert(len(output.scans.keys())==len(input.scans.keys())), \
           'The tod2map in and out lists have different numbers of CESs'

    # First run through the output list and exclude any diodes
    # that have already been excluded in the input list
    #si=set(input.scans.keys())
    for (scan,seg) in sorted(output.scans.keys()):
        #if (scan,seg) not in si: continue
        for module in range(quiet.N_MODULE_MAX):
            for diode in range(quiet.NDIODE_MODULE):
                if input.isExcluded(scan,seg,module,diode):
                    #if output.isIncluded(scan,seg,module,diode):
                    #    print scan,seg,module,diode
                    output.excludeSegmentDiode(scan,seg,module,diode)

    print 'Output (after chopping out already-cut diodes):'
    output.printSummary()

    # TEMPORARY? NEED TO SORT OUT BELOW HERE
    #for (scan,seg) in sorted(input.scans.keys()):
    #    if (scan,seg) not in output.scans.keys():
    #        for module in range(quiet.N_MODULE_MAX):
    #            for diode in range(quiet.NDIODE_MODULE):
    #                input.excludeSegmentDiode(scan,seg,module,diode)

    # Now run through removed and make it contain any diodes that
    # were in input but were not in output
    removed=copy.deepcopy(input)
    for (scan,seg) in sorted(removed.scans.keys()):
        for module in range(quiet.N_MODULE_MAX):
            for diode in range(quiet.NDIODE_MODULE):
                # WORKING HERE - NEED TO FIX THE O.KEYS() TEST
                #if (scan,seg) in output.scans.keys():
                if input.isIncluded(scan,seg,module,diode) and \
                       output.isExcluded(scan,seg,module,diode):
                    removed.includeSegmentDiode(scan,seg,module,diode)
                else:
                    removed.excludeSegmentDiode(scan,seg,module,diode)
                #else:
                #    removed.excludeSegmentDiode(scan,seg,module,diode)

    # These have to come after those two loops, because printSummary()
    # includes a call to removeEmptyScans()...
    print 'Input:'
    nin=input.printSummary()
    #nout0=output.printSummary()
    print 'Output:'
    nout1=output.printSummary()
    print 'Removed:'
    nremoved=removed.printSummary()

    if strict: assert(nremoved==nin-nout1), \
       'Discrepancy %i != %i - %i ' % (nremoved,nin,nout1)

    fracremoved=100.0*float(nin-nout1)/float(nin)
    fracleft=100.0*float(nout1)/float(nin)
    print 'Chi^2 cut summary:'
    print 'input: %i diodes (100 per cent)' % nin
    print 'output: %i diodes (%3.1f per cent)' % (nout1,fracleft)
    print 'removed: %i diodes (%3.1f per cent)' % (nremoved,fracremoved)

    return removed,output

#-----------------------------------------------------------

def CutMispointedScans(acceptlist=None,runlist=None,verbose=False,\
                       patch=None):

    """
    A synonym wrapper for CutMislabelledScans
    Ignore (unnecessary) runlist argument
    """

    hamlet=CutMislabelledScans(acceptlist=acceptlist,\
                               verbose=verbose,patch=patch)

    return hamlet

#-----------------------------------------------------------

def CutMislabelledScans(acceptlist=None,verbose=False,patch=None):

    """
    Cut specific scans that are mislabelled/pointed incorrectly
    List is complete for patch 6a, but
    incomplete for 2a,4a
    Absent for 7b and other patches
    """

    import copy

    # http://code.activestate.com/recipes/66062
    assert(patch is not None), 'patch is required for %s!' % whoami()

    # These are under the 'new' numbering scheme
    scans={}
    scans['patch_2a'] = [(377,6),(1302,1)] # Incomplete?
    scans['patch_4a'] = [(118,1),(1538,1),(1648,1)] # Incomplete?
    scans['patch_6a'] = [(783,1),(829,1),(829,2),(897,1),(1335,1),(1468,1),(1468,2),(1786,4)]
    # 1468.1 .2 should be 7b not 6a -> cut it
    # 783.1 - which patch? } -> cut them
    # 897.1 - which patch? }
    #scans['patch_7b'] = [] # Has none

    if patch not in ['patch_6a','patch_4a','patch_2a']:
        print '***No information for %s' % patch
        print '---> No changes to the acceptlist'
        return acceptlist

    if patch != 'patch_6a':
        print '***Mispointed scan list may be incomplete for %s' % patch

    # Set up default runlist and acceptlist if not user-supplied
    if acceptlist is None:
        rundict,acceptdict=GenerateCMBRunlists()
        acceptlist=acceptdict[patch]

    hamlet   = copy.deepcopy(acceptlist) # Fresh copy

    sa=set(acceptlist.scans.keys())
    for (scan,seg) in scans[patch]:
        if (scan,seg) in sa:
            #if not any(acceptlist.scans[scan,seg]): continue # Useful?
            hamlet.removeScan(scan,seg)
            if verbose:
                print '%i.%i cut from list' % (scan,seg)
        else:
            if verbose:
                print '%i.%i not in list ' % (scan,seg)

    hamlet.removeEmptyScans()

    return hamlet

#-----------------------------------------------------------

def CutPhaseSwitchHumidity(acceptlist=None,verbose=False,patch=None):

    """
    Cut specific scans for which the phase switches are unstable
    There is a cutlist for each of RQ11 and RQ12
    The two files contain information for all patches
    """

    import os,copy

    assert(patch is not None), 'patch is required for %s!' % whoami()

    cutsdir='/data4/quiet/nulltests/cuts/'

    modules=[11,12]
    files={11:'rq11_ps21_badscans.txt',12:'rq12_ps21_badscans.txt'}

    scans={}
    for module in sorted(files.keys()):
        fi=os.path.join(cutsdir,files[module])
        f=open(fi,'r')
        for line in f:
            if line.startswith('#') or line == '': continue
            scan,seg=line.split()
            if module not in scans.keys():
                scans[module]=[]
            scans[module].append((int(scan),int(seg)))

    # Set up default runlist and acceptlist if not user-supplied
    if acceptlist is None:
        rundict,acceptdict=GenerateCMBRunlists()
        acceptlist=acceptdict[patch]

    hamlet = copy.deepcopy(acceptlist) # Fresh copy

    sa=set(acceptlist.scans.keys())
    for module in modules:
        for (scan,seg) in scans[module]:
            if (scan,seg) in sa:
                hamlet.excludeSegmentModules(scan,seg,module)
                if verbose:
                    print '%i.%i RQ%2i cut from list' % (scan,seg,module)
            else:
                if verbose:
                    print '%i.%i RQ%2i not in list ' % (scan,seg,module)

    hamlet.removeEmptyScans()

    return hamlet

#-----------------------------------------------------------

def CutScanLength(acceptlist=None,tmin=None,runlist=None,\
                  verbose=False,patch=None):

    """
    For a given patch, cut from an acceptlist those CESs
    that are shorter than the specified duration in minutes
    """

    import quiet,copy
    import ephem

    assert(patch is not None), 'patch is required for %s!' % whoami()
    assert(tmin is not None), 'tmin is required for %s!' % whoami()

    # Set up default runlist and acceptlist if not user-supplied
    if runlist is None:
        rundict,acceptdict=GenerateCMBRunlists()
        runlist=rundict[patch]
        if acceptlist is None:
            acceptlist=acceptdict[patch]

    mjdstart={}
    mjdstop={}
    scanlength={}

    hamlet = copy.deepcopy(acceptlist) # Fresh copy
    for scan,seg in sorted(acceptlist.scans.keys()):
        az,el,dk,mjdstart[(scan,seg)],mjdstop[(scan,seg)]=\
             FetchPointing(scan,seg,runlist,patch)
        #print '%i.%i el %2.1f' % (scan,seg,el[(scan,seg)])
        scanlength[(scan,seg)] = (mjdstop[(scan,seg)] - mjdstart[(scan,seg)])\
            / ephem.minute # Put scanlength into minutes
        if scanlength[(scan,seg)] <= tmin: # cut scans shorter than this
            hamlet.removeScan(scan,seg)
            if verbose:
                print '%i.%i cut  %2.1f' % (scan,seg,scanlength[(scan,seg)])
        else:
            if verbose:
                print '%i.%i keep %2.1f' % (scan,seg,scanlength[(scan,seg)])

    hamlet.removeEmptyScans()

    return hamlet

#-----------------------------------------------------------

def WriteScanTimes(acceptmask=None,acceptlist=None,runlist=None,\
                  verbose=False,patch=None,module=None,diode=None):

    """
    Write to file CES scan id mjdstart mjdstop accept
    given a runlist, acceptlist and an acceptlist mask to apply
    Do this optionally for a specific diode only
    """

    import quiet

    assert(patch is not None), 'patch required!'

    # Set up default runlist and acceptlist if not user-supplied
    if runlist is None:
        rundict,acceptdict=GenerateCMBRunlists()
        runlist=rundict[patch]
        if acceptlist is None:
            acceptlist=acceptdict[patch]
        if acceptmask is None:
            acceptmask=acceptlist

    modules=[]; diodes=[];
    if module is None:
        modules = range(quiet.N_MODULES_MAX)
    else: modules=[module]
    if diode is None:
        diodes = range(quiet.NDIODE_MODULE)
    else: diodes=[diode]

    outfile='season-times.txt'
    outf=open(outfile,'w')

    linef='%i %i %i %i %10.5f %10.5f %i %s\n'
    sa=sorted(acceptlist.scans.keys())
    sm=set(acceptmask.scans.keys())
    for scan,seg in sa:
        az,el,dk,mjdstart,mjdstop=FetchPointing(scan,seg,runlist,patch)
        for mod in modules:
            for di in diodes:
                if (scan,seg) in sm:
                    if acceptmask.isIncluded(scan,seg,mod,di):
                        accept=1
                    else: accept=0
                else: accept=0
                line=linef % (scan,seg,mod,di,mjdstart,mjdstop,accept,patch)
                outf.write(line)

    outf.close()

    return None

#-----------------------------------------------------------

def WriteScanLength(acceptlist=None,runlist=None,\
                  verbose=False,patch=None):

    """
    For an acceptlist, fetch the scan lengths and write them to a file
    """

    import ephem
    import quiet

    # Set up default runlist and acceptlist if not user-supplied
    if runlist is None:
        rundict,acceptdict=GenerateCMBRunlists()
        runlist=rundict[patch]
        if acceptlist is None:
            acceptlist=acceptdict[patch]

    mjdstart={}
    mjdstop={}
    scanlength={}

    outfile='scanlengths.txt'
    outf=open(outfile,'w')

    for scan,seg in sorted(acceptlist.scans.keys()):
        az,el,dk,mjdstart[(scan,seg)],mjdstop[(scan,seg)]=\
             FetchPointing(scan,seg,runlist,patch)
        scanlength[(scan,seg)] = (mjdstop[(scan,seg)] - mjdstart[(scan,seg)])\
            / ephem.second # Put scanlength into seconds
        outf.write('%i %i %5.1f\n' % (scan,seg,scanlength[(scan,seg)]))

    outf.close()

    return None

#-----------------------------------------------------------


def CutSunMoonDistance(acceptlist=None,sunthresh=None,moonthresh=None,\
                       runlist=None,verbose=False,mods=None,patch='patch_6a'):

    """
    Cut from an acceptlist those CES-modules that are within
    some distance of the Sun or Moon at any time during the scan
    Remove all modules by default (can specify)
    """

    import quiet,copy,ephem

    # If no sunthresh supplied
    if sunthresh is None: sunthresh=-1.0 # A hard short-cut to subvert..
    # If no moonthresh supplied
    if moonthresh is None: moonthresh=-1.0 # A hard short-cut to subvert..

    # Set up default runlist and acceptlist if not user-supplied
    if runlist is None:
        rundict,acceptdict=GenerateCMBRunlists()
        runlist=rundict[patch]
        if acceptlist is None:
            acceptlist=acceptdict[patch]

    mjdstart={}
    mjdstop={}
    scanlength={}
    #print elmin,elmax

    hamlet   = copy.deepcopy(acceptlist) # Fresh copy

    objects=FetchSkyAngles(acceptlist=acceptlist,\
                           runlist=runlist,verbose=False)

    for (scan,seg,nmjd) in sorted(objects.keys()):
        sunsep=ang2deg(objects[(scan,seg,nmjd)]['Sunsep'])
        moonsep=ang2deg(objects[(scan,seg,nmjd)]['Moonsep'])
        if sunsep <= sunthresh or moonsep <= moonthresh:
            # To avoid double removals
            if (scan,seg) in hamlet.scans.keys():
                if mods is None: # Default remove all modules
                    hamlet.removeScan(scan,seg)
                else: # Remove only specified modules
                    hamlet.excludeSegmentModules(scan,seg,mods)
            if verbose:
                print '%i.%i %i cut  %3.1f %3.1f' \
                      % (scan,seg,nmjd,sunsep,moonsep)
        else:
            if verbose:
                print '%i.%i %i keep %3.1f %3.1f' \
                      % (scan,seg,nmjd,sunsep,moonsep)

    #if verbose: print hamlet

    hamlet.removeEmptyScans()

    return hamlet

#-----------------------------------------------------------

def CutElRange(acceptlist=None,mods=None,elmin=60.0,elmax=65.0,\
               runlist=None,verbose=False,patch=None):

    """
    Cut from an acceptlist those CES-modules that are within
    the specified el range
    Remove all modules by default (can specify)
    """

    import quiet,copy

    # Set up default runlist and acceptlist if not user-supplied
    if runlist is None:
        rundict,acceptdict=GenerateCMBRunlists()
        if patch is None: patch='patch_6a'
        runlist=rundict[patch]
        if acceptlist is None:
            acceptlist=acceptdict[patch]

    az={}
    el={}
    dk={}
    #print elmin,elmax

    hamlet   = copy.deepcopy(acceptlist) # Fresh copy
    for scan,seg in sorted(acceptlist.scans.keys()):
        az[(scan,seg)],el[(scan,seg)],dk[(scan,seg)],\
            mjdstart,mjdstop=FetchPointing(scan,seg,runlist,patch)
        #print '%i.%i el %2.1f' % (scan,seg,el[(scan,seg)])
        if el[(scan,seg)] >= elmin and el[(scan,seg)] <= elmax:
            if verbose: print '%i.%i cut  %2.1f' % (scan,seg,el[(scan,seg)])
            if mods is None: # Default remove all modules
                hamlet.removeScan(scan,seg)
            else: # Remove only specified modules
                hamlet.excludeSegmentModules(scan,seg,mods)
        else:
            if verbose: print '%i.%i keep %2.1f' % (scan,seg,el[(scan,seg)])

    if verbose: print hamlet

    hamlet.removeEmptyScans()

    return hamlet

#-----------------------------------------------------------

def ces_translate_old(scan,seg,direction=None):

    """
    Wrapper for Sigurd's ces_translate
    Would be better to wrap the c++ directly rather than via this route
    """

    import os

    exedir='/usit/titan/u1/jonathaz/quiet_svn/oslo/src/cpp/utils'
    exe='ces_translate'
    executable=os.path.join(exedir,exe)

    if direction=='-r':
        command = '%s %i %i' % (executable,scan,seg)
    else:
        command = '%s %s %i %i' % (executable,direction,scan,seg)

    newscan,newseg=os.popen(command).readlines()[0].split()

    return (int(newscan),int(newseg))

#-----------------------------------------------------------

def CutTTFailedMapScans(acceptlist=None,verbose=False,patch=None):

    """
    Cut from an acceptlist those CESs that
    feature in the following list of scans
    """

    import quiet
    import copy

    if patch != 'patch_6a':
        print '***Scans lists are only for patch_6a so far!!!'
        if acceptlist is not None:
            print '---> No changes to the acceptlist'
            return acceptlist
        else:
            return None

    # Set up default acceptlist if not user-supplied
    if acceptlist is None:
        rundict,acceptdict=GenerateCMBRunlists()
        if patch is None: patch='patch_6a'
        acceptlist=acceptdict[patch]

    TTlist=[(394,1),(1060,1),\
            (1071,1),(1082,2),\
            (1135,2),(1135,4),\
            (1141,2),(1151,2),\
            (1159,2),(1170,2),\
            (1180,2),(1281,1)]

    hamlet   = copy.deepcopy(acceptlist) # Fresh copy

    sh=set(hamlet.scans.keys())

    for (scan,seg) in TTlist:
        if (scan,seg) in sh:
            hamlet.removeScan(scan,seg)
            if verbose:
                print '%i.%i cut' % (scan,seg)

    hamlet.removeEmptyScans()

    return hamlet

#-----------------------------------------------------------

def CutFailedWeatherDiodes(acceptlist=None,verbose=False,\
                           patch=None,cuttype=None):

    """
    Cut from an acceptlist those CESs that
    fail the weather cut(s)
    There are now several versions of the weather cut (cuttype):
    1  30,60,300 1-sigma
    2  whole CESs (weather_cut_100114.txt)
    3  30,60 1.5-sigma
    4  whole CESs (TTrejectlist_100114.txt)
    5  30,60 1-sigma
    6  30,60 2-sigma (asymmetric)
    7  10,30 2-sigma (asymmetric)
    8  10,30 5-sigma (asymm with Gaussian)
    9  10,30 5-sigma (symm Gaussian)
    10 10,30 5-sigma (symm Gaussian) but with half-heaviside
    11 As 10, but return acceptlist of the diodes removed
    12 TP_3sigma_70percent_cutthese_2a
    13 As 10, but remove individual diodes not whole scans
    """

    import quiet
    import os,copy

    assert(patch is not None), 'patch is required for %s!' % whoami()
    assert(cuttype is not None), 'cuttype is required for %s!' % whoami()

    # Set up default acceptlist if not user-supplied
    if acceptlist is None:
        rundict,acceptdict=GenerateCMBRunlists()
        acceptlist=acceptdict[patch]

    # Specify the cuts directory
    cutsdir='/data4/quiet/nulltests/cuts/'

    hamlet = copy.deepcopy(acceptlist) # Fresh copy

    hamlet.title='columbia-%s-%i' % (patch,cuttype)

    if cuttype==1:
        # These are the patch_6a acceptlists for Laura's (3) weather cuts
        keepfs=['weather_acceptlist_asrun_300.txt',\
                'weather_acceptlist_asrun_30.txt',\
                'weather_acceptlist_asrun_60.txt']
        # Read Laura's three weather acceptlists
        # NB These files contain lists of diodes to KEEP
        # A scan is not in these files if all its diodes are rejected
        acceptd={}
        for weather_keepf_6a in keepfs:
            acceptd[weather_keepf_6a]=\
                quiet.AcceptedScanList(os.path.join(cutsdir,weather_keepf_6a))

        a1=acceptd[keepfs[0]]
        a2=acceptd[keepfs[1]]
        a3=acceptd[keepfs[2]]

        for (scan,seg) in sorted(acceptlist.scans.keys()):
            for module in range(quiet.N_MODULE_MAX):
                for diode in range(quiet.NDIODE_MODULE):
                    # NB boolean test is reversed because of the format of ai
                    if (scan,seg) in a1.scans.keys():
                        keep1=not(a1.isIncluded(scan,seg,module,diode))
                    else:
                        keep1=False # Reject if not in list
                    if (scan,seg) in a2.scans.keys():
                        keep2=not(a2.isIncluded(scan,seg,module,diode))
                    else:
                        keep2=False # Reject if not in list
                    if (scan,seg) in a3.scans.keys():
                        keep3=not(a3.isIncluded(scan,seg,module,diode))
                    else:
                        keep3=False # Reject if not in list

                    # Check whether the diode was already kept
                    alreadykept=hamlet.isIncluded(scan,seg,module,diode)

                    if (keep1 and keep2 and keep3 and alreadykept):
                        hamlet.includeSegmentDiode(scan,seg,module,diode)
                    else:
                        hamlet.excludeSegmentDiode(scan,seg,module,diode)

    elif cuttype==2:
        # This is an acceptlist of CESs to reject
        weatherf='weather_cut_100114.txt'
        f=os.path.join(cutsdir,weatherf)
        rejectlist=quiet.AcceptedScanList(f)

        # Run through the acceptlist, and
        # throw the scan out of hamlet if it's also in the rejectlist
        for (scan,seg) in sorted(acceptlist.scans.keys()):
            if (scan,seg) in rejectlist.scans.keys():
                hamlet.removeScan(scan,seg)
                if verbose:
                    print 'Removed CES %i.%i' % (scan,seg)

    elif cuttype==3:
        # These are the patch_6a acceptlists for Laura's (2) weather cuts
        # 1.5-sigma weather cuts
        keepfs=['weather_acceptlist_asrun_30_1.5sigma.txt',\
                'weather_acceptlist_asrun_60_1.5sigma.txt']
        # Read Laura's two 1.5-sigma weather acceptlists
        # NB These files contain lists of diodes to KEEP
        # A scan is not in these files if all its diodes are rejected
        acceptd={}
        for weather_keepf_6a in keepfs:
            acceptd[weather_keepf_6a]=\
                quiet.AcceptedScanList(os.path.join(cutsdir,weather_keepf_6a))

        a1=acceptd[keepfs[0]]
        a2=acceptd[keepfs[1]]

        for (scan,seg) in sorted(acceptlist.scans.keys()):
            for module in range(quiet.N_MODULE_MAX):
                for diode in range(quiet.NDIODE_MODULE):
                    # NB boolean test is reversed because of the format of ai
                    if (scan,seg) in a1.scans.keys():
                        keep1=not(a1.isIncluded(scan,seg,module,diode))
                    else:
                        keep1=False # Reject if not in list
                    if (scan,seg) in a2.scans.keys():
                        keep2=not(a2.isIncluded(scan,seg,module,diode))
                    else:
                        keep2=False # Reject if not in list

                    # Check whether the diode was already kept
                    alreadykept=hamlet.isIncluded(scan,seg,module,diode)

                    if (keep1 and keep2 and alreadykept):
                        hamlet.includeSegmentDiode(scan,seg,module,diode)
                    else:
                        hamlet.excludeSegmentDiode(scan,seg,module,diode)

    elif cuttype==4:
        # This is a list of CESs to reject
        weatherf='TTrejectlist_100114.txt'
        rejectlist=open(os.path.join(cutsdir,weatherf),'r')
        # Run through the list, and
        # throw the scan out of hamlet if it's there
        for line in rejectlist:
            strscan,strseg = line.split()
            scan=int(strscan)
            seg=int(strseg)
            if (scan,seg) in sorted(acceptlist.scans.keys()):
                hamlet.removeScan(scan,seg)
                if verbose: print 'Removed %i.%i' % (scan,seg)
        rejectlist.close()

    elif cuttype==5:
        # These are the patch_6a acceptlists for Laura's (2) weather cuts
        keepfs=['weather_acceptlist_asrun_30.txt',\
                'weather_acceptlist_asrun_60.txt']
        # Read Laura's three weather acceptlists
        # NB These files contain lists of diodes to KEEP
        # A scan is not in these files if all its diodes are rejected
        acceptd={}
        for weather_keepf_6a in keepfs:
            acceptd[weather_keepf_6a]=\
                quiet.AcceptedScanList(os.path.join(cutsdir,weather_keepf_6a))

        a1=acceptd[keepfs[0]]
        a2=acceptd[keepfs[1]]

        for (scan,seg) in sorted(acceptlist.scans.keys()):
            for module in range(quiet.N_MODULE_MAX):
                for diode in range(quiet.NDIODE_MODULE):
                    # NB boolean test is reversed because of the format of ai
                    if (scan,seg) in a1.scans.keys():
                        keep1=not(a1.isIncluded(scan,seg,module,diode))
                    else:
                        keep1=False # Reject if not in list
                    if (scan,seg) in a2.scans.keys():
                        keep2=not(a2.isIncluded(scan,seg,module,diode))
                    else:
                        keep2=False # Reject if not in list

                    # Check whether the diode was already kept
                    alreadykept=hamlet.isIncluded(scan,seg,module,diode)

                    if (keep1 and keep2 and alreadykept):
                        hamlet.includeSegmentDiode(scan,seg,module,diode)
                    else:
                        hamlet.excludeSegmentDiode(scan,seg,module,diode)

    elif cuttype==6:
        # These are the patch_6a acceptlists for Laura's (2) weather cuts
        keepfs=['weather_acceptlist_asrun_30_2sigmaasym.txt',\
                'weather_acceptlist_asrun_60_2sigmaasym.txt']
        # Read Laura's two weather acceptlists
        # NB These files contain lists of diodes to KEEP
        # A scan is not in these files if all its diodes are rejected
        acceptd={}
        for weather_keepf_6a in keepfs:
            acceptd[weather_keepf_6a]=\
                quiet.AcceptedScanList(os.path.join(cutsdir,weather_keepf_6a))

        a1=acceptd[keepfs[0]]
        a2=acceptd[keepfs[1]]

        for (scan,seg) in sorted(acceptlist.scans.keys()):
            for module in range(quiet.N_MODULE_MAX):
                for diode in range(quiet.NDIODE_MODULE):
                    # NB boolean test is reversed because of the format of ai
                    if (scan,seg) in a1.scans.keys():
                        keep1=not(a1.isIncluded(scan,seg,module,diode))
                    else:
                        keep1=False # Reject if not in list
                    if (scan,seg) in a2.scans.keys():
                        keep2=not(a2.isIncluded(scan,seg,module,diode))
                    else:
                        keep2=False # Reject if not in list

                    # Check whether the diode was already kept
                    alreadykept=hamlet.isIncluded(scan,seg,module,diode)

                    if (keep1 and keep2 and alreadykept):
                        hamlet.includeSegmentDiode(scan,seg,module,diode)
                    else:
                        hamlet.excludeSegmentDiode(scan,seg,module,diode)

    elif cuttype==7:
        # These are the patch_6a acceptlists for Laura's (2) weather cuts
        keepfs=['weather_acceptlist_asrun_10_2sigmaasym.txt',\
                'weather_acceptlist_asrun_30_2sigmaasym.txt']
        # Read Laura's two weather acceptlists
        # NB These files contain lists of diodes to KEEP
        # A scan is not in these files if all its diodes are rejected
        acceptd={}
        for weather_keepf_6a in keepfs:
            acceptd[weather_keepf_6a]=\
                quiet.AcceptedScanList(os.path.join(cutsdir,weather_keepf_6a))

        a1=acceptd[keepfs[0]]
        a2=acceptd[keepfs[1]]

        for (scan,seg) in sorted(acceptlist.scans.keys()):
            for module in range(quiet.N_MODULE_MAX):
                for diode in range(quiet.NDIODE_MODULE):
                    # NB boolean test is reversed because of the format of ai
                    if (scan,seg) in a1.scans.keys():
                        keep1=not(a1.isIncluded(scan,seg,module,diode))
                    else:
                        keep1=False # Reject if not in list
                    if (scan,seg) in a2.scans.keys():
                        keep2=not(a2.isIncluded(scan,seg,module,diode))
                    else:
                        keep2=False # Reject if not in list

                    # Check whether the diode was already kept
                    alreadykept=hamlet.isIncluded(scan,seg,module,diode)

                    if (keep1 and keep2 and alreadykept):
                        hamlet.includeSegmentDiode(scan,seg,module,diode)
                    else:
                        hamlet.excludeSegmentDiode(scan,seg,module,diode)

    elif cuttype==8:
        # These are the patch_6a acceptlists for Laura's (2) weather cuts
        keepfs=['weather_acceptlist_asrun_10_5sigma.txt',\
                'weather_acceptlist_asrun_30_5sigma.txt']
        # Read Laura's two weather acceptlists
        # NB These files contain lists of diodes to KEEP
        # A scan is not in these files if all its diodes are rejected
        acceptd={}
        for weather_keepf_6a in keepfs:
            acceptd[weather_keepf_6a]=\
                quiet.AcceptedScanList(os.path.join(cutsdir,weather_keepf_6a))

        a1=acceptd[keepfs[0]]
        a2=acceptd[keepfs[1]]

        for (scan,seg) in sorted(acceptlist.scans.keys()):
            for module in range(quiet.N_MODULE_MAX):
                for diode in range(quiet.NDIODE_MODULE):
                    # NB boolean test is reversed because of the format of ai
                    if (scan,seg) in a1.scans.keys():
                        keep1=not(a1.isIncluded(scan,seg,module,diode))
                    else:
                        keep1=False # Reject if not in list
                    if (scan,seg) in a2.scans.keys():
                        keep2=not(a2.isIncluded(scan,seg,module,diode))
                    else:
                        keep2=False # Reject if not in list

                    # Check whether the diode was already kept
                    alreadykept=hamlet.isIncluded(scan,seg,module,diode)

                    if (keep1 and keep2 and alreadykept):
                        hamlet.includeSegmentDiode(scan,seg,module,diode)
                    else:
                        hamlet.excludeSegmentDiode(scan,seg,module,diode)

    elif cuttype==9:
        # These are the patch_6a acceptlists for Laura's (2) weather cuts
        keepfs=['weather_acceptlist_asrun_10_5sigsym.txt',\
                'weather_acceptlist_asrun_30_5sigsym.txt']
        # Read Laura's two weather acceptlists
        # NB These files contain lists of diodes to KEEP
        # A scan is not in these files if all its diodes are rejected
        acceptd={}
        for weather_keepf_6a in keepfs:
            acceptd[weather_keepf_6a]=\
                quiet.AcceptedScanList(os.path.join(cutsdir,weather_keepf_6a))

        a1=acceptd[keepfs[0]]
        a2=acceptd[keepfs[1]]

        for (scan,seg) in sorted(acceptlist.scans.keys()):
            for module in range(quiet.N_MODULE_MAX):
                for diode in range(quiet.NDIODE_MODULE):
                    # NB boolean test is reversed because of the format of ai
                    if (scan,seg) in a1.scans.keys():
                        keep1=not(a1.isIncluded(scan,seg,module,diode))
                    else:
                        keep1=False # Reject if not in list
                    if (scan,seg) in a2.scans.keys():
                        keep2=not(a2.isIncluded(scan,seg,module,diode))
                    else:
                        keep2=False # Reject if not in list

                    # Check whether the diode was already kept
                    alreadykept=hamlet.isIncluded(scan,seg,module,diode)

                    if (keep1 and keep2 and alreadykept):
                        hamlet.includeSegmentDiode(scan,seg,module,diode)
                    else:
                        hamlet.excludeSegmentDiode(scan,seg,module,diode)

    elif cuttype==10:
        # This cut applies a heaviside to the below
        # These are the patch acceptlists for Laura's (2) weather cuts
        if patch=='patch_6a':
            keepfs=['weather_acceptlist_asrun_10_5sigsym.txt',\
                    'weather_acceptlist_asrun_30_5sigsym.txt']
        elif patch=='patch_2a':
            keepfs=['weather_acceptlist_2a_asrun_10_5sigsym.txt',\
                    'weather_acceptlist_2a_asrun_30_5sigsym.txt']
#        if patch=='patch_6a':
#            keepfs=['weather_acceptlist_asrun_10_patch_6a_100719.txt',\
#                    'weather_acceptlist_asrun_30_patch_6a_100719.txt']
#        elif patch=='patch_2a':
#            keepfs=['weather_acceptlist_asrun_30_patch_2a_100719.txt',\
#                    'weather_acceptlist_asrun_30_patch_2a_100719.txt']
            # Hacks for testing 10 v 30 s
            #keepfs=['weather_acceptlist_2a_asrun_10_5sigsym.txt',\
            #        'weather_acceptlist_2a_asrun_10_5sigsym.txt']
            #keepfs=['weather_acceptlist_2a_asrun_30_5sigsym.txt',\
            #        'weather_acceptlist_2a_asrun_30_5sigsym.txt']
            # Hack for testing all patches v 2a histograms
            #keepfs=['weather_acceptlist_asrun_10_allpatches.txt',\
            #        'weather_acceptlist_asrun_30_allpatches.txt']
        elif patch=='patch_4a':
            keepfs=['weather_4a_acceptlist_asrun_10_5sigma.txt',\
                    'weather_4a_acceptlist_asrun_30_5sigma.txt']
        elif patch=='patch_7b':
            keepfs=['weather_7b_acceptlist_asrun_10_5sigma.txt',\
                    'weather_7b_acceptlist_asrun_30_5sigma.txt']
#        elif patch=='patch_4a':
#            keepfs=['weather_acceptlist_asrun_10_patch_4a_100719.txt',\
#                    'weather_acceptlist_asrun_30_patch_4a_100719.txt']
#        elif patch=='patch_7b':
#            keepfs=['weather_acceptlist_asrun_10_patch_7b_100719.txt',\
#                    'weather_acceptlist_asrun_30_patch_7b_100719.txt']
        # Read Laura's two weather acceptlists
        # NB These files contain lists of diodes to KEEP
        # A scan is not in these files if all its diodes are rejected
        acceptd={}
        for weather_keepf in keepfs:
            acceptd[weather_keepf]=\
                quiet.AcceptedScanList(os.path.join(cutsdir,weather_keepf))

        a1=acceptd[keepfs[0]]
        a2=acceptd[keepfs[1]]

        threshold=0.7
        discount=6 # Modules and diodes already statically removed
        nbefore=len(hamlet.scans.keys())
        total=((quiet.N_MODULE_MAX-2)*quiet.NDIODE_MODULE)-discount
        for (scan,seg) in sorted(acceptlist.scans.keys()):
            ndiodes=0
            for module in range(quiet.N_MODULE_MAX-2):
                for diode in range(quiet.NDIODE_MODULE):
                    # NB boolean test is reversed because of the format of ai
                    if (scan,seg) in a1.scans.keys():
                        keep1=not(a1.isIncluded(scan,seg,module,diode))
                    else:
                        keep1=False # Reject if not in list
                    if (scan,seg) in a2.scans.keys():
                        keep2=not(a2.isIncluded(scan,seg,module,diode))
                    else:
                        keep2=False # Reject if not in list
                    # Check whether the diode was already kept
                    alreadykept=hamlet.isIncluded(scan,seg,module,diode)
                    # Keep track of how many diodes remain for this CES
                    if (keep1 and keep2 and alreadykept):
                        ndiodes+=1

            frac=float(ndiodes)/float(total)
            assert (0.0 <= frac <= 1.0)
            if frac < threshold:
                hamlet.removeScan(scan,seg)
                if verbose:
                    print 'Removed %i.%i with %3.1f per cent' \
                          %(scan,seg,100.0*float(ndiodes)/float(total))

        nafter=len(hamlet.scans.keys())
        print 'For cuttype %i, removing CESs with < %3.1f per cent of diodes' \
              % (cuttype,100.0*threshold)
        print '  and that removed %i scans' % (nbefore-nafter)

    elif cuttype==11:
        # SPECIAL: As 10, but return acceptlist of what has been cut
        # This cut applies a heaviside to the below
        # These are the patch acceptlists for Laura's (2) weather cuts
        if patch=='patch_6a':
            keepfs=['weather_acceptlist_asrun_10_5sigsym.txt',\
                    'weather_acceptlist_asrun_30_5sigsym.txt']
        elif patch=='patch_2a':
            keepfs=['weather_acceptlist_2a_asrun_10_5sigsym.txt',\
                    'weather_acceptlist_2a_asrun_30_5sigsym.txt']
        elif patch=='patch_4a':
            keepfs=['weather_4a_acceptlist_asrun_10_5sigma.txt',\
                    'weather_4a_acceptlist_asrun_30_5sigma.txt']
        elif patch=='patch_7b':
            keepfs=['weather_7b_acceptlist_asrun_10_5sigma.txt',\
                    'weather_7b_acceptlist_asrun_30_5sigma.txt']
        # Read Laura's two weather acceptlists
        # NB These files contain lists of diodes to KEEP
        # A scan is not in these files if all its diodes are rejected
        acceptd={}
        for weather_keepf in keepfs:
            acceptd[weather_keepf]=\
                quiet.AcceptedScanList(os.path.join(cutsdir,weather_keepf))

        a1=acceptd[keepfs[0]]
        a2=acceptd[keepfs[1]]

        threshold=0.7
        discount=6 # Modules and diodes already statically removed
        nbefore=len(hamlet.scans.keys())
        total=((quiet.N_MODULE_MAX-2)*quiet.NDIODE_MODULE)-discount
        for (scan,seg) in sorted(acceptlist.scans.keys()):
            ndiodes=0
            for module in range(quiet.N_MODULE_MAX-2): # Exclude TTs
                for diode in range(quiet.NDIODE_MODULE):
                    # NB boolean test is reversed because of the format of ai
                    if (scan,seg) in a1.scans.keys():
                        keep1=not(a1.isIncluded(scan,seg,module,diode))
                    else:
                        keep1=False # Reject if not in list
                    if (scan,seg) in a2.scans.keys():
                        keep2=not(a2.isIncluded(scan,seg,module,diode))
                    else:
                        keep2=False # Reject if not in list
                    # Check whether the diode was already kept
                    alreadykept=hamlet.isIncluded(scan,seg,module,diode)
                    # Keep track of how many diodes remain for this CES
                    if (keep1 and keep2 and alreadykept):
                        ndiodes+=1
                        hamlet.includeSegmentDiode(scan,seg,module,diode)
                    else:
                        hamlet.excludeSegmentDiode(scan,seg,module,diode)

            #ndiodes-=6 # Since the other Pol cuts have not been applied
            frac=float(ndiodes)/float(total)
            assert(0.0 <= frac <= 1.0)
            # Now retain only those with LESS THAN threshold per cent remaining
            if frac >= threshold:
                hamlet.removeScan(scan,seg)
                if verbose:
                    print 'Removed %i.%i with %3.1f per cent left in' \
                          %(scan,seg,100.0*frac)

        nafter=len(hamlet.scans.keys())
        print 'For cuttype %i, removing CESs with >= %3.1f per cent of diodes' \
              % (cuttype,100.0*threshold)
        print '  and that removed %i scans' % (nbefore-nafter)

    elif cuttype==12:
        if patch != 'patch_2a':
            print 'patch_2a only!'
            return None
        reject='TP_3sigma_70percent_cutthese_2a.txt'
        rejectlist=open(os.path.join(cutsdir,reject),'r')
        sa=set(acceptlist.scans.keys())
        for line in rejectlist:
            strscan,strseg=line.split()
            scan=int(strscan)
            seg=int(strseg)
            if (scan,seg) in sa:
                hamlet.removeScan(scan,seg)
                if verbose:
                    print 'Cut %i.%i' % (scan,seg)
        rejectlist.close()

    elif cuttype==13:
        # As 10, but remove individual diodes not whole scans
        # This cut applies a heaviside to the below
        # These are the patch acceptlists for Laura's (2) weather cuts
        if patch=='patch_6a':
            keepfs=['weather_acceptlist_asrun_10_5sigsym.txt',\
                    'weather_acceptlist_asrun_30_5sigsym.txt']
        elif patch=='patch_2a':
            keepfs=['weather_acceptlist_2a_asrun_10_5sigsym.txt',\
                    'weather_acceptlist_2a_asrun_30_5sigsym.txt']
        elif patch=='patch_4a':
            keepfs=['weather_4a_acceptlist_asrun_10_5sigma.txt',\
                    'weather_4a_acceptlist_asrun_30_5sigma.txt']
        elif patch=='patch_7b':
            keepfs=['weather_7b_acceptlist_asrun_10_5sigma.txt',\
                    'weather_7b_acceptlist_asrun_30_5sigma.txt']
        # Read Laura's two weather acceptlists
        # NB These files contain lists of diodes to KEEP
        # A scan is not in these files if all its diodes are rejected
        acceptd={}
        for weather_keepf in keepfs:
            acceptd[weather_keepf]=\
                quiet.AcceptedScanList(os.path.join(cutsdir,weather_keepf))

        a1=acceptd[keepfs[0]]
        a2=acceptd[keepfs[1]]

        threshold=0.7
        discount=6 # Modules and diodes already statically removed
        nbefore=len(hamlet.scans.keys())
        total=((quiet.N_MODULE_MAX-2)*quiet.NDIODE_MODULE)-discount
        for (scan,seg) in sorted(acceptlist.scans.keys()):
            ndiodes=0
            for module in range(quiet.N_MODULE_MAX-2):
                for diode in range(quiet.NDIODE_MODULE):
                    # NB boolean test is reversed because of the format of ai
                    if (scan,seg) in a1.scans.keys():
                        keep1=not(a1.isIncluded(scan,seg,module,diode))
                    else:
                        keep1=False # Reject if not in list
                    if (scan,seg) in a2.scans.keys():
                        keep2=not(a2.isIncluded(scan,seg,module,diode))
                    else:
                        keep2=False # Reject if not in list
                    # Check whether the diode was already kept
                    alreadykept=hamlet.isIncluded(scan,seg,module,diode)
                    # Keep track of how many diodes remain for this CES
                    if (keep1 and keep2 and alreadykept):
                        hamlet.includeSegmentDiode(scan,seg,module,diode)
                    else:
                        hamlet.excludeSegmentDiode(scan,seg,module,diode)

    else:
        print 'Invalid cuttype, %i' % cuttype
        return None

        nafter=len(hamlet.scans.keys())
        print 'For cuttype %i, removing CESs with < %3.1f per cent of diodes' \
              % (cuttype,100.0*threshold)
        print '  and that removed %i scans' % (nbefore-nafter)

    hamlet.removeEmptyScans()

    return hamlet

#-----------------------------------------------------------

def EnsureBothQUinModule(acceptlist=None,verbose=False,patch='patch_6a'):

    """
    For an acceptlist, check each CES-module, removing any for which
    there is not at least one Q AND one U diode
    Polarization modules only
    """

    import copy
    import quiet

    # Set up default acceptlist if not user-supplied
    if acceptlist is None:
        rundict,acceptdict=GenerateCMBRunlists()
        acceptlist=acceptdict[patch]

    hamlet   = copy.deepcopy(acceptlist) # Fresh copy

    #Diodes=['Q1','U1','U2','Q2']
    #DiodesLookup={'Q1':0,'U1':1,'U2':2,'Q2':3}

    for (scan,seg) in sorted(acceptlist.scans.keys()):
        for module in range(quiet.N_MODULE_MAX-2):
            noQ=False
            noU=False
            noQ=acceptlist.isExcluded(scan,seg,module,DiodesLookup['Q1']) \
                and acceptlist.isExcluded(scan,seg,module,DiodesLookup['Q2'])
            noU=acceptlist.isExcluded(scan,seg,module,DiodesLookup['U1']) \
                and acceptlist.isExcluded(scan,seg,module,DiodesLookup['U2'])
            if noQ or noU:
                hamlet.excludeModule(module)
                if verbose:
                    print 'Cut  %i.%i module %i' % (scan,seg,module)
            else:
                if verbose:
                    print 'Kept %i.%i module %i' % (scan,seg,module)

    hamlet.removeEmptyScans()

    return hamlet

#-----------------------------------------------------------

def CutStaticCESs(acceptlist=None,verbose=False,patch=None):

    """
    Cut from an acceptlist those CESs that
    feature in the following arrest warrant
    """

    import quiet
    import copy

    if patch is None:
        print 'patch required!'
        return None

    # Set up default acceptlist if not user-supplied
    if acceptlist is None:
        rundict,acceptdict=GenerateCMBRunlists()
        acceptlist=acceptdict[patch]

    hamlet = copy.deepcopy(acceptlist) # Fresh copy
    sa=set(acceptlist.scans.keys())

    # Whole runs:

    # Swapped preamp cables
    cesstart=329
    cesstop=341
    preamps=range(cesstart,cesstop+1) # inclusive

    # Groundscreen door came open
    cesstart=398
    cesstop=398
    screen=range(cesstart,cesstop+1) # inclusive

    # Swapped cables
    cesstart=563
    cesstop=596
    cables=range(cesstart,cesstop+1) # inclusive

    # Deck problems
    cesstart=937
    cesstop=937
    deck1=range(cesstart,cesstop+1) # inclusive

    # Deck problems
    cesstart=953
    cesstop=955
    deck2=range(cesstart,cesstop+1) # inclusive

    # Apply whole-scan cuts:
    wholecut=set(preamps+screen+cables+deck1+deck2)
    for (scan,seg) in sorted(acceptlist.scans.keys()):
        # It would be nice to be able to loop over the shorter lists
        if scan in wholecut:
            hamlet.removeScan(scan,seg)
            if verbose:
                print '%i.%i cut' % (scan,seg)

    # Blown fuse (selected modules)
    #mods=[7,8,9,12,13,16]
    #mods=[1,7,8,9,12,13,16] # UPDATED MODULES 100119
    mods=[7,8,9,12,13,16] # Reverted 100305
    cesstart=633
    cesstop=636
    blownfuse=set(range(cesstart,cesstop+1)) # inclusive
    #blownfuse2=set([(632,5),(637,1)])
    blownfuse2=set([(632,5)]) # Updated 100305

    for (scan,seg) in sorted(acceptlist.scans.keys()):
        if scan in blownfuse or (scan,seg) in blownfuse2:
            hamlet.excludeSegmentModules(scan,seg,mods)
            if verbose:
                print '%i.%i cut modules %s' % (scan,seg,str(mods))

    # Other cuts (see Colin's memo of 28.1.10):
    badlist=[\
        # Receiver off
        (529,2),\
        # 654.0 not moving; Oslo ignores this as a CES
        # Dome maintenance
        (759,1),\
        # Deck problems
        (951,1),\
        # T_cryo stabilizing
        (980,1),\
        # 1413.0 mount stall; Oslo ignores this as a CES
        ]

    hamlet.removeScansIn(badlist)
    if verbose:
        for (scan,seg) in badlist:
            print '%i.%i cut' % (scan,seg)

    hamlet.removeEmptyScans()

    return hamlet

#-----------------------------------------------------------

def fkneeListtoAccept(acceptlist=None,fkneef='fknee_chicagocut.txt',
                      verbose=False,patch=None):

    """
    Cut from an acceptlist those scans that fail the Chicago fknee cut,
    and return/write an acceptlist corresponding to what's been removed
    """

    import quiet
    import os,copy
    from ces_translate import ces_translate

    # Set up default acceptlist if not user-supplied
    if acceptlist is None:
        rundict,acceptdict=GenerateCMBRunlists()
        if patch is None: patch='patch_6a'
        acceptlist=acceptdict[patch]

    cutsdir='/data4/quiet/nulltests/cuts'
    fknees=open(os.path.join(cutsdir,fkneef),'r')

    fk=[]
    fkscans=[]
    for line in fknees:
        if line.startswith('#'): continue
        scan,seg,module,diode=line.split()
        # Translate Chicago -> Oslo scan-numbering scheme
        [(oscan,oseg)]=ces_translate([(int(scan),int(seg))],reverse=False)
        if oscan != -1:
            fk.append((oscan,oseg,int(module),int(diode)))
            fkscans.append((oscan,oseg))

    hamlet=copy.deepcopy(acceptlist)

    for (scan,seg) in sorted(acceptlist.scans.keys()):
        if (scan,seg) in fkscans:
            for module in range(quiet.N_MODULE_MAX):
                for diode in range(quiet.NDIODE_MODULE):
                    if (scan,seg,module,diode) in fk:
                        hamlet.excludeSegmentDiode(scan,seg,module,diode)
                        # And say so:
                        if verbose:
                            print '%i.%i cut %2i %1i for poor f_knee' \
                                  % (scan,seg,module,diode)
        else:
            hamlet.removeScan(scan,seg)

    hamlet.removeEmptyScans()

    outf='fknee.accept'
    hamlet.toFile(os.path.join(cutsdir,outf))

    return hamlet

#-----------------------------------------------------------

def BuildPatchDict(acceptdict):

    """
    For a dictionary of acceptlists keyed by patch,
    return a dictionary of patches keyed by (scan,seg)
    """

    # Now build dictionary of patches keyed by each scan ID
    patchdict={}
    for patch in acceptdict.keys():
        for (scan,seg) in acceptdict[patch].scans.keys():
            patchdict[(scan,seg)]=patch
    return patchdict                    

#-----------------------------------------------------------

def CutSunspikeScans(acceptlist=None,verbose=False,patch=None):

    """
    Synonym wrapper for original (list-based) sun-spike cut
    ***Replaced by more general CutSunSpikes
    """

    CutSunspikes(\
        acceptlist=acceptlist,verbose=verbose,patch=patch,\
        runlist=None,usealgorithm=False)

    return hamlet

#-----------------------------------------------------------

def CutSunspikes(acceptlist=None,verbose=False,patch=None,\
                 runlist=None,modules=None,usealgorithm=True,sunvers=None):

    """
    Cut from an acceptlist those CES-diodes that are near the sun

    Undertake for all modules, unless these are explicitly supplied in a list

    In algorithm mode,
    calculate on the fly

    In list mode,
    CESs are taken from a file (for now)
    Later, will add the actual algorithm to find them
    mysql query was
#mysql> select q19_scan.run_id,q19_scan.run_subid,q19_scan.sun_phi,q19_scan.sun_theta from q19_scan WHERE q19_scan.sun_phi BETWEEN 120 AND 180 AND q19_scan.sun_theta BETWEEN 44 AND 60;

    In list mode #2,
    CESs are taken from a file provided by Yuji

    """

    import quiet
    import os,copy,scanf
    from ces_translate import ces_translate

    if (sunvers is not None):
        uselist2=True
    else: uselist2=False
    assert(not (usealgorithm and uselist2)), 'Invalid choice of sun cut!'

    # Patch required
    assert(patch is not None), 'patch is required for %s!' % whoami()

    # Set up default acceptlist if not user-supplied
    if runlist is None:
        rundict,acceptdict=GenerateCMBRunlists()
        runlist=rundict[patch]
        acceptlist=acceptdict[patch]

    hamlet = copy.deepcopy(acceptlist) # Fresh copy
    sa=set(acceptlist.scans.keys())

    sunrejected={}
    if not usealgorithm and not uselist2:
        print '** (using CES list for sun-in-sidelobe cuts)'
        # Read the scans from file
        cutsdir='/data4/quiet/nulltests/cuts'
        badf='sunspike_rejected_by_chicago.txt'
        f=open(os.path.join(cutsdir,badf),'r')

        for line in f:
            if line.startswith('#'): continue
            (scan,seg,sun_phi,sun_theta) = line.split()
            # Translate Chicago -> Oslo scan-numbering scheme
            [(oscan,oseg)]=ces_translate([(int(scan),int(seg))],reverse=False)
            if oscan != -1:
                sunrejected[oscan,oseg]=[float(sun_phi),float(sun_theta)]
                # if verbose: print oscan,oseg,float(sun_phi),float(sun_theta)

        f.close()

    if not usealgorithm and uselist2:
        assert(modules is None), "**Why have you set [modules]??"
        assert(sunvers is not None), 'This should not happen!'
        print '** (using LIST (v%i) for sun-in-sidelobe cuts)' % sunvers
        cutsdir='/usit/titan/u1/jonathaz/sunspike/from-yuji-100623/'
        badf='q19_sun_spike.version%i.dat' % sunvers
        f=open(os.path.join(cutsdir,badf),'r')
#        verbose=True
        linef='%i %i %i %i'
        for line in f:
            if line.startswith('run_id'): continue
            try:
                (scan,seg,mod,has_spike)=scanf.sscanf(line,linef)
            except scanf.IncompleteCaptureError:
                print 'Skipped %s' % line
                continue
            # Translate Chicago -> Oslo scan-numbering scheme
            [(oscan,oseg)]=ces_translate([(scan,seg)],reverse=False)
            if oscan != -1:
                if ((oscan,oseg) in sa and has_spike > 0):
                    hamlet.excludeSegmentModules(oscan,oseg,mod)
                    if verbose:
                        print '%i.%i RQ_%02i cut (spike %i) %s' \
                              % (oscan,oseg,mod,has_spike,\
                                 str(hamlet.isIncluded(oscan,oseg,mod,0)))
        # Working here / CHECK / TEST HERE
        f.close()
        hamlet.removeEmptyScans()
        return hamlet

    elif usealgorithm:
        print '** (using algorithm for sun-in-sidelobe cuts)'
        # Loop over CESs
        for (scan,seg) in sorted(acceptlist.scans.keys()):
            # call FetchPointing
            #az[(scan,seg)],el[(scan,seg)],dk[(scan,seg)],\
            #    mjdstart[(scan,seg)],mjdstop[(scan,seg)]=\
            #    FetchPointing(scan,seg,runlist,patch)
            #mjd=(mjdstart[(scan,seg)],mjdstop[(scan,seg)])/2.0
            #theta,phi=CalculateSunCoords(az[(scan,seg)],el[(scan,seg)],\
            #                             dk[(scan,seg)],mjd,verbose=verbose)
            theta,phi=CalculateCESSunCoords(\
                scan,seg,runlist=runlist,patch=patch,\
                verbose=False,translate=False)

            # Populate sunrejected
            # These limits are from Yuji's memo of 9 Feb 2010
            cut_lunch   = 44.0 < theta < 60.0 and 120.0 < phi < 180.0
            cut_dinner1 = 44.0 < theta < 60.0 and -75.0 < phi < -15.0
            cut_dinner2 = 60.0 < theta < 65.0 and -90.0 < phi <  0.0
            cut = cut_lunch or cut_dinner1 or cut_dinner2
            #cut = cut_lunch
            if cut: sunrejected[(scan,seg)]=[theta,phi]

    # Run through the acceptlist and apply the cut
    for (scan,seg) in sorted(sunrejected.keys()):
        if (scan,seg) in sa:
            if modules is None:
                hamlet.removeScan(scan,seg)
                if verbose:
                    print '%i.%i cut' % (scan,seg)
            elif type(modules)==type(list()):
                #modules=[0,1,2]
                for module in modules:
                    hamlet.excludeSegmentModules(scan,seg,module)
                    if verbose:
                        print '%i.%i RQ%2i cut' % (scan,seg,module)
            else:
                print 'modules argument in incorrect format %s', str(modules)
                return None

    hamlet.removeEmptyScans()

    return hamlet

#-----------------------------------------------------------

def CalculatePatchSunCoords(acceptlist=None,runlist=None,patch=None,\
                          verbose=False,translate=False,outfile=None):

    """
    For a given patch, for an acceptlist,
    return suntheta,sunphi keyed by CES
    """

    assert(patch is not None), 'patch is required for %s!' % whoami()

    # Set up runlist if not user-supplied
    if runlist is None:
        rundict,acceptdict=GenerateCMBRunlists()
        runlist=rundict[patch]
        if acceptlist is None:
            acceptlist=acceptdict[patch]

    tofile=False
    if outfile is not None:
        outf=open(outfile,'w')
        outf.write('# CES suntheta sunphi\n')
        tofile=True # Could use dir()

    coords={}
    for (scan,seg) in sorted(acceptlist.scans.keys()):
        suntheta,sunphi=CalculateCESSunCoords(\
            scan,seg,runlist=runlist,patch=patch,\
            verbose=verbose,translate=False)
        coords[(scan,seg)]=[suntheta,sunphi]
        #print '%i.%i: %3.2f %3.2f' % (scan,seg,suntheta,sunphi)
        if tofile:
            line='%i.%i %3.2f %3.2f\n' % (scan,seg,suntheta,sunphi)
            outf.write(line)

    if tofile: outf.close()

    return coords

#-----------------------------------------------------------

def CalculateCESSunCoords(scan,seg,runlist=None,patch=None,\
                          verbose=False,translate=False):

    """
    For a given CES, return sun theta,phi coordinates
    """

    assert(patch is not None), 'patch is required for %s!' % whoami()

    # Set up runlist if not user-supplied
    if runlist is None:
        rundict,acceptdict=GenerateCMBRunlists()
        runlist=rundict[patch]

    if translate:
        from ces_translate import ces_translate
        oldscan=scan
        oldseg=seg
        [(scan,seg)]=ces_translate([(scan,seg)])
        if scan == -1:
            print 'Scan %i.%i not available in Oslo' % (scan,seg)
            return
    else:
        oldseg=seg

    try:
        type='cmb'
        scan in runlist.targetTypes[type].targets[patch].scans.keys()
    except KeyError:
        print 'scan not in runlist, %s %s %s' % (scan,patch,type)
        return None

    (az,el,dk,mjd1,mjd2)=FetchPointing(scan,seg,runlist,patch=patch)
    mjd=(mjd1+mjd2)/2.0

    if verbose:
        print 'mjd (start, stop): %f (%f, %f)' % (mjd,mjd1,mjd2)

    suntheta,sunphi=CalculateSunCoords(az,el,dk,mjd,verbose=verbose)

    return suntheta,sunphi

#-----------------------------------------------------------

def CalculateSunCoords(az,el,dk,mjd,verbose=False,returnsunpos=False):

    """
    For a given telescope az, el, dk, and an MJD,
    Return boresight theta and phi coordinates for Sun
    """

    # Fetch source info

    # MJD -> AZsource,ELsource
    # NB sunaz,sunel are pyephem angle objects
    [sunaz,sunel]=CalculatePlanetPositions(mjd,verbose=False)['Sun']

    # Arguments and outputs in radians
    suntheta,sunphi=Horiz2BoresightSystem(\
        deg2rad(az),deg2rad(el),deg2rad(dk),float(sunaz),float(sunel),
        verbose=verbose)

    #if verbose: print 'theta,phi', rad2deg(suntheta),rad2deg(sunphi)

    if returnsunpos:
        return rad2deg(suntheta),rad2deg(sunphi),\
               ang2deg(sunaz),ang2deg(sunel)

    return rad2deg(suntheta),rad2deg(sunphi)

#-----------------------------------------------------------

def ang2rad(ang):

    """
    Convert pyephem angle object to radians
    Really just a reminder that float(ang) is the operation ang2rad(ang)
    """

    return float(ang)


#-----------------------------------------------------------

SidelobeCoordsTriple={'theta':52.0,'phi':150.0}

#-----------------------------------------------------------

def exclusive(seta,setb):
    "Test whether two sets are mutually exclusive"
    return (seta & setb) == set()

#-----------------------------------------------------------

def CalculateSidelobeElevation(az,el,dk,theta,phi,degrees=False):

    """
    Calculate the evelation of a sidelobe, given:
    AZ,EL,DK of boresight
    theta,phi coordinates of sidelobe
    All angles in radians
    Optionally take and return angles in degrees
    """

    from math import sin,cos,asin

    if degrees:
        az=deg2rad(az)
        el=deg2rad(el)
        dk=deg2rad(dk)
        theta=deg2rad(theta)
        phi=deg2rad(phi)

    el_sidelobe = asin(sin(el)*cos(theta) \
                       + cos(el)*cos(dk)*sin(theta)*cos(phi) \
                       - cos(el)*sin(dk)*sin(theta)*sin(phi))

    if degrees:
        el_sidelobe=rad2deg(el_sidelobe)

    return el_sidelobe

#-----------------------------------------------------------

def Horiz2BoresightSystem(az,el,dk,azsource,elsource,verbose=False,\
                          degrees=False):

    """
    At a given MJD, take AZ, EL and DK coords of telescope in radians
    Return thetasource, phisource in radians
    Optionally take and return angles in degrees
    """

    from math import cos,sin,acos,atan2,atan,pi
    from numpy import cross,dot

    if degrees:
        az=deg2rad(az)
        el=deg2rad(el)
        dk=deg2rad(dk)
        azsource=deg2rad(azsource)
        elsource=deg2rad(elsource)

    # Some scans fail these assertions...
    strict=False
    if strict:
        assert(0.0 <= el <= pi/2.0), el
        assert(deg2rad(-135.0) <= az <= deg2rad(265.0)), az
        assert(deg2rad(-265.0) <= dk <= deg2rad(176.0)), dk
    assert(-2.0*pi <= azsource <= 2.0*pi), azsource
    assert(-pi/2.0 <= elsource <= pi/2.0), elsource


    #print rad2deg(az),rad2deg(el),rad2deg(dk),rad2deg(azsource),rad2deg(elsource)

    #if degrees:
    #    for angle in (az,el,dk,azsource,elsource):
    #        angle=deg2rad(angle)

    # Notation from memo of 28 January 2010
    # NB p0, ps should have a minus in the second component

    p0=[cos(az)*cos(el),-sin(az)*cos(el),sin(el)]
    s0=[-cos(az)*sin(el)*cos(dk)-sin(az)*sin(dk),\
        sin(az)*sin(el)*cos(dk)-cos(az)*sin(dk),\
        cos(el)*cos(dk)]
    r0=list(cross(p0,s0))

    v=[cos(azsource)*cos(elsource),-sin(azsource)*cos(elsource),sin(elsource)]

    # Calculate the separation of source and boresight
    thetasource=acos(float(dot(v,p0)))
    phisource=atan2(float(dot(v,r0)),float(dot(v,s0)))

    if verbose:
        print 'az,el,dk,azs,els/deg',rad2deg(az),rad2deg(el),rad2deg(dk),\
              rad2deg(azsource),rad2deg(elsource)
        print 'p0',p0
        print 's0',s0
        print 'r0',r0
        print 'v',v
        print 'theta,phi/deg',rad2deg(thetasource),rad2deg(phisource)

    assert(abs(phisource) <= pi), 'phisource (%3.2f) out of range' % phisource

    if degrees:
        thetasource=rad2deg(thetasource)
        phisource=rad2deg(phisource)

    return thetasource,phisource

#-----------------------------------------------------------

def CutTypeBBadness(acceptlist=None,badnesses=None,\
                    verbose=False,patch=None):

    """
    Cut from an acceptlist those CES-diodes that have poor enough
    Type-B parameters to be in the read-in badness file
    """

    import quiet
    import os,copy

    assert(patch is not None), 'patch is required for %s!' % whoami()

    # Set up acceptlist if not user-supplied
    if acceptlist is None:
        rundict,acceptdict=GenerateCMBRunlists()
        acceptlist=acceptdict[patch]

    # If no badnesses supplied, read them from file
    if badnesses is None:
        cutsdir='/data4/quiet/nulltests/cuts'
        badf='typeb_badness_outofrange_0.01_10.0.txt' # In file because slow
        badnesses=ReadTypeBBadness(os.path.join(cutsdir,badf))
    bs=set(badnesses.keys())

    hamlet = copy.deepcopy(acceptlist) # Fresh copy
    # Run through the acceptlist
    for (scan,seg) in sorted(acceptlist.scans.keys()):
        # If that scan is listed in the list of poor type-B fits:
        if (scan,seg) in bs:
            # Run through all badness diode entries for that scan:
            for (modnum,diode,precorr,postcorr) in badnesses[(scan,seg)]:
                # If that diode entry is in the acceptlist,
                # and the badness is not < 0.01 (i.e., is > 10):
                # (though Read's minthresh may have stripped out < 0.01)
                if hamlet.isIncluded(scan,seg,modnum,diode) and \
                       postcorr >= 1.0:
                    # Then exclude the entry:
                    hamlet.excludeSegmentDiode(scan,seg,modnum,diode)
                    # And say so:
                    if verbose:
                        print '%i.%i cut %2i %1i for poor fit %4.1f' \
                              % (scan,seg,modnum,diode,postcorr)

    hamlet.removeEmptyScans()

    return hamlet

#-----------------------------------------------------------

def ReadTypeBBadness(filename,thresh=None,minthresh=None,verbose=False):

    """
    Read Type-B badnesses from file, to be cut
    Convert to Oslo-style CES numbering scheme
    Optionally, only when post-correction Chi^2 is below thresh
    Default to cutting out Chi^2 BELOW some minthresh (see below)
    The file comes from qufs SQL (query at top of file)
    Chi^2 < 0.01 has removed anything with a zero value for that quantity
    Chi^2 > 10.0 is what Chicago uses for its Chi^2 cut
    """

    from ces_translate import ces_translate

    if minthresh is None: minthresh=0.01 # To remove missing Chi^2 values

    file=open(filename)

    badness={}
    print 'reading (slow) %s for type-B Chi^2s' % filename

    # If thresh not supplied
    if thresh is None:
        for line in file:
            if line.startswith('#'): continue
            scan,seg,module,diode,precorr,postcorr = line.split()
            if postcorr >= minthresh:
                # Translate Chicago -> Oslo scan-numbering scheme
                [(oscan,oseg)]=ces_translate([(int(scan),int(seg))],reverse=False)
                if oscan != -1:
                    if (oscan,oseg) not in badness.keys():
                        badness[oscan,oseg]=[]
                    di=DiodesLookup[diode]
                    badness[oscan,oseg].append(\
                        (int(module),di,float(precorr),float(postcorr)))
                    if verbose: print oscan,oseg,int(module),di,postcorr
                else:
                    pass

    else:
        # Otherwise only read if postcorr >= supplied threshold
        for line in file:
            if line.startswith('#'): continue
            scan,seg,module,diode,precorr,postcorr = line.split()
            if float(postcorr) <= thresh and float(postcorr) >= minthresh:
                # Translate Chicago -> Oslo scan-numbering scheme
                [(oscan,oseg)]=\
                    ces_translate([(int(scan),int(seg))],reverse=False)
                if oscan != -1:
                    if (oscan,oseg) not in badness.keys():
                        badness[oscan,oseg]=[]
                    di=DiodesLookup[diode]
                    badness[oscan,oseg].append((\
                        int(module),di,float(precorr),float(postcorr)))
                    if verbose: print oscan,oseg,int(module),di,postcorr
            else:
                pass

    file.close()
    print '(reading is complete)'

    return badness

#-----------------------------------------------------------

def FetchTypeBBadness(arg):

    """
    """

    import os
    import ces_translate

    #thresh=10.0

    # Read goodnesses-of-fit from file, if not done already
    if type(arg)==type(''):
        # This step is slow because the file is huge
        badnesses_chicago=ReadTypeBBadness(file)
    else:
        badnesses_chicago=arg

    # Test case
    #scan=414
    #seg=0
    #[(newscan,newseg)]=\
    #     ces_translate.ces_translate([(scan,seg)],reverse=False)
    #print newscan,newseg

    # In progress

    print badnesses_chicago.keys()[0]
    badnesses={}
    for (scan,seg) in badnesses_chicago.keys():
        [(newscan,newseg)]=\
             ces_translate.ces_translate([(scan,seg)],reverse=False)
        if scan != -1:
            if [(scan,seg)] in badnesses.keys():
                badnesses
        else:
            pass
   #     if scan==-1:
   #         del badnesses_chicago[(oldscan,oldseg)]
   #     else:
   #         pass
#            print scan,seg

    print badnesses.keys()[0]

    return
    #badnesses=badnesses_chicago
    badnesses={}
    return badnesses

#-----------------------------------------------------------

def QueryPointing(scan,seg,patch=None,runlist=None,translate=False):

    """
    Print/return the pointing centres of a given scan
    """

    assert(patch is not None), 'patch is required for %s!' % whoami()

    from ephem import minute

    if translate:
        from ces_translate import ces_translate
        oldscan=scan
        oldseg=seg
        [(scan,seg)]=ces_translate([(scan,seg)])
        if scan == -1:
            print 'Scan %i.%i not available in Oslo' % (scan,seg)
            return
    else:
        oldseg=seg

    if runlist is None:
        rundict,accdict=GenerateCMBRunlists()
        runlist=rundict[patch]

    az,el,dk,mjdstart,mjdstop=\
        FetchPointing(scan,seg,runlist,patch,type='cmb')

    scanlength=(mjdstop-mjdstart)/minute

    print '%i.%i (.%i chic) %s: %3.1f %3.1f %3.1f %f %f (%3.1f mins)' \
          % (scan,seg,oldseg,patch,\
             az,el,dk,mjdstart,mjdstop,scanlength)

    return az,el,dk,mjdstop,mjdstart

#-----------------------------------------------------------

def FetchSkyAngles(verbose=True,acceptlist=None,runlist=None,patch=None):

    """
    Calculate distances from Sun, Moon for each scan,seg in an acceptlist
    """

    # Set up default runlist and acceptlist if not user-supplied
    if runlist is None:
        rundict,acceptdict=GenerateCMBRunlists()
        runlist=rundict[patch]
        if acceptlist is None:
            acceptlist=acceptdict[patch]

    az={}
    el={}
    dk={}
    mjdstart={}
    mjdstop={}
    objects={}
    # Fetch pointing and date information from the runlist
    for scan,seg in sorted(acceptlist.scans.keys()):
        az[(scan,seg)],el[(scan,seg)],dk[(scan,seg)],\
            mjdstart[(scan,seg)],mjdstop[(scan,seg)]=\
            FetchPointing(scan,seg,runlist,patch)
        mjd0=mjdstart[(scan,seg)]
        mjd1=(mjdstart[(scan,seg)]+mjdstop[(scan,seg)])/2.0
        mjd2=mjdstop[(scan,seg)]
        for nmjd,mjd in enumerate([mjd0,mjd1,mjd2]):
            obj=CalculatePlanetPositions\
                (mjd,verbose,az[(scan,seg)],el[(scan,seg)])
            if verbose: print scan,seg,nmjd,obj
            objects[(scan,seg,nmjd)]=obj

    if verbose:
        print 'This acceptlist (for %s) has %i CESs' \
            % (patch,acceptlist.acount())
        for (scan,seg,nmjd) in sorted(objects.keys()):
            print '%4i %1i %1i %3.1f %4.1f %4.1f' \
                  % (scan,seg,nmjd,objects[(scan,seg,nmjd)]['El'],\
                     ang2deg(objects[(scan,seg,nmjd)]['Sunsep']),\
                     ang2deg(objects[(scan,seg,nmjd)]['Moonsep']),)

    return objects

#-----------------------------------------------------------

def CalculatePlanetPositions(mjd,verbose=False,az=None,el=None):

    """
    Calculate positions of Sun and Moon (and others, to be added later)
    given the MJD
    If az and el are supplied, return separations
    """

    import ephem

#+++ quiet/oslo/src/f90/include/quiet_pointing_mod.f90   2009-12-06 05:08:24 UTC (rev 787)
#!  real(dp), parameter :: QUIET_GEODETIC_LONGITUDE  = -67.76166667d0  ! 67 degrees, 45 minutes, 42.0 seconds WEST
#!  real(dp), parameter :: QUIET_GEODETIC_LATITUDE   = -23.02822222d0 ! 23 degrees, 1 minute, 41.6 seconds SOUTH
#!  real(dp), parameter :: QUIET_ALTITUDE            =  5020.d0        ! meters
#real(dp), parameter :: QUIET_GEODETIC_LONGITUDE  = -67.7667d0      ! Keith's longitude
#  real(dp), parameter :: QUIET_GEODETIC_LATITUDE   = -23.0333d0      ! Keith's latitude
#real(dp), parameter :: QUIET_ALTITUDE            =  5080.d0        ! Keith's altitude

# real(dp), parameter :: QUIET_GEOCENTRIC_LATITUDE = -22.890103588d0 ! degrees
# real(dp), parameter :: QUIET_RADIUS              = 1.000277634d0   ! units of earth radii

    CO = ephem.Observer()
    CO.long = '-67.7667'
    CO.lat = '-23.0333'
    CO.elevation = 5080.0 # Assumed to be metres
    #CO.date='2010/01/01 19:35:16'
    djd=mjd2djd(mjd)
    CO.date=ephem.Date(djd)

    s=ephem.Sun(CO)
    m=ephem.Moon(CO)

    #s=ephem.Sun()
    #m=ephem.Moon()

    #s.compute(mjd)
    #m.compute(mjd)

    if verbose:
        print 'Sun RADEC: ', s.ra, s.dec
        print 'Moon RADEC:', m.ra, m.dec

        print 'Sun AZEL:', s.az,s.alt
        print 'Moon AZEL:', m.az,m.alt

    # Set up the observatory
    objects={}

    objects['Lat']    = CO.lat # Latitude
    objects['Long']   = CO.long # Longitude
    objects['Height'] = CO.elevation # Height
    objects['Date']   = mjd # MJD
    objects['Sun']    = [s.az,s.alt] # Sun AZ, EL
    objects['Moon']   = [m.az,m.alt] # Moon AZ, EL
    objects['Sepn']   = ephem.separation(s,m) # Separation of Moon and Sun...

    # Keep the pyephem objects in case we need them
    objects['Obs']    = CO
    objects['M']      = m
    objects['S']      = s

    if az is not None:
        objects['Az'] = az
    if el is not None:
        objects['El'] = el

    if az is not None and el is not None:
        objects['Sunsep'],objects['Moonsep']=\
            SpecificSepns(az,el,objects['S'],objects['M'])

    return objects

#-----------------------------------------------------------

def ang2deg(ang):

    """
    Convert a pyephem angle object into a float angle in degrees
    """

    from math import pi

    return (float(ang)*180.0/pi)

#-----------------------------------------------------------

    #if test:
    #    q=45.0
    #    q==secateur.ang2deg(secateur.deg2ang(q))

#-----------------------------------------------------------

def deg2ang(deg):

    """
    Convert an angle (in degrees) into a pyephem angle object [radians]
    """

    from math import pi
    from ephem import degrees

    return degrees(deg*pi/180.0)

#-----------------------------------------------------------

def deg2rad(deg):

    """
    Convert float angle in degrees into radians
    """

    from math import pi

    return deg*pi/180.0

#-----------------------------------------------------------

def rad2ang(rad):

    """
    Convert float angle in radians into pyephem angle
    """

    return deg2ang(rad2deg(rad))

#-----------------------------------------------------------

def rad2deg(rad):

    """
    Convert float angle in radians into degrees
    """

    from math import pi

    return rad*180.0/pi

#-----------------------------------------------------------

def SpecificSepns(az,el,sun,moon,verbose=False):

    """
    Return CES angular separations from Sun and Moon
    az,el come from e.g. secateur.FetchPointing()
    sun, moon objects come from e.g. secateur.CalculatePlanetPositions()
    """

    import ephem

    # Fairly unrobust check of inputs
    # Need to strengthen, and add check on sun and moon arguments
    if az is None:
        print 'az is not set in calling routine'
        return 1
    if el is None:
        print 'el is not set in calling routine'
        return 2
    else:
        pass

    # Calculate separations via pyephem
    sunsep=sepn((deg2ang(az),deg2ang(el)),(sun.az,sun.alt))
    moonsep=sepn((deg2ang(az),deg2ang(el)),(moon.az,moon.alt))

    # And via our own trig:
    sunsep2=sepn2(deg2rad(az),deg2rad(el),float(sun.az),float(sun.alt))
    moonsep2=sepn2(deg2rad(az),deg2rad(el),float(moon.az),float(moon.alt))

    if float(sunsep) != sunsep2:
        # Should be a machine EPS test
        print 'Ambiguity1', sunsep, sunsep2
        return 3

    if float(sunsep) != sunsep2:
        # Should be a machine EPS test
        print 'Ambiguity2', sunsep, sunsep2
        return 4

    if verbose:
        print 'Via pyephem: ',sunsep,moonsep
        print 'Via sph trig:',rad2ang(sunsep2),rad2ang(moonsep2)

    return sunsep,moonsep

#-----------------------------------------------------------

def sepn(ang1,ang2):

    """
    Calculate pyEphem angle object between two loci on a sphere
    ***ang1 and ang2 must be pyephem angle objects
    *** Or I think they can be a sequence of coord pairs
    """

    from ephem import separation

    # To be got working - boils down to modifying of list in place..
    #if type(ang1[0]) != type(ephem.degrees(0.0)):
    #    for ang in ang1:
    #        ang=deg2ang(ang)
    #else:
    #    pass

    #if type(ang2[0]) != type(ephem.degrees(0.0)):
    #    ang2=deg2ang(ang2)
    #else:
    #    pass
    
    return separation(ang1,ang2)

#-----------------------------------------------------------

def FetchPointing(scan,seg,runlist,patch=None,type='cmb'):

    """
    Extract pointing information for a ces from a L2 runlist
    """

    import quiet

    assert(patch is not None), 'patch is required for %s!' % whoami()

    #scan,seg=str(ces).split('.')

    #for scan,seg in sorted(runlist[patch].scans.keys()):
    #if type(runlist)==type(quiet.L2Runlist): # How to test?

    # Check arguments robustly:
    #if type(runlist)==type(dict()): # str object not callable error..!
    #    r=runlist[patch]
    #    runlist=r # Not sure this trick will work!

    #if type(scan)==type(''):
    #    scan=int(scan)

#    if scan not in runlist.targetTypes[type].targets[patch].scans.keys():
#        print 'scan not in runlist,', scan
#        return None

    try:
        scan in runlist.targetTypes[type].targets[patch].scans.keys()
    except KeyError:
        print 'scan not in runlist,', scan, patch, type
        return None

    #if type(seg)==type(''):
    #    seg=int(seg)

    nseg = runlist.targetTypes[type].targets[patch].scans[scan].nSegment
    
    assert (seg >= 1), 'segment %i out of range (< 1)' % seg
    assert (seg <= nseg), 'segment %i out of range (> %i)' % (seg,nseg)

    #segid = runlist.targetTypes[type].targets[patch].scans[scan].id

    nsg=seg-1
    #print seg,nsg
    #print runlist.targetTypes[type].targets[patch].scans[scan].segments[0].index
    #print runlist.targetTypes[type].targets[patch].scans[scan].segments[:].index(seg)

    #print type,patch,scan,seg
    #print sorted(runlist.targetTypes[type].targets[patch].scans.keys())

    az = runlist.targetTypes[type].targets[patch].scans[scan].segments[nsg].az
    el = runlist.targetTypes[type].targets[patch].scans[scan].segments[nsg].el
    dk = runlist.targetTypes[type].targets[patch].scans[scan].segments[nsg].dk
    mjdstart = runlist.targetTypes[type].targets[patch].scans[scan].segments[nsg].startDate
    mjdstop = runlist.targetTypes[type].targets[patch].scans[scan].segments[nsg].endDate
    id = runlist.targetTypes[type].targets[patch].scans[scan].segments[nsg].id

    #print nseg

    #if unfixaz: az = 360.0 - az # Sigurd emails @ 100113

    return az, el, dk, mjdstart, mjdstop

#-----------------------------------------------------------

def FetchAPEXData(UTdate):

    """
    Fetch one day's worth of APEX data
    """

    import string,time,ephem

    #if ('/' in UTdate): # Can't get this syntax right, so just ignore
        #startdate = string.replace(UTdate,'-','/')
    startdate = UTdate
    stopdate  = float(ephem.Date(startdate)) + 1.0 # Stop at the next day

    startdate=ephem.Date(startdate)
    stopdate=ephem.Date(stopdate-0.5*ephem.minute) # To avoid repeated midnight

    print 'Fetching', startdate, '->', stopdate

    dataList = APEX_weather(startdate, stopdate)

    return dataList

#-----------------------------------------------------------

def ReadAPEXData(filename="/usit/titan/u1/jonathaz/quiet_data/calib_data/jonathaz/apex/apex.dat"):

    """
    Read APEX data from file into list of lists
    """

    import numpy

    #catalogue={}

    lines=ReadList(filename)
    #print 'LL',lines[0]

    apex_data=[]
    for index,line in enumerate(lines):
        words=re.split('\s+',line)

        #print line
        #[date,time,mjd,pwv,dew_point,humidity,pressure,T_apex,\
        # wind_dirn,wind_vel]\
        # = '%s %s %f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f' %(line)

        date      = str(words[0])    #
        time      = str(words[1])    #
        mjd       = float(words[2])  #
        pwv       = float(words[3])  #
        dew_point = float(words[4])  #
        humidity  = float(words[5])  #
        pressure  = float(words[6])  #
        T_apex    = float(words[7])  #
        wind_dirn = float(words[8])  #
        wind_vel  = float(words[9])  #

        date_string = date+' '+time
        #catalogue[mjd]=pwv
        #apex_data.append([mjd,date_string,pwv])
        apex_data.append([mjd,pwv])

    #print_dict(catalogue)

    #print numpy.array(apex_data,dtype=numpy.float64)[0][0]
    #return numpy.array(apex_data)

    return apex_data

#-----------------------------------------------------------

def WriteAPEXData(band,filename):

    """
    Write APEX data to file
    NEED TO HANDLE GAPS
    """

    import os,ephem,math

    # Check for and open file
    if os.path.exists(filename):
        os.rename(filename,filename+'.old')

    output=open(filename,'w')

    # Write data file header
    line = '# Time (UTC) (MJD)'
    for key in APEX_weather_keyList:    
        line += '  %s (%s)' %(key, APEX_weather_unit[key])
    #print line
    output.write(line+'\n')

    #UTdate = '2009-05-05'

    startdate,stopdate = QUIETBookends(band)
    startdate = math.floor(float(ephem.Date(startdate))) # float optional here
    stopdate = math.ceil(float(ephem.Date(stopdate))) # float optional here

    #stopdate=startdate+1.0
    #print 'Fetching', startdate , '->', stopdate

    # Write date ranges to data file
    line = '# '+str(ephem.Date(startdate))+' -> '+str(ephem.Date(stopdate))
    output.write(line+'\n')

    # Fetch APEX data, one day at a time, and write to data file
    for n in xrange(int(startdate),int(stopdate)+1):
        #print ephem.Date(n-0.5)
        dataList = FetchAPEXData(ephem.Date(n-0.5)) # 0.5 needed for midnight

        for data in dataList:
            #print data
            line = time.strftime('%Y/%m/%d %H:%M:%S', data['tm_tuple'])
            line+=' '+ '%f' %(ephem.julian_date(line))
            for key in APEX_weather_keyList:
                #line += '\t%6.2f [%s]' %(data[key], APEX_weather_unit[key])
                line += '\t%6.2f' %(data[key])
            #print line
            output.write(line+'\n')

    # Close data file
    output.close()

    return

#-----------------------------------------------------------

# Define static cuts on data
# List is on wiki
# cf, Chicago wiki list (to be done)
# Times to be updated eventually using Laura's parsing of CESs
# !!! * = needs to be checked !!!
StaticData=\
             [('2008/11/28 16:40', '2008/11/30 17:20', '329-341', 'all', 'all', 'Swapped preamps for testing','[Dorothea]'),\
('2009/01/09 11:00', '2009/01/15 20:00', '563-593','all', 'all', 'Swapped cables', '[Rob]'),\
('2009/02/04 15:20', '2009/02/04 17:40', 'none', 'none', 'none', 'Did Side-lobe tests', '[Glen]'),\
('2009/02/14 06:00', '2009/02/16 18:22', 'none', 'none', 'none', 'SJIM full (a few data files from this time, any zero-data files not processed to L1)','[]'),\
('2009/03/26 16:20', '2009/03/30 00:35', 'none', 'none', 'none', 'Two instances of ADC servers running','[]'),\
('2009/04/24 14:30', '2009/04/24 20:30', 'none', 'none', 'none', '* Running wrong version of RCS', '[Glen]'),\
('1998/01/01', '1998/01/02', '632.5-636', '7-9,12-13,16', 'all', '* Blown fuse', '[]')]

#-----------------------------------------------------------

# Numbering scheme to be checked [JTLZ]
DeadModulesQ=[('RQ16 all diodes',16,'all'),('RQ08 U2',8,2)]

#-----------------------------------------------------------

def CutDead(outdir='.',file='dead.cut',band='Q'):

    """
    Cut dead modules and diodes
    """

    import os

    mjdstart,mjdstop=QUIETBookends(band,'MJD')

    file = os.path.join(outdir,file)
    if os.path.exists(file):
        os.rename(file,file+'.old')

    outfile=OpenFile(file)

    WriteHeader(file,'dead','cut',outfile)

    patch = 'all'
    ces = 'all'
    for corpse in DeadModulesQ:
        modulelist = str(corpse[1])
        diodelist = str(corpse[2])
        comment = corpse[0]
        WriteCutLine(outfile,mjdstart,mjdstop,ces,patch,modulelist,\
                     diodelist,comment)

    outfile.close()
    print '  Cuts written to', file


    return

#-----------------------------------------------------------


def CutStatic(outdir='.',StaticData=StaticData,file='static.cut'):

    """
    Cut data based on standard static cuts
    """

    import os

    file = os.path.join(outdir,file)
    if os.path.exists(file):
        os.rename(file,file+'.old')

    outfile=OpenFile(file)

    WriteHeader(file,'static','cut',outfile)

    #ces = 'all'
    #modulelist='all'
    #diodelist='all'
    #print file
    patch = 'all'

    for static in StaticData:
        mjdstart=UTC2MJD(static[0])
        mjdstop=UTC2MJD(static[1])
        ces=static[2]
        modulelist=static[3]
        diodelist=static[4]
        comment=static[-1]+' '+static[-2]
        WriteCutLine(outfile,mjdstart,mjdstop,ces,patch,modulelist,\
                     diodelist,comment)

    outfile.close()
    print '  Cuts written to', file

    return


#-----------------------------------------------------------

def CutPWV(cuts,outdir='.',patch='all',modulelist='all',diodelist='all'):

    "Cut data based on PWV threshold(s)"

    import os, numpy

    cutargs=re.split('\s+',cuts)

    pwv_data=ReadAPEXData()
    minthresh=0.0   # mm
    maxthresh=100.0 # mm

    pwv_data=array(pwv_data)

    #print pwv_data[:,0],pwv_data[:,1]

    if len(cutargs) <= 3:
        # --cut out data with pwv greater than minthresh
        if cutargs[0:2] == ['pwv','>']:
            minthresh=float(cutargs[-1]) # cut threshold in mm
            cutcriterion = (pwv_data[:,1]>minthresh)
        # --cut out data with pwv less than maxthresh
        if cutargs[0:2] == ['pwv','<']:
            maxthresh=float(cutargs[-1]) # cut threshold in mm
            cutcriterion = (pwv_data[:,1]<maxthresh)
    elif len(cutargs) > 3:
        # --keep data for which pwv in range given
        if cutargs[1:4] == ['<','pwv','<']:
            minthresh=float(cutargs[0])
            maxthresh=float(cutargs[-1])
            keepcriterion = (minthresh<pwv_data[:,1]<maxthresh)
        # --keep data for which pwv in range given
        if cutargs[1:4] == ['>','pwv','>']:
            minthresh=float(cutargs[-1])
            maxthresh=float(cutargs[0])
            keepcriterion = (minthresh<pwv_data[:,1]<maxthresh)

    #pwv_mask=where(criterion,1,0)
    #diff_mask=numpy.diff(pwv_mask[:])
    #print pwv_mask[1:1000]

    # Differentiate mask and identify non-zero elements (i.e. changes)
    #if cutcriterion: # HOW TO TEST FOR ARRAY EXISTENCE IN NUMPY?
    diff_mask=numpy.diff(where(cutcriterion,1,0))
    #elif keepcriterion:
    #    diff_mask=numpy.diff(where(keepcriterion,0,1)) # CHECK THIS

    steps=diff_mask.nonzero()

    #print steps
    #print diff_mask[600:800]
    #print pwv_data[600:800,:]
    #print pwv_data[0,0],pwv_data[0,1]
    #print pwv_data[1,0],pwv_data[1,1]

    #mask[:,1]=diff_mask[:]*pwv_data[:,1]

    #mask=numpy.zeros((len(pwv_data[:,0]),2))
    mask = numpy.where(cutcriterion,1,0)#[:,numpy.newaxis]

    #for index,value in enumerate(mask[:]):
    #    print index,mask[index],pwv_data[index,1]

    file = os.path.join(outdir,'pwv_'+str(minthresh)+'_'+str(maxthresh)+'.cut')
    if os.path.exists(file):
        os.rename(file,file+'.old')

    outfile=OpenFile(file)

    WriteHeader(file,'pwv','cut',outfile)

    patch = patch
    ces = 'all'
    modulelist='all'
    diodelist='all'
    comment=cuts

    nblocks=0
    if mask[0] > 0:
        print 'Starting cut'
        nblocks+=1
        #print 'cut %12.5f' %(pwv_data[0,0]),
        mjdstart = pwv_data[0,0]
    else:
        print 'Starting kept'

    for index,value in enumerate(mask[1:-1]):
        if mask[index] > 0 and mask[index-1] < 1:
            nblocks+=1
            #print 'cut %12.5f' %(pwv_data[index,0]),
            mjdstart = pwv_data[index,0]
            if mask[index+1] < 1: # Handle lone cut samples
                #print '%12.5f' %(pwv_data[index,0])
                mjdstop = pwv_data[index,0]
                WriteCutLine(outfile,mjdstart,mjdstop,ces,patch,modulelist,\
                             diodelist,comment)

        elif mask[index] > 0 and mask[index+1] < 1:
            #print '%12.5f' %(pwv_data[index,0])
            mjdstop = pwv_data[index,0]
            WriteCutLine(outfile,mjdstart,mjdstop,ces,patch,modulelist,\
                         diodelist,comment)

    if mask[-1] > 0:
        #print '%12.5f' %(pwv_data[-1,0])
        mjdstop = pwv_data[-1,0]
        WriteCutLine(outfile,mjdstart,mjdstop,ces,patch,modulelist,\
                     diodelist,comment)
        print 'Finishing cut'
    else:
        print 'Finishing kept'

    print '[Cut', nblocks, 'contiguous blocks]'

    outfile.close()
    print '  Cuts written to', file

    return pwv_data

#-----------------------------------------------------------

def CutHK(doparallel,cuts,outdir='.',patch='all',modulelist='all',\
          diodelist='all'):

    "Cut data based on L1 HK (T_encl, T_cryo) threshold(s)"

    import os, numpy

    verbose=False
    #verbose=True

    cutargs=re.split('\s+',cuts)
    # Strip out the hk tag:
    cutargs=cutargs[1:]

    data_rootdir='/data4/quiet'

    data={}
    data['T2'],data['T3']=FetchMultiTData(doparallel,data_rootdir)

    # At the moment this has been going round the loop more than once...
    # Bug - being called once per np
    # Use np=1? (but then no point in MPI..)

    # Sort T entries into time order...
    # Reshape T dictionaries and split into mjd,T2 and mjd,T3
    #T2_vec=[]
    #T3_vec=[]
    #for f in sorted(T2_data.keys()):
    #    #print f
    #    T2_vec.append(T2_data[f][:])
    #    T3_vec.append(T3_data[f][:])

    #print 'T data:'
    #for sample in T2_vec[:]:
    #    for el in sample
    #        print 's:',sample[0],sample[1]
    #print len(T2_vec[:])
    #print T2_vec

    # Then calculate the MJDs where the threshold is passed
    # Then write it to file a la PWV

    T2minthresh=0.0   # K
    T2maxthresh=100.0 # K

    T3minthresh=0.0   # K
    T3maxthresh=100.0 # K

    # Split the requested temperature cuts up
    
    #f='level1/2008/10/B-101808164442-0.fits'
    types={}
    thresh={}
    if 't2' in cutargs:
        Ttype='T2'
        print cutargs
        if len(cutargs) > 3:
            types[Ttype]=cutargs[0:len(cutargs)/2] # T2 comes first
        else:
            types[Ttype]=cutargs
        thresh[Ttype]=[T2minthresh,T2maxthresh]
        #data[Ttype]=array(T2_data[f])
    if 't3' in cutargs:
        Ttype='T3'
        #T_data=array(T3_data[f])
        print cutargs
        if len(cutargs) > 3:
            types[Ttype]=cutargs[len(cutargs)/2:] # and T3 second...
        else:
            types[Ttype]=cutargs
        thresh[Ttype]=[T3minthresh,T3maxthresh]
        #data[Ttype]=array(T3_data[f])

    #print T3_data[f]
    # Select average required for cutting
    median=True
    mean=not(median)
    if median: avcol=1
    elif mean: avcol=2

    files={}
    outfiles={}
    filesopened={}

    # Process the thresholds and write the cut lines
    for Ttype in types:
        cutargs=types[Ttype]
        print 'c',cutargs
        minthresh,maxthresh=thresh[Ttype]
        print minthresh,maxthresh
        filesopened[Ttype]=False

        # Loop over L1 files
        for L1File in sorted(data[Ttype].keys()):
            T_data=array(data[Ttype][L1File])

            if len(cutargs) <= 3:
                # --cut out data with T greater than minthresh
                if cutargs[0:2] == [Ttype.lower(),'>']:
                    minthresh=float(cutargs[-1]) # cut threshold in mm
                    cutcriterion = (T_data[:,avcol]>minthresh)
                # --cut out data with T less than maxthresh
                if cutargs[0:2] == [Ttype.lower(),'<']:
                    maxthresh=float(cutargs[-1]) # cut threshold in mm
                    cutcriterion = (T_data[:,avcol]<maxthresh)
            elif len(cutargs) > 3:
                # This doesn't seem ever to have worked? Need to check
                # --keep data for which T in range given
                if cutargs[1:4] == ['<',TType.lower(),'<']:
                    minthresh=float(cutargs[0])
                    maxthresh=float(cutargs[-1])
                    keepcriterion = (minthresh<T_data[:,avcol]<maxthresh)
                # --keep data for which T in range given
                if cutargs[1:4] == ['>',TType.lower(),'>']:
                    minthresh=float(cutargs[-1])
                    maxthresh=float(cutargs[0])
                    keepcriterion = (minthresh<T_data[:,avcol]<maxthresh)

            # Differentiate mask and identify non-zero elements (i.e. changes)
            # if cutcriterion: # HOW TO TEST FOR ARRAY EXISTENCE IN NUMPY?
            diff_mask=numpy.diff(where(cutcriterion,1,0))
            steps=diff_mask.nonzero()
            mask = numpy.where(cutcriterion,1,0)#[:,numpy.newaxis]

            # Set up the output files (once per Ttype))
            print 'm',minthresh,maxthresh
            if not filesopened[Ttype]:
                files[Ttype] = os.path.join(outdir,Ttype+'_'+str(minthresh)+'_'+str(maxthresh)+'.cut')
                file=files[Ttype]

                if os.path.exists(files[Ttype]):
                    os.rename(files[Ttype],files[Ttype]+'.old')

                outfiles[Ttype]=OpenFile(files[Ttype])
                outfile=outfiles[Ttype]
                extras=cuts+' on '+str(sorted(data[Ttype].keys()))
                WriteHeader(files[Ttype],extras,'cut',outfile)
                filesopened[Ttype]=True

            outfile=outfiles[Ttype]
            patch = patch
            ces = 'all'
            modulelist='all'
            diodelist='all'
            comment=cuts+' '+StripExtension(os.path.basename(L1File))[0]
            print L1File, comment

            nblocks=0
            if mask[0] > 0:
                print 'Starting cut'
                nblocks+=1
                if verbose: print 'cut %12.5f' %(T_data[0,0]),
                mjdstart = T_data[0,0]
            else:
                print 'Starting kept'

            for index,value in enumerate(mask[1:-1]):
                if mask[index] > 0 and mask[index-1] < 1:
                    nblocks+=1
                    if verbose: print 'cut %12.5f' %(T_data[index,0]),
                    mjdstart = T_data[index,0]
                    if mask[index+1] < 1: # Handle lone cut samples
                        if verbose: print '%12.5f' %(T_data[index,0])
                        mjdstop = T_data[index,0]
                        WriteCutLine(outfile,mjdstart,mjdstop,ces,patch,modulelist,\
                                     diodelist,comment)

                elif mask[index] > 0 and mask[index+1] < 1:
                    if verbose: print '%12.5f' %(T_data[index,0])
                    mjdstop = T_data[index,0]
                    WriteCutLine(outfile,mjdstart,mjdstop,ces,patch,modulelist,\
                                 diodelist,comment)

            if mask[-1] > 0:
                if verbose: print '%12.5f' %(T_data[-1,0])
                mjdstop = T_data[-1,0]
                WriteCutLine(outfile,mjdstart,mjdstop,ces,patch,modulelist,\
                             diodelist,comment)
                print 'Finishing cut'
            else:
                print 'Finishing kept'

            print '[Cut', nblocks, 'contiguous blocks]'

    # Close the .cut files
    for Ttype in types:
        outfiles[Ttype].close()
        print '  Cuts written to', files[Ttype]

    return data

#-----------------------------------------------------------

def ReadTPData(runlist):

    "Read TP data for a runlist of L2 CESs"

    print "TP multiple fetch not yet implemented"
    tp_data = None
    return None

    return tp_data

#-----------------------------------------------------------

def CutTP(cuts,patch='all',modulelist='all',diodelist='all'):

    """
    Cut data based on TP threshold(s)
    """

    import os, numpy

    var = 'tp'

    cutargs=re.split('\s+',cuts)

    # Need to change this for TP data
    print "TP cut not yet implemented"
    return None
    tp_data=ReadTPData(runlist) # NB NOT REALLY TP DATA YET

    minthresh=0.0   # K
    maxthresh=100.0 # K

    tp_data=array(tp_data)

    #print pwv_data[:,0],pwv_data[:,1]

    if len(cutargs) <= 3:
        # --cut out data with tp greater than minthresh
        if cutargs[0:2] == [var,'>']:
            minthresh=float(cutargs[-1]) # cut threshold in mm
            cutcriterion = (tp_data[:,1]>minthresh)
        # --cut out data with tp less than maxthresh
        if cutargs[0:2] == [var,'<']:
            maxthresh=float(cutargs[-1]) # cut threshold in mm
            cutcriterion = (tp_data[:,1]<maxthresh)
    elif len(cutargs) > 3:
        # --keep data for which tp in range given
        if cutargs[1:4] == ['<',var,'<']:
            minthresh=float(cutargs[0])
            maxthresh=float(cutargs[-1])
            keepcriterion = (minthresh<tp_data[:,1]<maxthresh)
        # --keep data for which tp in range given
        if cutargs[1:4] == ['>',var,'>']:
            minthresh=float(cutargs[-1])
            maxthresh=float(cutargs[0])
            keepcriterion = (minthresh<tp_data[:,1]<maxthresh)

    diff_mask=numpy.diff(where(cutcriterion,1,0))

    steps=diff_mask.nonzero()

    mask = numpy.where(cutcriterion,1,0)

    #for index,value in enumerate(mask[:]):
    #    print index,mask[index],pwv_data[index,1]

    file=os.path.join(outdir,var+'_'+str(minthresh)+'_'+str(maxthresh)+'.cut')
    if os.path.exists(file):
        os.rename(file,file+'.old')

    outfile=OpenFile(file)

    WriteHeader(file,'TP','cut',outfile)

    patch = patch
    ces = 'all'
    modulelist='all'
    diodelist='all'
    comment=cuts

    nblocks=0
    if mask[0] > 0:
        print 'Starting cut'
        nblocks+=1
        #print 'cut %12.5f' %(pwv_data[0,0]),
        mjdstart = tp_data[0,0]
    else:
        print 'Starting kept'

    for index,value in enumerate(mask[1:-1]):
        if mask[index] > 0 and mask[index-1] < 1:
            nblocks+=1
            #print 'cut %12.5f' %(tp_data[index,0]),
            mjdstart = tp_data[index,0]
            if mask[index+1] < 1: # Handle lone cut samples
                #print '%12.5f' %(tp_data[index,0])
                mjdstop = tp_data[index,0]
                WriteCutLine(outfile,mjdstart,mjdstop,ces,patch,modulelist,\
                             diodelist,comment)

        elif mask[index] > 0 and mask[index+1] < 1:
            #print '%12.5f' %(tp_data[index,0])
            mjdstop = tp_data[index,0]
            WriteCutLine(outfile,mjdstart,mjdstop,ces,patch,modulelist,\
                         diodelist,comment)

    if mask[-1] > 0:
        #print '%12.5f' %(tp_data[-1,0])
        mjdstop = tp_data[-1,0]
        WriteCutLine(outfile,mjdstart,mjdstop,ces,patch,modulelist,\
                     diodelist,comment)
        print 'Finishing cut'
    else:
        print 'Finishing kept'

    print '[Cut', nblocks, 'contiguous blocks]'

    outfile.close()
    print '  Cuts written to', file

    return tp_data

#-----------------------------------------------------------

def SplitTime(acceptlist,flavour='season',outdir='.',writefiles=True):

    """
    Jackknife L2 acceptlist in time by CES, and output to acceptlist files
    Inelegant!
    Operates on dictionary of acceptlists keyed by patch
    (outputted by e.g. GenerateCMBRunlists()).
    NB Jackknife halves are not currently forced to be the same length.
    Optionally write acceptlists to file
    """

    import quiet, copy

    assert(writefiles or \
           (len(acceptlist.keys())<=1) and not writefiles), \
           "Can't run silently for more than one patch"

    backup=copy.deepcopy(acceptlist) # NB acceptlist is a dict

    for patch in acceptlist.keys():

        sharks = copy.deepcopy(acceptlist[patch]) # Fresh copy
        jets   = copy.deepcopy(acceptlist[patch]) # Fresh copy

        # Jackknife by alternate CESs
        if flavour=='interleaved':
            toggle=False
            for scan,seg in sorted(acceptlist[patch].scans.keys()):
                if toggle: sharks.removeScan(scan,seg)
                toggle=not(toggle)

            toggle=True
            for scan,seg in sorted(acceptlist[patch].scans.keys()):
                if toggle: jets.removeScan(scan,seg)
                toggle=not(toggle)
            flavour='inter'

        # Jackknife by alternate pairs of CESs
        # NOT YET IMPLEMENTED
        if flavour=='chequerboard':
            print flavour,' not yet implemented'
            return 'notyet', 'notyet'
            #sys.exit(0)
#            toggle=False
#            for n,(scan,seg) in enumerate(sorted(acceptlist[patch].scans.keys())):
#                if toggle: sharks.removeScan(scan,seg)
#                toggle=not(toggle)

#            toggle=True
#            for n,(scan,seg) in enumerate(sorted(acceptlist[patch].scans.keys())):
#                if toggle: jets.removeScan(scan,seg)
#                toggle=not(toggle)

        # Jackknife whole season in two contiguous halves of CES
        if flavour=='season':
            nscans=len(acceptlist[patch].scans.keys())
            for n,(scan,seg) in enumerate(sorted(acceptlist[patch].scans.keys())):
                if n < nscans/2: sharks.removeScan(scan,seg)
            #for n,(scan,seg) in enumerate(sorted(acceptlist[patch].scans.keys())):
                if (scan,seg) in sharks.scans.keys():
                #if n >= nscans/2:
                    jets.removeScan(scan,seg)

        # Jackknife whole season by alternate quarters of CESs
        if flavour=='battenberg':
            nscans=len(acceptlist[patch].scans.keys())
            for n,(scan,seg) in enumerate(sorted(acceptlist[patch].scans.keys())):
                # JZ This line needs to be robustified!
                # if nscans % 4 == 0, < is ok before 3; otherwise <=!
                if n < (nscans-1)/4 or ((nscans-1)/2 <= n <= 3*(nscans-1)/4):
                    sharks.removeScan(scan,seg)
            #for n,(scan,seg) in enumerate(sorted(acceptlist[patch].scans.keys())):
                if (scan,seg) in sharks.scans.keys():
                #if (nscans/4 <= n <= nscans/2) or n >= 3*nscans/4:
                    jets.removeScan(scan,seg)
            flavour='quart'
            print nscans
        # Jackknife whole season by quarters of CESs
        # The first and final quarters are paired against the middle pair
        # (a la tandem rigging in a rowing four)
        if flavour=='tandem':
            nscans=len(acceptlist[patch].scans.keys())
            for n,(scan,seg) in enumerate(sorted(acceptlist[patch].scans.keys())):
                if n < nscans/4 or (n >= 3*nscans/4):
                    sharks.removeScan(scan,seg)
            #for n,(scan,seg) in enumerate(sorted(acceptlist[patch].scans.keys())):
                if (scan,seg) in sharks.scans.keys():
                #if (nscans/4 <= n <= 3*nscans/4):
                    jets.removeScan(scan,seg)

        # Check pairs are mutually exclusive in their CESs
        assert(exclusive(set(sharks.scans.keys()),set(jets.scans.keys())))

        # Check pairs are within one CES of each other in length
        ls=len(sharks.scans.keys())
        lj=len(jets.scans.keys())
        assert(abs(ls-lj)<=1),\
           'jackknife pairs should be within one CES in length of each other: %i /= %i' % (ls,lj)

        # Write splits to file, using title to carry name
        splitname=flavour
        sharks.title = '%s_%s.a' % (patch,splitname)
        jets.title   = '%s_%s.b' % (patch,splitname)

        # Assert original acceptlist is unchanged
        assert(acceptlist==backup), 'acceptlist changed!'

        if writefiles:
            return WriteSplits(sharks,jets,patch,outdir)
        else:
            return sharks,jets

    return

#-----------------------------------------------------------

def ListToHistogram(list,nbins=10,shift=False):

    """
    Return a histogram from a list of bins,
    as two lists: (i) frequencies and (ii) bins
    """

    import numpy

    array = numpy.array(list)

    if shift:
        zeropoint=array[0]
        array-=zeropoint
    else:
        zeropoint=0.0

    (n,bins) = numpy.histogram(array,bins=nbins,normed=False)

    return n,bins,zeropoint

#-----------------------------------------------------------

def SplitBy(acceptlist,quantity,label='',unit='',diodes=False):

    """
    Take an acceptlist, and split it into two after sorting by the quantity
    Quantity must be a dictionary keyed by (scan,seg)
    Includes assertions on exclusivity and approx. length match
    Currently assumes all CESs in the acceptlist feature in the quantity list
    If diodes is True, allow split by CES-diode rather than CES
    """

    import copy
    from operator import itemgetter

    # Remove any irrelevant entries in the quantity list
    for key in quantity.keys():
        (scan,seg)=key[0:2]
        if (scan,seg) not in acceptlist.scans.keys(): del quantity[key]

    sharks = copy.deepcopy(acceptlist) # Fresh copy
    jets   = copy.deepcopy(acceptlist) # Fresh copy

    # Sort the quantity list by ascending value
    quantitys=sorted(quantity.items(), key=itemgetter(1))
    nkeys=len(quantitys)
    for n,(key,qu) in enumerate(quantitys):
        if not diodes:
            (scan,seg)=key; nscans=nkeys;
            if n < nscans/2:
                sharks.removeScan(scan,seg)
            if (scan,seg) in sharks.scans.keys():
                jets.removeScan(scan,seg)
        elif diodes:
            (scan,seg,mod,di)=key; ndiodes=nkeys;
            if n < ndiodes/2:
                sharks.excludeSegmentDiode(scan,seg,mod,di)
            elif n >= ndiodes/2:
                jets.excludeSegmentDiode(scan,seg,mod,di)

    quantonly=[qq for (key,qq) in quantitys]
    print '%s range: %3.2f -> %3.2f, median ~ %3.2f /%s' \
          % (label,quantonly[0],quantonly[-1],quantonly[int(nkeys/2)],unit)

    # For intra-CES splits, need to add some assertions
    # assert
    # assert
    # assert
    if diodes:
        # I think this will fail once cuts are in place
        # -> Need to mask using input acceptlist
        from numpy import bitwise_and, bitwise_not
        for (scan,seg) in acceptlist.scans.keys():
            masked_sharks=bitwise_and(acceptlist.scans[(scan,seg)],\
                                  sharks.scans[(scan,seg)])
            masked_jets_not=bitwise_and(acceptlist.scans[(scan,seg)],\
                                  bitwise_not(jets.scans[(scan,seg)]))
            assert(masked_sharks.all()==masked_jets_not.all()),\
                'discrepancy in %i.%i (\n%s \nv \n%s)' % \
                (scan,seg,str(masked_sharks),str(masked_jets_not))
        print 'Warning: skipping jackknife assertions'
        sharks.removeEmptyScans(); jets.removeEmptyScans();
        return sharks, jets

    # Check pairs are mutually exclusive in their CESs
    assert(exclusive(set(sharks.scans.keys()),set(jets.scans.keys())))

    # Check pairs are within one CES of each other in length
    assert(abs(len(sharks.scans.keys())-len(jets.scans.keys()))<=1),\
           'jackknife pairs should be within one CES in length of each other'

    return sharks,jets

#-----------------------------------------------------------

def WriteSplits(sharks,jets,patch,outdir):

    """
    Write jackknife halves to file
    """

    sharkfile = os.path.join(outdir,sharks.title)
    jetfile   = os.path.join(outdir,jets.title)

    if os.path.exists(sharkfile):
        os.rename(sharkfile,'%s.old'%sharkfile)
    if os.path.exists(jetfile):
        os.rename(jetfile,'%s.old'%jetfile)

    sharkheader=WriteHeader(sharkfile,patch,'jk')
    jetsheader=WriteHeader(sharkfile,patch,'jk')

    sharks.toFile(sharkfile,sharkheader)
    jets.toFile(jetfile,jetsheader)

    print 'Jackknives written to %s.x' % '.'.join(sharkfile.split('.')[:-1])

    return sharkfile,jetfile

#-----------------------------------------------------------

def SplitGain(acceptlist,flavour=None,outdir='.',writefiles=False):

    """
    !!! LOGIC TO BE CHECKED IF ANY OTHER FLAVOURS ADDED
    Jackknife L2 acceptlist in gain proxy by CES
    Operates on dictionary of acceptlists keyed by patch
    (outputted by e.g. GenerateCMBRunlists()).
    Jackknife halves are checked to be within one in length of each other
    Optionally write acceptlist halves to file
    """

    import copy

    assert(writefiles or \
           (len(acceptlist.keys())<=1) and not writefiles), \
           "Can't run silently for more than one patch at a time"

    assert(flavour in ['reg']), '** Invalid jackknife flavour!'

    # Jackknife by pwv high v low
    if flavour=='reg':
        # Fetch the T_enc regulating v not fraction
        reg_frac = ReadRegulatingFraction()

    backup=copy.deepcopy(acceptlist)
    hamlet=copy.deepcopy(acceptlist) # NB acceptlist is a dict of acceptlists

    for patch in acceptlist.keys():
        # Jackknife by regulating fraction, 1 v not (having removed -1's)
        if flavour=='reg':
            splitname='tencreg'
            # >  half the values are 100% regulating, so do these against
            # the others
            for (scan,seg) in acceptlist[patch].scans.keys():
                # For now just remove any scans for which no value
                if reg_frac[(scan,seg)] < 0: # i.e. nan
                    hamlet[patch].removeScan(scan,seg)

            sharks = copy.deepcopy(hamlet[patch]) # Fresh copy
            jets   = copy.deepcopy(hamlet[patch]) # Fresh copy
            for (scan,seg) in hamlet[patch].scans.keys():
                if reg_frac[(scan,seg)] < 1.0:
                    sharks.removeScan(scan,seg)
                else:
                    jets.removeScan(scan,seg)

            # Check pairs are mutually exclusive in their CESs
            assert(exclusive(set(sharks.scans.keys()),set(jets.scans.keys())))

            # Check pairs are within one CES of each other in length
            # [ Invalid assertion in this case - threshold division ]

        # Write splits to file, using title to carry name
        sharks.title = '%s_%s.a' % (patch,splitname)
        jets.title   = '%s_%s.b' % (patch,splitname)

        # Check acceptlist has not been modified
        assert(acceptlist==backup), 'acceptlist modified!'

        if writefiles:
            return WriteSplits(sharks,jets,patch,outdir)
        else:
            return sharks,jets

    return

#-----------------------------------------------------------

def SplitFrequency(acceptlist,ces_stats=EXAMPLE_CES_STATS,\
                   flavour=None,outdir='.',writefiles=False,skipTT=True):

    """
    Operates on dictionary of acceptlists keyed by patch
    (outputted by e.g. GenerateCMBRunlists()).
    Jackknife L2 acceptlist in ces_stats by CES, and output to acceptlist files
    Jackknife halves are checked to be within one in length of each other
    Optionally write acceptlist halves to file
## WORKING HERE    Optionally skip the TTs
    """

    import os,copy

    assert(writefiles or \
           (len(acceptlist.keys())<=1) and not writefiles), \
           "Can't run silently for more than one patch at a time"

    assert(flavour in ['fknee','chisq_ss']), '** Invalid jackknife flavour!'

    # Fetch the ces_stats
    #ces_stats = ReadCESStats(f=EXAMPLE_CES_STATS)
    #statsdir='/data4/quiet/jonathaz4/xx'

    backup=copy.deepcopy(acceptlist)

    for patch in acceptlist.keys():
        ifknee=4; ichisqss=6; ichisqhigh=8;
        fknees={}; chisqss_pos={}; chisqss_neg={};
        chisqhigh_pos={}; chisqhigh_neg={};
        # Retrieve CES-diode stats
        for (scan,seg,mod,di) in ces_stats.keys():
            if skipTT and mod > 16: continue
            # But skip any diodes already excluded from jackknife
            if acceptlist[patch].isExcluded(scan,seg,mod,di): continue
            fknees[(scan,seg,mod,di)]=ces_stats[(scan,seg,mod,di)][ifknee]
            chisqss=ces_stats[(scan,seg,mod,di)][ichisqss]
            chisqhigh=ces_stats[(scan,seg,mod,di)][ichisqhigh]
            if chisqss >= 0.0:
                chisqss_pos[(scan,seg,mod,di)]=chisqss
            elif chisqss < 0.0:
                chisqss_neg[(scan,seg,mod,di)]=chisqss
            if chisqhigh >= 0.0:
                chisqhigh_pos[(scan,seg,mod,di)]=chisqhigh
            elif chisqhigh < 0.0:
                chisqhigh_neg[(scan,seg,mod,di)]=chisqhigh

        # Jackknife by fknee high v low
        if flavour=='fknee':
            splitname='fknee'
            quantity=fknees; label=flavour; unit='Hz'
        # Jackknife by chisq_ss_pos high v low
        elif flavour=='chisqsspos':
            splitname='chisqsspos'
            quantity=chisqss_pos; label=flavour; unit='sigma'
        # Jackknife by chisq_ss_neg high v low
        elif flavour=='chisqssneg':
            splitname='chisqssneg'
            quantity=chisqss_neg; label=flavour; unit='sigma'
        # Jackknife by chisq_high_pos high v low
        elif flavour=='chisqhighpos':
            splitname='chisqhighpos'
            quantity=chisqhigh_pos; label=flavour; unit='sigma'
        # Jackknife by chisq_high_neg high v low
        elif flavour=='chisqhighneg':
            splitname='chisqhighneg'
            quantity=chisqhigh_neg; label=flavour; unit='sigma'

        # Assertions are handled inside SplitBy
        sharks,jets=SplitBy(acceptlist[patch],quantity,label,unit,diodes=True)

        # Write splits to file, using title to carry name
        sharks.title = '%s_%s.a' % (patch,splitname)
        jets.title   = '%s_%s.b' % (patch,splitname)

        # Check acceptlist has not been modified
        assert(acceptlist==backup), 'acceptlist modified!'

        if writefiles:
            return WriteSplits(sharks,jets,patch,outdir)
        else:
            return sharks,jets

    return

#-----------------------------------------------------------

def SplitWeather(acceptlist,flavour=None,outdir='.',writefiles=False):

    """
    Operates on dictionary of acceptlists keyed by patch
    (outputted by e.g. GenerateCMBRunlists()).
    Jackknife L2 acceptlist in weather by CES, and output to acceptlist files
    Jackknife halves are checked to be within one in length of each other
    Optionally write acceptlist halves to file
    """

    import copy

    assert(writefiles or \
           (len(acceptlist.keys())<=1) and not writefiles), \
           "Can't run silently for more than one patch at a time"

    assert(flavour in ['pwv','humidity']), '** Invalid jackknife flavour!'

    # Fetch the weather summary data
    apex_summary = ReadAPEXSummary()

    hamlet={}
    backup=copy.deepcopy(acceptlist)

    for patch in acceptlist.keys():
        hamlet[patch]=copy.deepcopy(acceptlist[patch])
        # Jackknife by pwv high v low
        if flavour=='pwv':
            splitname='apexpwv'
            index=0; label=flavour.upper(); unit='mm'
        # Jackknife by humidity high v low
        elif flavour=='humidity':
            splitname='apexhum'
            index=2; label=flavour; unit='%'

        quantity={}
        for (scan,seg) in acceptlist[patch].scans.keys():
            # For now just remove any scans for which no APEX data
            # Alternative is to include + assume poor weather/high PWV/hum
            # Another is to include but check using *local* humidity
            quantity_mean = apex_summary[(scan,seg)][index]
            if quantity_mean < 0: # i.e. nan
                hamlet[patch].removeScan(scan,seg)
            else:
                quantity[(scan,seg)] = quantity_mean

        # Assertions are handled inside SplitBy
        sharks,jets=SplitBy(hamlet[patch],quantity,label,unit)

        # Write splits to file, using title to carry name
        sharks.title = '%s_%s.a' % (patch,splitname)
        jets.title   = '%s_%s.b' % (patch,splitname)

        # Check acceptlist has not been modified
        assert(backup==acceptlist), 'acceptlist modified by %s!' % flavour

        if writefiles:
            return WriteSplits(sharks,jets,patch,outdir)
        else:
            return sharks,jets

    return

#-----------------------------------------------------------

def SplitPointing(acceptlist,runlist,flavour='elevation',outdir='.',\
                  writefiles=True):

    """
    acceptlist and runlist args are dictionaries of those keyed by patch
    Jackknife L2 acceptlist in pointing by CES, and output to acceptlist files
    Operates on dictionary of acceptlists keyed by patch
    (outputted by e.g. GenerateCMBRunlists()).
    NB Jackknife halves are not currently forced to be the same length.
    Jackknife halves are checked to be within one in length of each other
    Optionally write acceptlist halves to file
    """

    # Useful! http://code.activestate.com/recipes/425397
    # [seq[i::num] for i in range(num)]

    import quiet, copy
    from operator import itemgetter
    from numpy import array

    assert(writefiles or \
           (len(acceptlist.keys())<=1 and not writefiles)), \
           "Can't run silently for more than one patch at a time"

    backup=copy.deepcopy(acceptlist)

    buckets={}
    for patch in acceptlist.keys():

        sharks = copy.deepcopy(acceptlist[patch]) # Fresh copy
        jets   = copy.deepcopy(acceptlist[patch]) # Fresh copy

        az={}; el={}; dk={}
        if flavour=='elsidelobe': elsidelobe={}

        # Extract all pointing info from runlist but discard mjd info
        for (scan,seg) in sorted(acceptlist[patch].scans.keys()):
            az[(scan,seg)],el[(scan,seg)],dk[(scan,seg)],mjdstart,mjdstop=\
                FetchPointing(scan,seg,runlist[patch],patch)
            #print scan,seg,el[(scan,seg)]
            if flavour=='elsidelobe':
                elsidelobe[scan,seg]=CalculateSidelobeElevation(\
                    az[(scan,seg)],el[(scan,seg)],dk[(scan,seg)],\
                    SidelobeCoordsTriple['theta'],SidelobeCoordsTriple['phi'],\
                    degrees=True)

        if flavour=='elsidelobe':
            splitname='pointside'
            sels=sorted(elsidelobe.items(), key=itemgetter(1))
            nscans=len(sels)
            for n,((scan,seg),el) in enumerate(sels):
                if n < nscans/2:
                    sharks.removeScan(scan,seg)
                if (scan,seg) in sharks.scans.keys():
                    jets.removeScan(scan,seg)
            selsonly=[ees for ((scan,seg),ees) in sels]
            print 'Sidelobe El range: %2.0f -> %2.0f, median ~ %2.0f /deg' % (selsonly[0],selsonly[-1],selsonly[int(nscans/2)])
            
        # Jackknife by elevation high v low
        if flavour=='elevation':
            splitname='pointel'
            # Sort the acceptlist by ascending el into a list called els
            els=sorted(el.items(), key=itemgetter(1))
            nscans=len(els)
            for n,((scan,seg),el) in enumerate(els):
                if n < nscans/2:
                    sharks.removeScan(scan,seg)
                if (scan,seg) in sharks.scans.keys():
                    jets.removeScan(scan,seg)
            elsonly=[ee for ((scan,seg),ee) in els]
            print 'Elevation range: %2.0f -> %2.0f, median ~ %2.0f /deg' \
                  % (elsonly[0],elsonly[-1],elsonly[int(nscans/2)])

        # Jackknife by elevation azimuth east v west
        if flavour=='azimuth':
            splitname='pointaz'
            # Sort the acceptlist by ascending az into a list called azs
            azs=sorted(az.items(), key=itemgetter(1))
            nscans=len(azs)
            for n,((scan,seg),az) in enumerate(azs):
                if n < nscans/2:
                    sharks.removeScan(scan,seg)
                if (scan,seg) in sharks.scans.keys():
                    jets.removeScan(scan,seg)
            azsonly=[aa for ((scan,seg),aa) in azs]
            print 'Azimuth range: %2.0f -> %2.0f, median ~ %2.0f /deg' \
                  % (azsonly[0],azsonly[-1],azsonly[int(nscans/2)])

        if flavour == 'sunsep':
            # Working here
            splitname='pointsun'
            obj=secateur.FetchSkyAngles(acceptlist=acceptlist[patch],\
                runlist=runlist[patch],patch=patch,verbose=False)
            
        # Jackknife by deck angle (90-deg phase)
        if flavour=='deck90':
            splitname='pointdk90'
            phase=90.0 # DK test phase in deg
            dk0=0.0 # DK offset in deg
            tol=5.0 # Tolerance in target
            # Colin @ 100412
            # Think fans
            for (scan,seg) in sorted(acceptlist[patch].scans.keys()):
                if abs((dk[(scan,seg)] % (2.0*phase) - dk0)) < phase:
                    sharks.removeScan(scan,seg)
                else:
                    jets.removeScan(scan,seg)

        # Jackknife by deck angle (45-deg phase)
        if flavour=='deck45':
            splitname='pointdk45'
            phase=45.0 # DK test phase in deg
            dk0=0.0 # DK offset in deg
            tol=-5.0 # Tolerance in target (affects 225deg only -> agree w/ Ch)
            # Colin @ 100412
            # Think fans
            for (scan,seg) in sorted(acceptlist[patch].scans.keys()):
                if abs((dk[(scan,seg)] % (2.0*phase) - dk0)) < phase+tol:
                    sharks.removeScan(scan,seg)
                else:
                    jets.removeScan(scan,seg)

        # Jackknife by deck angle
        if flavour=='deck1':
            splitname='pointdk1'

            # Sort the acceptlist by ascending dk into a list called dks
            dks=sorted(dk.items(), key=itemgetter(1))

            # Produce a second list without the scans
            dksonly=[dd for ((scan,seg),dd) in dks]

            verbose=True
            #verbose=False
            # Get a rough list of angles for this patch
            if verbose:
                angles = list(set([('%2d' % int(round(dd,-1))) for dd in dksonly]))
                angles.sort(key=int)
                print patch,angles

            # Bin the DKs to establish what angles to jackknife
            nbins=array([-180,-135,-90,-45,0,45,90,135,180,225])
            # nbins=8 # Or set the number of bins
            freqs,bins,zp=ListToHistogram(dksonly,nbins,shift=False)

            # Set minimum number of scans in each bin
            # This leaves four angles for all patches but three for 7a!
            #7b was not included in this though
            if patch != 'patch_7a': # 7a is the exception...
                thresh=50
            else:
                thresh=0
            bucket=[]
            for nbin,freq in enumerate(freqs):
                if freq > thresh:
                    bucket.append(nbin)
                buckets[patch]=bucket

            if verbose:
                print '# nbin n bin: high low centre width'
                for i in buckets[patch]:
                    print 'x',i+1,freqs[i],bins[i],bins[i+1],\
                          (bins[i+1]+bins[i])/2.0,(bins[i+1]-bins[i])
                # for nbin,freq in enumerate(freqs):
                #    if freq > thresh:
                #        print '%i %3i % 4.2f % 4.2f % 4.2f % 4.2f' % (nbin+1,freq,\
                #            bins[nbin],bins[nbin-1],bins[nbin]+zp,bins[nbin-1]+zp)
                # Use of round is very rough in this context

            # Make the jackknife
            # Now we have the bins, need to identify scans in bins
            perm=[1,2]
            # Working here (NB FixPerms for dirs)
            j=0 # .a = first bin
            j2=1 # .b = second bin
            bin=buckets[patch]
            bin2=buckets[patch]
            nscans=len(dk)
            for n,((scan,seg),dk) in enumerate(dks):
                print dk, bin[j], bin[j+1]
                # Working here - to fix
                if dk < bin[j] or dk > bin[j+1]:
                    sharks.removeScan(scan,seg)
                if dk < bin2[j] or dk > bin2[j+1]:
                    jets.removeScan(scan,seg)

        # Jackknife by deck angle
        if flavour=='deck2':
            splitname='pointdk2'
            return 'Not implemented!'
            inner=[4,5,8,9,10,13,14]
            outer=[0,1,2,3,6,7,11,12,15]
            a=inner
            b=outer

        # Jackknife left v right
        if flavour=='deck3':
            splitname='pointdk3'
            return 'Not implemented!'
            left=[7,8,9,10,11,12,13,14,15]
            right=[0,1,2,3,4,5,6]
            a=left
            b=right

        # Check pairs are mutually exclusive in their CESs
        assert(exclusive(set(sharks.scans.keys()),set(jets.scans.keys())))

        # Check pairs are within one CES of each other in length
        if not flavour.startswith('deck'):
            assert(abs(len(sharks.scans.keys())-len(jets.scans.keys()))<=1),\
           '  jackknife pairs should be within one CES in length of each other'

        # Write splits to file, using title to carry name
        sharks.title = '%s_%s.a' % (patch,splitname)
        jets.title   = '%s_%s.b' % (patch,splitname)

        # Check acceptlist has not been modified
        assert(acceptlist==backup), 'acceptlist has been modified!'

        if writefiles:
            return WriteSplits(sharks,jets,patch,outdir)
        else:
            return sharks,jets

    return

#-----------------------------------------------------------

def FetchLeakages(f='/data4/quiet/calib_data/ingunn/i2qu_160709.dat'):

    """
    Read module leakage values from file
    """

    import scanf

    file=open(f,'r')

    # Working here

    linef= 'i2qu_leak   0   0   -1   -1   %i   %i    %f'
    leakages={}
    for line in file:
        mod,diode,leakages[(mod,diode)]=scanf.sscanf(line,linef)

    file.close()

    return leakages

#-----------------------------------------------------------

def ConvertLeakages(leakages):

    """
    Take a dictionary of diode leakages from FetchLeakages and
    return a dictionary of mean module Q leakage
    """

    meanleakages={}
    meanleak={}

    for (mod,di) in sorted(leakages.keys()):
#        print mod,di,leakages[(mod,di)]
        if 'Q' in Diodes[di]:
#            print mod,di,leakages[(mod,di)]
            if mod not in meanleakages.keys():
                meanleakages[mod]=[]
            meanleakages[mod].append(abs(leakages[(mod,di)]))

    for mod in sorted(meanleakages.keys()):
        meanleak[mod]=sum(meanleakages[mod])/float(len(meanleakages[mod]))

    return meanleak

#-----------------------------------------------------------

def SplitModule(acceptlist,flavour='alternate',outdir='.',patch=None,\
                writefiles=True,verbose=False):

    """
    Jackknife L2 acceptlist by module, and output to acceptlist files
    Inelegant!
    Operates on dictionary of acceptlists keyed by patch
    (outputted by e.g. GenerateCMBRunlists()).
    NB Jackknife halves are not currently forced to be the same length.
    NB These splits are currently Q-band specific
    Optionally write half-acceptlists to file
    """

    import quiet, copy

    assert(writefiles or \
           (len(acceptlist.keys())<=1) and not writefiles), \
           "Can't run silently for more than one patch"

    band='Q'

    # Don't do anything with TTs for now!
    TT=[17,18]
    # Or with the dead module...
    dead=[16]

    patches=[]
    if patch is None: # If patch not specified explicitly
        patches=acceptlist.keys()
    else: # If patch passed explicitly
        patches.append(patch)

    backup=copy.deepcopy(acceptlist)

    for patch in patches:

        sharks = copy.deepcopy(acceptlist[patch]) # Fresh copy
        jets   = copy.deepcopy(acceptlist[patch]) # Fresh copy

        mab1=[0,1,2,3,4,5,6]
        mab2=[10,11,14,15]
        mab3=[7,8,9,12,13]

        # Jackknife by alternate modules
        if flavour=='alternate':
            splitname='modalt'
            odd =[1,3,5,7,9,11,13,15]
            even=[0,2,4,6,8,10,12,14]
            a=odd
            b=even

        # Jackknife by central frequency of module
        if flavour=='frequency':
            splitname='modfreq'
            lower=[15,3,14,6,13,8,4,0]
            upper=[2,5,11,7,12,9,10,1]
            lower.sort()
            upper.sort()
            a=lower
            b=upper

        # Jackknife inner v outer
        if flavour=='radius':
            splitname='modrad'
            inner=[4,5,8,9,10,13,14]
            outer=[0,1,2,3,6,7,11,12,15]
            a=inner
            b=outer

        # Jackknife left v right
        if flavour=='horizontal':
            splitname='modhoriz'
            left=[7,8,9,10,11,12,13,14,15]
            right=[0,1,2,3,4,5,6]
            a=left
            b=right

       # Jackknife in stripes type 1
        if flavour=='stripes1':
            splitname='modstrip1'
            black=[0,2,3,5,7,9,11,13,15]
            white=[1,4,6,8,10,12,14]
            a=black
            b=white

       # Jackknife in stripes type 2
        if flavour=='stripes2':
            splitname='modstrip2'
            black=[0,2,4,6,7,9,11,12,14]
            white=[1,3,5,8,10,13,15]
            a=black
            b=white

       # Jackknife in stripes type 3
        if flavour=='stripes3':
            splitname='modstrip3'
            black=[0,1,2,7,8,9,10,11]
            white=[3,4,5,6,12,13,14,15]
            a=black
            b=white

        # Jackknife top v bottom type 1
        if flavour=='vertical1':
            splitname='modvert1'
            top=[0,1,3,4,7,8,12,13]
            bottom=[2,5,6,9,10,11,14,15]
            a=top
            b=bottom

        # Jackknife top v bottom type 1
        if flavour=='vertical2':
            splitname='modvert2'
            top=[0,3,4,7,8,9,12,13]
            bottom=[1,2,5,6,10,11,14,15]
            a=top
            b=bottom

        # Jackknife by MAB 1v23
        if flavour=='mab1':
            splitname='modmab1'
            mab23=mab2+mab3
            mab23.sort()
            a=mab1
            b=mab23

        # Jackknife by MAB 2v13
        if flavour=='mab2':
            splitname='modmab2'
            mab13=mab1+mab3
            mab13.sort()
            a=mab2
            b=mab13

        # Jackknife by MAB 3v12
        if flavour=='mab3':
            splitname='modmab3'
            mab12=mab1+mab2
            mab12.sort()
            a=mab3
            b=mab12

        # Jackknife by MAB 1v2
        if flavour=='mab12':
            splitname='modmab12'
            a=mab1
            b=mab2

        # Jackknife by MAB 1v3
        if flavour=='mab13':
            splitname='modmab13'
            a=mab1
            b=mab3

        # Jackknife by MAB 2v3
        if flavour=='mab23':
            splitname='modmab23'
            a=mab2
            b=mab3

        # Jackknife by mean module Q leakage
        if flavour=='modleak':
            splitname='modleak'
            modules=GenerateModuleList('Q')
            leakagefile='/data4/quiet/calib_data/ingunn/i2qu_160709.dat'
            leakages=FetchLeakages(leakagefile)
            # Also exclude the two dead diodes
            for (mod,diode) in leakages.keys():
                if mod not in modules or \
                       (mod,diode)==(4,0) or (mod,diode)==(8,2):
                    del leakages[(mod,diode)]
            # Take the mean Q leakage per module, after these exclusions
            meanleak=ConvertLeakages(leakages)
            from operator import itemgetter
            mleaks=sorted(meanleak.items(), key=itemgetter(1))
            a=[]; b=[]
            nmod=len(mleaks)
            # And split on it
            for n,(mod,mleak) in enumerate(mleaks):
                if n <= nmod/2:
                    a.append(mod)
                else:
                    b.append(mod)

        # Exclude the modules from the acceptlists
        if flavour != 'diodes' and flavour != 'leakage' and flavour != 'typea':
            for module in a:
                sharks.excludeModule(module)
            for module in b:
                jets.excludeModule(module)


        # ASIDE: The following jackknives are by diode rather than by module

        # Jackknife by leakage value
        elif flavour == 'leakage':
            splitname='leak'
            modules=GenerateModuleList('Q')
            # Exclude dead #16 and TT:
            for mod in modules[:]:
                if mod in dead or mod in TT: modules.remove(mod)
            leakagefile='/data4/quiet/calib_data/ingunn/i2qu_160709.dat'
            leakages=FetchLeakages(leakagefile)
            # Also exclude the two dead diodes
            for (mod,diode) in leakages.keys():
                if mod not in modules or \
                       (mod,diode)==(4,0) or (mod,diode)==(8,2):
                    del leakages[(mod,diode)]
            from operator import itemgetter
            sleaks=sorted(leakages.items(), key=itemgetter(1))
            ndi=len(sleaks)
            for n,((mod,di),leakage) in enumerate(sleaks):
                if n <= ndi/2:
                    sharks.excludeDiode(mod,di)
                    if verbose:
                        print 'sharks loses %i %i (leakage %4.2f)' \
                              % (mod,di,leakage)
                else:
                    jets.excludeDiode(mod,di)
                    if verbose:
                        print 'jets loses %i %i (leakage %4.2f)' \
                              % (mod,di,leakage)

        elif flavour == 'diodes':
            splitname='modQU'
            modules=GenerateModuleList('Q')
            # [0,1,2,3]=Q1,U1,U2,Q2
            for module in modules[:-3]: # Exclude dead #16 and TT
                sharks.excludeDiode(module,0) # Q1
                jets.excludeDiode(module,1)   # U1
                jets.excludeDiode(module,2)   # U2
                sharks.excludeDiode(module,3) # Q2

        elif flavour == 'typea':
            splitname='modtypa'
            modules=GenerateModuleList('Q')
            # [0,1,2,3]=Q1,U1,U2,Q2
            # email from Colin: 
            # RQ-00 U2, RQ-02 U2, RQ-03 U1, RQ-06 Q1 and Q2,\
            # RQ-08 Q2, RQ-09 U1, RQ-10 Q2, RQ-12 U2, RQ-14 all
            # WORKING HERE
            glitchers=set([(0,U2),(2,U2),(3,U1),(6,Q1),(6,Q2),(8,Q2),(9,U1),\
                       (10,Q2),(12,U2),(14,Q1),(14,U1),(14,U2),(14,Q2)])
            for module in modules[:-3]: # Exclude dead #16 and TT
                for diode in range(quiet.NDIODE_MODULE):
                    if (module,Diodes[diode]) in glitchers:
                        sharks.excludeDiode(module,diode)
                    else:
                        jets.excludeDiode(module,diode)

        # RETURNED: Write out the jackknives

        # It would be good to implement the appropriate assertions here

        # Write splits to file, using title to carry name
        sharks.title = '%s_%s.a' % (patch,splitname)
        jets.title   = '%s_%s.b' % (patch,splitname)

        # Check acceptlist has not been corrupted
        assert(acceptlist==backup), 'acceptlist modified!'

        if writefiles:
            return WriteSplits(sharks,jets,patch,outdir)
        else:
            return sharks,jets

    return

#-----------------------------------------------------------


def makesafe(filename):

    """
    If a file x exists, rename it to x.old to avoid overwriting
    """

    import os

    if os.path.exists(filename):
        os.rename(filename,filename+'.old')

    return filename

#-----------------------------------------------------------


def GetTP(filename="/data4/quiet/level2_tp+b/cmb/patch_2a/cmb_patch_2a_1000_seg001.unf"):

    "Import TP data from L2 file"

    import quiet

    return quiet.l2_read(filename,TP=True)

#-----------------------------------------------------------

def TestTP(tpdata):

    "Print standard deviation of TP data read in from L2 by GetTP"

    for module in tpdata:
        for diode in [0,1,2,3]:
            print module.TP[diode].std()

    return

#-----------------------------------------------------------

def ApplySingleCut(dompi,cut,patch,outdir='.'):

    import os

    """
    Apply a single cut based on the (command-line style) argument supplied
    Write out the .cut file
    """

    # Handle .cut output directory here
    if not outdir: outdir='.'
    #if not outdir.endswith('/'): outdir+='/'

    cutok = True # status - not really used (yet)

    cut = cut.lower() # Enforce lowercase
    cutargs=re.split('\s+',cut)
    print 'Applying cut: '+cut

    if 'pwv' in cutargs:
        CutPWV(cut,outdir)

    if 'tp' in cutargs:
        CutTP(cut,outdir)

    if 'hk' in cutargs:
        CutHK(dompi,cut,outdir)

    if 'tencl' in cutargs:
        pass

    if 'bookends' in cutargs:
        pass

    if 'static' in cutargs:
        if 'time' in cutargs:
            # Write out static cuts
            CutStatic(outdir)
        # elif rather than if needed in next line, or got both!
        elif 'modules' in cutargs or 'diodes' in cutargs or 'dead' in cutargs:
            CutDead(outdir)
        else:
            print 'No cuts made - need to specify type of static cut required'

    if 'test' in cutargs:
        file=os.path.join(outdir,'test.cut')
        testfile=OpenFile(file)
        mjdstart=34236454.1
        mjdstop=13413413.1
        patch = patch
        ces = 432.1
        modulelist='all'
        diodelist='all'
        comment='test'
        WriteCutLine(testfile,mjdstart,mjdstop,ces,patch,modulelist,\
                     diodelist,comment)
        testfile.close()

    return cutok

#-----------------------------------------------------------

def SplitAcceptlist(acceptlist=None,patch=None,flavourdict={},\
                    verbose=False):

    """
    Take a single acceptlist, and apply a set of jackknives to it
    Return pairs of half-acceptlists, keyed by flavour
    """

    import copy

    assert(acceptlist is not None)
    assert(patch is not None)

    # Fetch pointing information if needed
    if 'point' in flavourdict.keys():
        rundict,junkaccdict=GenerateCMBRunlists()

    # Fetch ces_stats if needed
    if 'frequency' in flavourdict.keys():
        ces_stats=ReadCESStats()

    pairs={}

    acceptdict={patch:acceptlist}
    backup=copy.deepcopy(acceptdict)

    for kind in flavourdict.keys():
        flavours=flavourdict[kind]
        if type(flavours)!= list: flavours=[flavours] # Ensure value is a list
        for flavour in flavours:
            print 'Running at %s/%s' % (kind,flavour)
            if kind=='time':
                a,b=SplitTime(acceptdict,flavour=flavour,writefiles=False)
            elif kind=='module':
                a,b=SplitModule(acceptdict,flavour=flavour,patch=patch,\
                                writefiles=False)
            elif kind=='point':
                a,b=SplitPointing(acceptdict,rundict,flavour=flavour,\
                                  writefiles=False)
            elif kind=='weather':
                a,b=SplitWeather(acceptdict,flavour=flavour,\
                                  writefiles=False)
            elif kind=='gain':
                a,b=SplitGain(acceptdict,flavour=flavour,\
                                  writefiles=False)
            elif kind=='frequency':
                a,b=SplitFrequency(acceptdict,ces_stats=ces_stats,\
                              flavour=flavour,writefiles=False)
            if verbose:
                print a.title,b.title
                print len(a.scans.keys()),len(b.scans.keys())
            pairs[flavour]=[a,b]
            assert(backup==acceptdict), '*acceptdict changed by %s %s' \
                                        % (kind,flavour)

    return pairs

#-----------------------------------------------------------

def ApplySingleJackknife(split,outdir=None,acceptlist=None,runlist=None):

    import os

    """
    Apply a single split based on the (command-line style) argument supplied
    Write out the .x files
    """

    # Handle .cut output directory here
    if outdir is None: outdir = '.'

    #splitok = True # status - not really used (yet)

    split = split.lower() # Enforce lowercase
    splitargs=re.split('\s+',split)
    print 'Applying jackknife '+split

    # Enforce jackknife
    # Read in existing cut files
    # Specify how to split data
    # Write out a new cut file

    # Read runlists if not passed as explicit arguments
    # Split into required patches
    if runlist is None or acceptlist is None:
        print 'Generating runlist and acceptlist from master'
        runlist,acceptlist=GenerateCMBRunlists(silent=True)
    else: pass

    if 'time' in splitargs:
        # Split into two by...
        flavours = ['season',       # ...contiguous halves\
                    'interleaved',  # ...interleaved halves\
                    'battenberg',   # ...alternate quarters\
                    'tandem',       # ...end quarters against central pair\
                    'chequerboard'] # ... [NOT IMPLEM] alternate CES pairs

        for flavour in flavours:
            if flavour in splitargs:
                aname,bname=SplitTime(acceptlist,flavour,outdir)

    if 'module' in splitargs:
        # Split into two by...
        flavours = ['alternate',    # ... alternate modules \
                    'radius',       # ... inner v outer \
                    'horizontal',   # ... left v right \
                    'vertical1',    # ... top v bottom (type 1) \
                    'vertical2',    # ... top v bottom (type 2) \
                    'stripes1',     # ... diagonal stripes (type 1) \
                    'stripes2',     # ... diagonal stripes (type 2)\
                    'stripes3',     # ... diagonal stripes (type 3) \
                    'mab1',         # ... MAB 1 v 23 \
                    'mab2',         # ... MAB 2 v 13 \
                    'mab3',         # ... MAB 3 v 12 \
                    'mab12',        # ... MAB 1 v 2 \
                    'mab13',        # ... MAB 1 v 3 \
                    'mab23',        # ... MAB 2 v 3 \
                    'diodes',       # ... Q v U diodes \
                    'frequency',    # ... mean central frequency
                    'typea']        # ... Type A glitchers v not

        for flavour in flavours:
            if flavour in splitargs:
                aname,bname=SplitModule(acceptlist,flavour,outdir)

    if 'pointing' in splitargs:
        # Split into two by...
        flavours = ['elevation',    # ... elevation high v low \
                    'dk1',       # ...  \
                    'dk2',        # ... \
                    'dk3']    # ... 

        for flavour in flavours:
            if flavour in splitargs:
                aname,bname=SplitPointing(acceptlist,runlist,flavour,outdir)

    return aname,bname

#-----------------------------------------------------------

def ApplyMultipleOperations(opersfile,patch,acceptlist=False,runlist=False):

    "Apply multiple (cut/split) operations read from .cmd file"

    import os

    print 'Read command file', opersfile
    operations=ReadList(opersfile)

    # Room for trivial parallelization here:
    cutargs=[]
    splitargs=[]
    for operation in operations[:]:
        operation = operation.lower()
        # Split operation to determine its type
        opargs=re.split('\s+',operation)
        # Remove operation type from arguments
        opertype=opargs.pop(0)
        # And reassemble for passing to Apply function
        operation=' '.join(opargs)
        if opertype == 'cutdir':
            cutdir = opargs[0]
            if not os.path.isabs(cutdir):
                cutdir=os.path.abspath(cutdir)
            print 'Cuts --> dir', cutdir
        elif opertype == 'jackdir':
            jackdir = opargs[0]
            if not os.path.isabs(jackdir):
                jackdir=os.path.abspath(jackdir)
            print 'Jackknives --> dir', jackdir
        elif opertype == 'cut':
            cutargs = opargs
            #print opertype, cutargs
            try: cutdir # declare if does not already exist
            except NameError: cutdir=None
            ApplySingleCut(dompi,operation,options.patch,cutdir)
        elif opertype == 'split':
            splitargs = opargs
            try: jackdir # declare if does not already exist
            except NameError: jackdir=None
            ApplySingleJackknife(operation,jackdir)
        else:
            print 'Ignoring illegal operation', opertype

    return


#-----------------------------------------------------------

def P3t_VtoC(value):

    """
    Convert cold temperatures from V to K [Osamu]
    Make sure PYTHONPATH has /usit/titan/u1/jonathaz/src/L1Read/trunk/python
    Make sure LD_LIBRARY_PATH has /usit/titan/u1/sigurdkn/local/lib
    """

    import conversion.murata

    converter_p3t = conversion.murata.murata('conversion/Murata.txt')

    return converter_p3t.get_T_from_R_interpolation(abs(value*100))

#-----------------------------------------------------------

def P2t_VtoK(value):

    """
    Convert warm temperatures from V to K [Osamu]
    Make sure PYTHONPATH has /usit/titan/u1/jonathaz/src/L1Read/trunk/python
    Make sure LD_LIBRARY_PATH has /usit/titan/u1/sigurdkn/local/lib
    """

    import conversion.curve10

    curve10 = conversion.curve10.curve10('conversion/Curve10.txt')

    return curve10.get_T_from_V(abs(value))

#-----------------------------------------------------------

def PhaseSwitchV2A(input,reverse=False):

    """
    Convert phase-switch volts to amps
    or do the reverse if requested
    """

    # From l1_read.pro
    HK_GAIN=-1.0145
    HK_OFFSET=-0.006

    if reverse:
        return (float(input)-HK_OFFSET)/HK_GAIN
    else:
        return HK_OFFSET+(HK_GAIN*float(input))

#-----------------------------------------------------------

def CheckPhaseSwitches(filelist='/usit/titan/u1/jonathaz/phsw.txt'):

    """
    """

    f=open(filelist,'r')

    threshold=0.5
    extn='HQ_A11PS11I'
    for line in f:
        fi=line.split()[0]
        FlagSingleHKExtension('/data4/quiet/'+fi,extn,threshold)

    return

#-----------------------------------------------------------

def FlagSingleHKExtension(f,extn,threshold):

    """
    For one L1 file, and one extension, flag the file if any of the data
    in that HK extension
    Useful for extn = HQ_A11PS21I | HQ_A12PS21I (threshold=0.5 mA)
    """

    import numpy
    import pyfits

    if extn is None:
        extn='HQ_A11PS21I'

    if not extn.startswith('HQ_'):
        print 'This routine is only for housekeeping data!'
        return None

    #data,hdr=pyfits.getdata(filename,extn=extn,header=True)
    #data=pyfits.getdata(filename,extn,header=False)

    vector=pyfits.getdata(f,extn,header=False).field('Value')
    #print min
    #print max
    print '%s %f %f %f' % (f,PhaseSwitchV2A(vector.max()),PhaseSwitchV2A(vector.min()),threshold)

    # Working here

    return
    #return pyfits.getdata(f,extn,header=False).field('Value').max()>threshold

#-----------------------------------------------------------

def ReadL1HK(filename='/data4/quiet/level1/2008/12/B-120108002935-0.fits'):

    """
    Read housekeeping data from L1 file
    """

    import os, sys, re, pyfits

    # Experiment to use pyfits.getdata to only retrieve selected extns..
    #extensions=range(28,521)
    #for extn in extensions:
    #    extensions[extn].remove()
    #    data,hdr=getdata(filename,extn,header=True)

    hduList = pyfits.open(filename)

    # Fetch list of binary-table extensions
    extensionslist=[]
    for extnum in range(0,len(hduList)):
        extname = hduList[extnum].name
        #print extnum,extname
        extensionslist.append(extname)

    # Delete unwanted data extensions (by list of regular expressions)
    # Delete the data extensions (or how to just read them?)
    removals=['^RQ',         # Remove: # L1 data itself
              '^PRIMARY',              # The (not very useful) FITS header
              '^HQ\_A\d{2}G\d{2}V$',   # GVs
              '^HQ\_A\d{2}D\d{2}I$',   # DIs
              '^HQ\_A\d{2}S\d{2}B$',   # SBs
              '^HQ\_A\d{2}S\d{2}A$',   # SAs
              '^HQ\_M0\dSENSE$',       # MOx senses
              '^HQ\_A\d{2}PS\d{2}I$',  # PSW currents
              '^RUN_MODE',             # ?
              '^TELESCOPE'             # ?
              ]

    for extname in extensionslist[:]:
        for removal in removals:
            if re.search(removal,extname):
                #print removal, extname, re.search(removal,extname)
                extensionslist.remove(extname)

    return hduList,extensionslist

#-----------------------------------------------------------

def ReadL1Runlist(runfile='/data4/quiet/runlist_Q.txt'):

    """
    Read L1 runlist for CMB patches
    """

    import quiet

    # Read L1 runlist
    runlist=quiet.L1Runlist(runfile)

    # Remove any non-CMB patches
    cmblist=runlist.targetTypes['cmb']
    for patch in cmblist.targets.keys():
        if patch not in CMBPatches: cmblist.removeTarget(patch)

    return cmblist

#-----------------------------------------------------------

def FetchMultiTData(parallel=False,data_rootdir='/data4/quiet'):

    """
    Fetch aggregated L1 HK temperature data for many L1 files
    """

    import os, quiet

    nmax=4 # Temporary!

    host = os.getenv('HOSTNAME').split('.')[0]
    if re.search('login-0-\d',host) and parallel:
        from mpi4py import MPI
        print "Enabling MPI... use:"
        print "qlogin --account=quiet --constraint=ib --tasks-per-node=8 --nodes=10"
        world=MPI.COMM_WORLD
        rank=world.rank
        size=world.size
        print 'MPI info:',rank, size
        dompi=True
    else:
        print "MPI disabled..."
        print "Have you thought about using it? I know I do"
        rank=0
        dompi=False

    # Handle incoming file information:
    #if type(filelist) == type(str()): # One file, string
    #    filelist=[filelist]            #          string -> list
    # Wanted to make this a little more all-singing, but no:
    # IN FACT it's no longer relevant - fetch by L1 file NOT L2 CES
    #if type(ceslist[0]) != type(str()): # then assume it's scan,seg
    #    for cespair in ceslist:
    #        ces='.',join(str(s) for s in [cespair]) # -> list of strings

    # This step is slow - only do once

    # Fetch the L1 runlist
    L1run = ReadL1Runlist() # Need to add explicit CMB option

    mapping={}
    filelist=[]
    # For each patch (assumes CMB only has been fetched)
    for patch in sorted(L1run.targets.keys()):
        level1run=L1run.targets[patch]
        # For each scan
        for scan in sorted(level1run.scans.keys()):
            thisscan=level1run.scans[scan]
            mapping[scan]=thisscan.files
            for f in thisscan.files:
                filelist.append(f)

    # Collate HK temperature data
    T2={}
    T3={}
    # Parallelize here..
    if not dompi:
        for f in filelist[0:nmax]: # Temporary
            T2[f],T3[f]=FetchTData(rank,f,data_rootdir)
            #print T2, T3
    elif dompi:
        for f in filelist[0:nmax]:
        #for f in filelist[rank::size]:
            print "Processor %d handling file %s" % (rank,f)
            T2[f],T3[f]=FetchTData(rank,f,data_rootdir)

    return T2,T3

#-----------------------------------------------------------

def FetchTData(myid,L1file='',data_rootdir='/data4/quiet'):

    """
    Fetch aggregated L1 HK temperature data for one L1 file
    """

    import os

    #if not data_rootdir.endswith('/'): data_rootdir+='/'

    filename = os.path.join(data_rootdir,L1file)

    # Read in the data
    hduList,extnlist=ReadL1HK(filename)

    print myid,': Aggregating', filename

    # Aggregate the T data
    T = AggregateHKTemperatureData(hduList)

    #if error: print 'Error in file', filename

    # Split temperatures into cold and warm:
    T2=T['P2T'] # } How to do safely on one line?
    T3=T['P3T'] # }

    return T2,T3

#-----------------------------------------------------------

def AggregateHKTemperatureData(inhduList,
                               avg='median',TempType='both'):

    """
    Combine temperature data for this CES, to allow for cutting
    """

    import pyfits, re, numpy

    #problem=False

    hduList=inhduList[:] # Make a copy so that the input list is not modified

    regexps=[]
    if TempType == 'both':
        regexps.append('^HQ\_P2T\d{1,2}') # P2Ts
        regexps.append('^HQ\_P3T\d{1,2}') # P3Ts        
    elif TempType == 'cold':
        regexps.append('^HQ\_P2T\d{1,2}') # P2Ts
    elif TempType == 'warm':
        regexps.append('^HQ\_P3T\d{1,2}') # P3Ts

    # Make copies of hdulist, one for each regexp/data type
    hduLists={}
    for regexp in regexps:
        #print 'Handling', regexp
        hduLists[regexp]=hduList[:]

        # For each copy, remove irrelevant extension tables
        for extn in hduLists[regexp]:
            tablelist=hduLists[regexp]
        #print extn.name, regexp, re.search(regexp,extn.name)
        #extname = hduList[extnum].name
            if not re.search(regexp,extn.name):
                hduLists[regexp].remove(tablelist[extn.name])
            #print 'removing', extn.name
            #else:
            #    #print 'keeping ', extn.name
            #    pass

    # Now run through each copy and aggregate the temperature data
    T_dict={}
    for regexp in regexps:
        hduList=hduLists[regexp]
        T_vector=[]
        T_med  = -99.0
        T_mean = -99.0
        T_rms  = -99.0 # Don't do anything else with this..
        iter=0
        #for iter,(mjd1,mjd2) in enumerate(hduList[0].data.field('Time')):
        # I want to rewrite this loop so it matches strictly by MJD
        nmin=min(extn.header['NAXIS2'] for extn in hduList[:])
        nmax=max(extn.header['NAXIS2'] for extn in hduList[:])
        if nmin != nmax:
            print 'discrepancy of ',nmin, nmax, ' in ', hduList[0]._file
        nsamples=nmin
        for iter in range(0,nsamples):
            T_i=[]
            for extn in hduList[:]:
                T_i.append(hduList[extn.name].data.field('Value')[iter])
                #except IndexError:
                # *Sometimes* index out of bounds
                # e.g. '/data4/quiet/level1/2008/10/B-102608155848-0.fits'
                # NAXIS2  =                 1578 / number of rows in table
                # BUT >>> e.data.field('Value')[1578]
                # IndexError: index out of bounds
                #    print 'Error:', extn.name, iter
                #    print 'File:', hduList[extn.name]._file
                #    print 'Header:', hduList[extn.name].header
                #    print '      ', len(hduList[0].data.field('Value')[:])

            # Technically we shouldn't do any averages over the MUXed HK
            #     i.e. we are assuming T variations >> 1 second
            T_med = numpy.median(T_i)
            T_mean = numpy.mean(T_i)

            Tc_med=ConvertT(T_med,regexp)
            Tc_mean=ConvertT(T_mean,regexp)

            mjd1,mjd2=hduList[0].data.field('Time')[iter]
            time=ConvertMJD(mjd1,mjd2)

            # Could have generated a dict keyed by mjd time
            T_vector.append([time,Tc_med,Tc_mean])

        key=re.findall('P\dT',regexp)[0]
        T_dict[key]=T_vector

    return T_dict

#-----------------------------------------------------------

def ConvertMJD(mjd1,mjd2):

    """
    Convert mjd1/days + mjd2/msec into a single MJD in days
    """

    return mjd1 + mjd2/(24.0*60.0*60.0*1000.0)

#-----------------------------------------------------------

def ConvertT(T,regexp):

    """
    Handle the temperature-conversion curves
    """

    import re

    if re.search('P2T',regexp):
        return P2t_VtoK(T)
    elif re.search('P3T',regexp):
        return P3t_VtoC(T)
    else:
        return None

    return

#-----------------------------------------------------------

##############################
#  Main program starts here  #
##############################

if __name__ == '__main__':

    usage = 'modes are: \n ./secateur.py -h \n'
    parser = OptionParser(usage=usage)

    # For help (no arguments):
    #if len(parser.parse_args()) < 3 :
    #    __name__='secateur'
    #    help(__name__)
    #    sys.exit(0)

    # Specify options

    parser.add_option("-c", "--cut", dest="cut",
                      help="make a cut",type="string",metavar='[pwv | static[dead|time] | hk]')

    parser.add_option('-s', '--split', dest='split',
                      help='make a jackknife:',type='string',metavar='[ time [season|interleaved|battenberg|tandem] | module[alternate|radius|horizontal|vertical1|vertical2|stripes1|stripes2|stripes3|mab1|mab2|mab3|mab12|mab13|mab23|diodes|frequency] ] NB NEED ONE PAIR OF QUOTES AROUND SENTENCE')

    # Generate CMB L2 runlists from master runlist, without editing
    
    parser.add_option('-g', '--gen-runlists', action='store_true',\
                      dest='genrun',default=False, help='generate L2 runlists')

    # Use MPI where appropriate
    parser.add_option('-m', '--mpi', action='store_true',\
                      dest='mpi',default=False, help='use MPI where appropriate')

    # Specify output dir for .cut files
    parser.add_option("-o", "--output-cut", dest="cutdir",
                      help="specify output path for .cut files", type="string",metavar='.')

    # Specify output dir for .a .b jackknife files
    parser.add_option("-j", "--output-jack", dest="jackdir",
                      help="specify output path for .a .b jackknife files", type="string",metavar='.')

    # Specify APEX data file
    parser.add_option("-a", "--fetch-pwv", dest="apexfile",
                      help="download APEX PWV data to file", type="string",metavar='apex.dat')

    # Specify QUIET band
    parser.add_option('-b', '--band', dest='band',
                      help='specify QUIET band', type='string', default='Q',metavar='[Q|W]')

    # Specify CMB patch
    parser.add_option("-p", "--patch", dest="patch",
                      help="specify QUIET CMB patch", type="string", default="patch_2a",metavar='[patch_2a|patch_4a|patch_6a|patch_7a|patch_7b]')

    # Specify operations (cuts/splits) commands file
    parser.add_option("-f", "--opers-file", dest="operscmd",
                      help="specify cut/split commands", type="string",metavar='cuts.cmd')

    # Toggling this flag may help with diagnostics
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      default=False, help="make lots of noise")

    # Reassert defaults
    parser.set_defaults(verbose=False,band="Q",patch="patch_2a")
    #parser.set_defaults(verbose=False,apexfile="apex.dat",band="Q")

    # Parse command line
    (options, args) = parser.parse_args()

    # Select patches
    if options.patch:
        if options.patch not in CMBPatches:
            print 'Patch is not a CMB science patch'
            print 'Exiting'
            sys.exit(1)
        else:
            pass
            #patches=option.patch

    if options.cutdir:
        #if not options.cutdir.endswith('/'):
        #    options.cutdir+='/'
        if not os.path.isabs(options.cutdir):
            options.cutdir=os.path.abspath(options.cutdir)
            
    else:
        options.cutdir=os.path.abspath('.')

    if options.mpi:
        dompi=True
    else:
        dompi=False

    if options.jackdir:
        #if not options.jackdir.endswith('/'):
        #    options.jackdir+='/'
        if not os.path.isabs(options.jackdir):
            options.jackdir=os.path.abspath(options.jackdir)
    else:
        options.jackdir=os.path.abspath('.')

    # Generate CMB L2 runlists from master runlist, without editing
    if options.genrun:
        GenerateCMBRunlists(silent=False)

    # Fetch APEX data if required
    if options.apexfile:
        if options.band:
            print 'Fetching', options.band, '-band data to', options.apexfile
            WriteAPEXData(options.band,options.apexfile)

     # Apply single cut from command line
    if options.cut:
        ApplySingleCut(dompi,options.cut,options.patch,options.cutdir)

    # Apply single jackknife from command line
    if options.split:
        ApplySingleJackknife(options.split,options.jackdir)

    # Apply multiple cut/split operations from file
    if options.operscmd:
        ApplyMultipleOperations(options.operscmd,options.patch)

    # Backstop
    sys.exit(0)


############################
#  Main program ends here  #
############################

#-----------------------------------------------------------
