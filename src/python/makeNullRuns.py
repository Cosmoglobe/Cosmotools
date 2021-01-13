#!/usr/bin/python

##!/usr/bin/python2.4

"""
Collection of utilities for making null tests, including directory structure \n
$Author: laura $
$Date: 2009-11-19 22:27:42 +0100 (Thu, 19 Nov 2009) $
Usage:
      developed on python 2.4-2.5
      ./makeNullRun.py cuttype jack-knife \n
Footnotes: \n
0. The final output will assume you want to remove anything in dead/static cuts
1. I assume this is going into /data4/quiet/nulltests

"""

import sys,platform,os,time,re
from numpy import *
import datetime, string

# -------------------------------------------------------------

def getdirname(maindir):

    """
    get date
    check if that date exists already
    if it does, increment the number from the last run
    """

    os.chdir(maindir)
    now=datetime.datetime.now()
    base = 'pcl-%04d%02d%02d' %(now.year,now.month,now.day)
    yn = 0
    i = 1
    print base, maindir
    while yn == 0:
       dirtomake = '%s-%03d' % (base,i)

       print dirtomake
       if not os.path.isdir(dirtomake):
           yn=1
       else: i=i+1

    return os.path.join(maindir,dirtomake)

# -------------------------------------------------------------

def makedirstructure(dirname):

    """
    make a directory structure:
          cuts
          params
          maps
          powspec
    """

    from secateur import FixPerms

    cuts = os.path.join(dirname,'cuts')
    params = os.path.join(dirname,'params')
    maps = os.path.join(dirname,'maps')
    powspec = os.path.join(dirname,'powspec')
    offsetdir = os.path.join(maps,'saved_offsets')

    # fixperms=False
    fixperms=True
    for directory in [dirname,cuts,params,maps,powspec,offsetdir]:
        os.mkdir(directory)
        if fixperms:
            FixPerms(directory)

    print 'Made Directory: ', dirname, ' and subdirectories: cuts, params, maps, powspec, maps/saved_offsets'

#-----------------------------------------------------------

##############################
#  Main program starts here  #
##############################

if __name__ == '__main__':

    if len(sys.argv) !=3:
          print >> sys.stderr, '%s <cutlist> <jacktype>' % sys.argv[0]
          print 'cutlist = file with all cuts, YOU MUST INPUT FULL PATH'
          print ' (will improve this to take multiple cuts, standardize path)'
          print 'jacktype = type of jack-knife, current options are: '
          print '      season, interleaved, battenberg, tandem'
          print 'if jacktype = all, do all jackknives'
          print 'to add: optional directory input (currently uses /data4/quiet/nulltests'
          sys.exit(2)

    import re
    import quiet, secateur
    import makeRun
    import makeJackknives_091130

    cutlist = sys.argv[1]
    jackknife = sys.argv[2]

    trunkdir = '/data4/quiet/nulltests/pcl'
    directory = getdirname(trunkdir)
    makedirstructure(directory)
    
    os.chdir(string.join([directory,'cuts'],'/'))
    # secateur.CMBPatches()

    patch='patch_6a'
    dompi=False
    #cutdir='.'
    cutdir=string.join([directory,'cuts'],'/')
    jack_outdir=string.join([directory,'cuts'],'/')
    for cut in ['static dead','static time']:
        # produces static.cut, dead.cut in cuts directory
        secateur.ApplySingleCut(dompi,cut,patch,cutdir)

    #jack_outdir='.'
    #cutfiles=['dead.cut','static.cut']

    os.chdir('/'.join([directory, 'cuts']))

    flavours={}
    if jackknife == 'all':
        flavours['time'] = ['season','battenberg','interleaved','tandem']
        flavours['module'] = ['alternate','radius','horizontal','vertical1',\
                              'vertical2','stripes1','stripes2','stripes3',\
                              'mab1','mab2','mab3','mab12','mab13','mab23',\
                              'diodes','frequency','elevation']
        flavours['pointing'] = ['elevation','dk1','dk2','dk3']
    else:
        key,flame=re.split(' ',jackknife)
        flavours[key]=[flame]

    # Cut the cake
    for kind in sorted(flavours.keys()):
        for flavour in flavours[kind]:
            jackknife = kind +' '+ flavour
            cutfiles=[string.join([directory,'/cuts/dead.cut'],''),string.join([directory,'/cuts/static.cut'],'')]
            halves = \
                   makeJackknives_091130.makeJackknife(cutlist,\
                   jack_outdir,cutdir,cutfiles,patch,jackknife)
            for split in halves:
                print split, directory
                makeRun.makeRun(split,directory)

    sys.exit(0)


############################
#  Main program ends here  #
############################


