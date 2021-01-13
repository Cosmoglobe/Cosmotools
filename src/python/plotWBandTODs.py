#!/usr/bin/env python

"""
Read an L2/L3 file for a given CES and module - mjd, tod, tp, housekeeping
Convert to numpy arrays
Plot to file using pyplot
"""

import matplotlib
matplotlib.use("Agg")

#-------------------------------------------------------------------

# Define ordering of diodes
Diodes=['Q1','U1','U2','Q2']
DiodesLookup={'Q1':0,'U1':1,'U2':2,'Q2':3}
Q1='Q1'; U1='U1'; U2='U2'; Q2='Q2'

#-------------------------------------------------------------------

biasMapping = {
    'GA3' :0,\
    'DB2I':1,\
    'DB2V':2,\
    'GB3' :3,\
    'GB2' :4,\
    #'06A3GA2' :5,\
    'DA1I':5,\
    'DA1V':6,\
    'GA1' :7,\
    'DB1I':8,\
    'DB1V':9,\
    'GB1' :10,\
    'DA2I':11,\
    'DA2V':12}

#-------------------------------------------------------------------

def fetchModuleMapping(f=\
   '/usit/titan/u1/jonathaz/quiet_svn/oslo/src/f90/include/quiet_l1_defs.f90'):

    """
    Extract channel numbering from quiet_l2_defs.f90
    """

    import re

    oldlines=open(f,'r').readlines()

    lines=[]

    for line in oldlines:
        if 'M' in line and '=' in line:
            lines.append(line)

    #print lines

    all=' '.join([e.rstrip('\n') for e in lines])

    #print all

    p = re.compile('M\d\dA\d[A-Z]+\d\s=\s+\d+')
    q = re.compile('M\d\dA\d[A-Z]+\d[A-Z]\s=\s+\d+')

    single=p.findall(all)+q.findall(all)

    #print single

    mapping={}

    for el in single:
        name,mapping[name]=el.split(' = ')[0],int(el.split(' = ')[1])

    #print mapping

    return mapping

#-------------------------------------------------------------------

def constructMapping(bias_map,module_map,mod_num,field_name):

    """
    Given a module number RWnn, and a variable name (e.g. DA2I),
    return the index of that bias variable
    The mapping dictionary is extracted from quiet_l1_defs.f90
    """

    mod_name=module_map[mod_num] # Fetch the obscure name of that module

    req_name=mod_name+field_name # Construct the full field name

    index=bias_map[req_name] # Find the index for that full field name

    return index

#-------------------------------------------------------------------

def readModuleMapping(\
    f='/projects/quiet/jonathaz/wband/cuts/hk/module-mapping.txt'):

    """
    Read module mapping from file into a dictionary
    """

    data=open(f,'r')

    mapping={}
    for line in data:
        try:
            module,mapping[int(module)]=line.split()
        except:
            pass

    data.close()

    return mapping

#-------------------------------------------------------------------

def moduleMapping(mapping,module):

    """
    Retrieve the string corresponding to the requested module number (0-)
    Read from file each time (pretty fast though...)
    """

    return mapping[module]

#-------------------------------------------------------------------

def mod2dinum(mod):

    """
    Convert a module number (0-) to a range of diode numbers (0-)
    """

    down=4*mod
    up=down+3

    return down,up

#-------------------------------------------------------------------

def plotTOD(todd,module,ces,accept,otodd,ochannel,\
            dohist,doall,target,downsample):

    """
    Take todd object and for a given module, plot DD,TP to file
    2 up: 2 x [DD,TP] per page
    Replace TP with a single bias variable if required
    #1 up, one per diode (to allow room for housekeeping)
    """

    import numpy
    import matplotlib.pyplot as plt
    from scipy import stats

    V2mV=1.0e3
    dpi=100 # => 800  x 600  default
    dpi=150 # => 1200 x 900  pixels
    dpi=200 # => 1600 x 1200 pixels

    print numpy.shape(todd),numpy.shape(otodd)

    #mjdstart=min(todd[0,:])
    #mjdstop=max(todd[0,:])

    mjdstart=todd[0,0]
    mjdstop=todd[0,-1]

    deltat=86400.0*(mjdstop-mjdstart)
    print 'TOD (%f-%f) lasts %i seconds' % (mjdstart,mjdstop,deltat)

    zp=numpy.floor(mjdstart*10.0)/10.0

    deltat_start=int(86400.0*(mjdstart-zp))
    deltat_stop=int(86400.0*(mjdstop-zp))

    zp_sec=numpy.floor(deltat_start/1000.0)*1000.0
#    print zp_sec

    xlabel = 'Time after MJD %f plus %i sec / seconds' % (zp,zp_sec)
    span = '%i sec -> %i sec' % (deltat_start-zp_sec,deltat_stop-zp_sec)

    print xlabel
    print span
    #####

    mjd=todd[0,:]
    if otodd is not None: omjd=otodd[0,:]

    r=readModuleMapping()  # Fetch the module mapping
    f=fetchModuleMapping() # Fetch the bias mapping

    if doall:
        npanels=3
    else:
        npanels=2

    down,up=mod2dinum(module)

    for di in range(down,up+1):
        if otodd is None:
            flav=''
        elif doall:
            flav='_all'
        else:
            flav='_%s' % ochannel
        fign='ML%i_RW%02i%s_%i_st%i.png'\
              %(ces,module,flav,di%4,accept[(module,di%4)])
        title='*%i* - %s ML%i RW%02i di%i %03i %s - %f-%f' \
               % (accept[(module,di%4)],target[-2:],ces,module,di%4,di,\
                  Diodes[di%4],mjdstart,mjdstop)
        plt.figure(1,dpi=dpi)

        plt.subplot(npanels,1,1)
        style=','
        data=todd[di+1,:] # +1 for the MJD
        plt.plot(mjd,data*V2mV,style)
        plt.xlabel('MJD')
        plt.ylabel('DD/mV')
        plt.title(title)

        plt.subplot(npanels,1,2)
        if otodd==None or doall:
            data=todd[di+1+364,:] # 1 for the MJD and 364 for the DDs
            plt.plot(mjd,data*V2mV,style)
            plt.xlabel('MJD')
            if downsample is None:
                plt.ylabel('TP/mV')
            else:
                rms=data[:int(downsample/2.0)].std()
                plt.ylabel('DD_%i (%4.2e)/mV' % (downsample,rms))
        else:
            style=','
            oindex=constructMapping(f,r,module,ochannel)
            print 'index %i for module %i channel %s' \
                  % (oindex,module,ochannel)
            index=oindex

            odata=otodd[index,:]

            skew=stats.skew(odata)
            kurtosis=stats.kurtosis(odata)
            print 'ML%i %i x %i %s ind %i skew %5.2f kurt %5.2f st%i' \
                  % (ces,module,di%4,ochannel,index,skew,kurtosis,accept[(module,di%4)])

            plt.plot(omjd,1000.0*odata/65536.0,style)
            plt.xlabel('MJD (skew %5.2f kurt %5.2f)' % (skew,kurtosis))
            plt.ylabel('%s/mV (index %i)' % (ochannel,oindex))


        if doall:
            plt.subplot(npanels,1,3)
            style=','
            oindex=constructMapping(f,r,module,ochannel)
            print 'index %i for module %i channel %s' \
                  % (oindex,module,ochannel)
            index=oindex

            odata=otodd[index,:]

            skew=stats.skew(odata)
            kurtosis=stats.kurtosis(odata)
            print 'ML%i %i x %i %s ind %i skew %5.2f kurt %5.2f st%i' \
                  % (ces,module,di%4,ochannel,index,skew,kurtosis,accept[(module,di%4)])

            plt.plot(omjd,1000.0*odata/65536.0,style)
            plt.xlabel('MJD (skew %5.2f kurt %5.2f)' % (skew,kurtosis))
            plt.ylabel('%s/mV (index %i)' % (ochannel,oindex))

        plt.savefig(fign,dpi=(dpi))
        print 'Look in %s' % fign
        plt.close('all')

    # Relying on the last cycle of the last loop..
    if otodd is not None and dohist:
        figm='ML%i_RW%02i%s_hist.png' % (ces,module,flav)
        title='ML%i - RW%02i (skew %5.2f kurt %5.2f)' \
               % (ces,module,skew,kurtosis)
        plt.hist(1000.0*odata/65536.0,bins=20,histtype='bar')
        plt.xlabel('%s/mV (index %i)' % (ochannel,oindex))
        plt.ylabel('Counts')
        plt.title(title)
        plt.savefig(figm,dpi=(dpi))
        print 'Look in %s' % figm
        plt.clf()
        plt.close('all')
    ### NEW ###

    return

#-------------------------------------------------------------------

if __name__ == '__main__':

    from optparse import OptionParser

    print
    # Set up options
    usage = "./plotWBandTODs.py -h\n"
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose", action="store_false",\
                      dest="verbose",help="make lots of noise")

    parser.add_option("-p", "--plot", action="store_true",\
                      dest="plot",help="plot TODs")

    parser.add_option("-b", "--bias", action="store_true",\
                      dest="bias",help="Plot bias rather than TP")

    parser.add_option("-a", "--all", action="store_true",\
                      dest="all",help="Plot bias AND TP")

    parser.add_option("-s", "--stats", action="store_true",\
                      dest="dostats",help="Compute stats")

    parser.add_option("-t", "--dohist", action="store_true",\
                      dest="dohist",help="plot histograms")

    parser.add_option("-i", "--infile", dest="statsfile",\
           help="specify input ces_stats file", type="string",\
           metavar='ces_statistics.txt')

    parser.add_option("-f", "--field", dest="ochannel",\
           help="specify bias field", type="string",\
           metavar='DA2I')

    parser.add_option("-c", "--ces", dest="ces",\
                      help="specify ces", type="int",\
                      metavar='412')

    parser.add_option("-m", "--module", dest="module",\
                      help="specify module", type="int",\
                      metavar='4')

    parser.add_option("-d", "--downsample", dest="downsample",\
                      help="downsample TOD by this factor", type="int",\
                      metavar='500')

    # Assert defaults
    parser.set_defaults(statsfile=None,ces=None,module=None,\
                        verbose=False,plot=False,dostats=False,\
                        bias=False,ochannel='',dohist=False,\
                        all=False,downsample=None)

    # Parse command line
    (options, args) = parser.parse_args()

    verbose=options.verbose
    statsfile=options.statsfile
    ces=options.ces
    module=options.module
    plot=options.plot
    dostats=options.dostats
    bias=options.bias
    ochannel=options.ochannel
    dohist=options.dohist
    all=options.all
    downsample=options.downsample

    assert(ces is not None), 'ces required!'

    import sys,os,string
    import numpy,scanf
    from scipy import stats
    import h5py,quietW

    #data=quietW.l3_read(ces=ces)  #can also read by filename

    #print 'Available fields are - %s' % ' '.join(data.fields)
    #prints:
    # ['alpha', 'apex', 'bias', 'coord_sys', 'corr', 'cryo', 'decimation', 'diode_stats', 'encl', 'fknee', 'gain', 'nside', 'orig_point', 'peri', 'pixels', 'point', 'samprate', 'scanfreq', 'sigma0', 'stats', 'time', 'time_gain', 'tod', 'tp']

    #mjd=data.time #mjd is numpy array
    #print data.tp.shape
    #prints (364,33799)
    # so it is organized by absolute diode number (0-3 = module 0,
    # 4-7=module 1, etc .) etc.

    #if you don't want all the fields it is faster to read in like this:
    fields=['time','tod','tp']
    data=quietW.l2_read(ces=ces,fields=fields)
    #print data.time.shape, data.tod.shape, data.tp.shape

    mjd=numpy.array(data.time)
    tod=numpy.array(data.tod)
    tp=numpy.array(data.tp)

    print 'm',numpy.shape(tod)

    # Downsample the data if required
    if downsample is not None:
        rate=int(downsample/2.0) # window -> +/- window/2
        mjd_temp=mjd.copy()
        tod_temp=tod.copy()
        tp_temp=tp.copy()

        nsamp=numpy.shape(mjd_temp)[0]
        ndiodes=numpy.shape(tod)[0]

        # NB This skips the bookends
        print 'Diode',
        for di in range(mod2dinum(module)[0],mod2dinum(module)[1]+1):
            print di,
            for s in range(rate,nsamp-rate):
                tod_temp[di,s]=tod[di,s-rate:s+rate].mean()
                tp_temp[di,s]=tp[di,s-rate:s+rate].mean()
            #print sample,mjd_temp[sample],tod_temp[:,sample],tp_temp[:,sample]
        print

        # Now arrange for
        #   (i)  tod = original tod (DD)
        #   (ii) tp  = downsampled tod
        # (or whatever you want.. - set it up here)
        mjd=mjd_temp
        tod=tod
        tp=tod_temp

    todd=numpy.vstack([mjd,tod,tp])
    #print numpy.shape(todd)
    # todd now a numpy array thus:
    # mjd tod(di0-di4) tp(di0-di4)
    #print len(todd[1][:])

    if bias or all:
        ofields=['bias']
        odata=quietW.l3_read(ces=ces,fields=ofields)
        print odata.bias.time.shape, odata.bias.value.shape

        omjd=numpy.array(odata.bias.time)
        obias=numpy.array(odata.bias.value)
        otodd=numpy.vstack([omjd,obias])
        print numpy.shape(otodd)
    else:
        otodd=None

    #print module,' '.join([str(x) for x in mod2dinum(module)])

#    x=numpy.ones(5,9)

    #for t in range(5):
    #    print mjd[t],tod[t,4*module:4*module+4],tp[t,4*module:4*module+4]
    #    x[t,:]=mjd,tod[t,4*module:4*module+4],tp[t,4*module:4*module+4]


    # Gather status
    if True:
        if statsfile is None:
            statsfile='./ces_statistics.txt'
        indata=open(statsfile,'r')
        accept={}
        print 'Reading %s (for the statuses)...' % statsfile
        for line in indata:
            if line.startswith('#'): continue
            try:
#                print line
#                sys.exit(0)
                cesid,obj,mod,di,status,r,sigma0,alpha,fknee,\
                  chisq_025,sscan,chisq_002,CMB,chisq_1025,\
                  TODchisq,TODabs,chi_az,typeb,wthr10,wther30,jump,d1,d2,d3\
                  =line.split()
            except:
                print 'Skipping %s' % line
                pass

            if cesid != str(ces):
                continue
            else:
                accept[(int(mod),int(di))]=int(status)

        indata.close()

    if False:
        dumpTOD(todd,module,ces)


    # Gather parent patch
    # Slow because must be reread every time (as with ces_stats)
    target=quietW.Runlist().targetForCES(ces)
    print 'target is %s' % target

    if plot:
        plotTOD(todd,module,ces,accept,\
                otodd,ochannel,dohist,all,target,downsample)
#        sys.exit(0)

    #print '%i x %i skew %5.2f kurt %5.2f' \
    #      % (module,di%4,stats.skew(todd[di+1+364,:]),stats.kurtosis(todd[di+1+364,:]))

    if dostats:
        statsf='stats_ML%i.txt' % ces
        stat=open(statsf,'w')
        stat.write('# ces mod di status skew kurt\n')
        print 'Computing TP stats for all modules, diodes for this CES'

        dis=range(364)

        for di in dis:
            line = '%i %i %i %i %5.2f %5.2f\n' \
                   % (ces,di/4,di%4,accept[(di/4,di%4)],\
                      # Add 365 to shift to the TPs
                      stats.skew(todd[di+1+364,:]),\
                      stats.kurtosis(todd[di+1+364,:]))
            if verbose: print line
            stat.write(line)

        stat.close()


    print '\nFinished!\n'
    sys.exit(0)

