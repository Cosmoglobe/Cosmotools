#!/usr/bin/python

# TO ADD:
# Check for pre-existing files to avoid spurious appending

import sys,os,string,time
import numpy,scipy,pylab,matplotlib
import quiet,secateur
from numpy import *
from matplotlib import *
from scipy import stats
from quiet import N_MODULE_MAX_Q,N_MODULE_MAX_W


#LEVEL2_DIR_DEFAULT='/data4/quiet/level2/'
LEVEL2_DIR='/usit/titan/u1/sigurdkn/wruns/level2_ref/'
LEVEL2_DIR_DEFAULT=LEVEL2_DIR

#-------------------------------------------------------------------

def constructOutputString(run,seg,time,nmodules,array,index1,index2):

    """
    From the rmsarr, and various other variables,
    construct and return the output string to be written to file
    """

    s2=' %d %d %e ' % (run, seg, time)
    for module in range(nmodules):
        s2 += '%e ' % array[index1,index2,module]
    s2+='\n'

    return s2

#-------------------------------------------------------------------

if __name__ == '__main__':

    from optparse import OptionParser

    # Set up options
    usage = 'This runs through all runs in runlist_patch.txt\n and calculates, for each module, the avg(rms) and rms(rms) for binned data\n and outputs 12 files, 6 for each Q diode (binned, and binned detrended):\n      average(rms): Qx_Mean_STD_6a_binsz_.txt\n      rms(rms): Qx_STD_STD_6a_binsz_.txt\n      average(average): Qx_Mean_Mean_binsz_.txt\n using the N-1 convention for the normalization of the st. deviation\n'

    parser = OptionParser(usage=usage)            

    parser.add_option('-p', '--patch', dest='patch',\
                      help='specify QUIET CMB patch', type='string',\
                      metavar='(2a|4a|6a|7b)')

    parser.add_option('-r', '--runfile', dest='runfile',\
                      help='specify L2 runfile', type='string',\
                      metavar='runlist_l2.txt')

    parser.add_option('-t', '--timescale', dest='timescale',\
                      help='specify time scale for binning',\
                      type='int',metavar='(10|30)')

    parser.add_option('-l', '--l2dir', dest='l2dir',\
                      help='specify L2 dir', type='string',\
                      metavar='/path/to/level2/')

    parser.add_option('-b', '--band', dest='band',\
                      help='specify band', type='string',\
                      metavar='(Q|W)')

    parser.set_defaults(patch=None,runfile=None,timescale=None,\
                        band=None,l2dir=LEVEL2_DIR_DEFAULT)

    # Parse command line
    (options, args) = parser.parse_args()

    patch=options.patch
    runfile=options.runfile
    timescale=options.timescale
    l2dir=options.l2dir
    band=options.band

    #set the name of the level 2 dir (so we can work out where the files are)
    #level2 = "/data4/quiet/level2/"
    level2 = l2dir

    bins = str(timescale)
    binsize = timescale
    patch = patch
    runfile = runfile

    if band == 'Q':
        nmodules=N_MODULE_MAX_Q
    elif band == 'W':
        nmodules=N_MODULE_MAX_W
    else:
        nmodules=None

    print '\nBand %s has %i modules\n' % (band,nmodules)


#if len(sys.argv) !=4:
#      print >> sys.stderr,' %s <bin size> <patch (2a,4a,6a,7b)> ' \
#            % os.path.basename(sys.argv[0])
#      print 'This runs through all runs in runlist_patch.txt'
#      print 'and calculates, for each module, the avg(rms) and rms(rms) for binned data'
#      print 'and outputs 12 files, 6 for each Q diode (binned, and binned detrended): '
#      print '      average(rms): Qx_Mean_STD_6a_binsz_.txt'
#      print '      rms(rms): Qx_STD_STD_6a_binsz_.txt'
#      print '      average(average): Qx_Mean_Mean_binsz_.txt'
#      print ' using the N-1 convention for the normalization of the st. deviation'
#      sys.exit(2)

#bins = sys.argv[1]
#binsize = int(bins)
#patch = sys.argv[2]
#runfile = sys.argv[3]

#load the runlist into a structure
#theseruns = string.join(['/usit/titan/u1/newburgh/repository/quiet_svn/oslo/src/python/weather_cutting/runlist_',patch,'.txt'],'')
#theseruns='runlist_6103.txt'
#theseruns='runlist_l2.txt'
#runlist=quiet.L2Runlist(theseruns)

    print 'Reading %s...' % runfile
    run=quiet.L2Runlist(runfile)

    # Remove from runlist all patches except those specified
    print 'Cleaning runlist...'
    runlist=secateur.CleanL2Runlist(run,'patch_%s'%patch)

    # get the full table of segments from the runlist
    segments = runlist.segments()
    nsegments=len(segments)
    print '\nThere are %i CESs to consider for patch %s\n' %(nsegments,patch)
    t0=time.time()

    # Mapping for output filenames
    outputname={1:('Q1_Mean_Mean',''),\
                2:('Q2_Mean_Mean',''),\
                3:('Q1_Mean_STD',''),\
                4:('Q2_Mean_STD',''),\
                5:('Q1_STD_STD',''),\
                6:('Q2_STD_STD',''),\
                7:('Q1_Mean_Mean','_detr'),\
                8:('Q2_Mean_Mean','_detr'),\
                9:('Q1_Mean_STD','_detr'),\
                10:('Q2_Mean_STD','_detr'),\
                11:('Q1_STD_STD','_detr'),\
                12:('Q2_STD_STD','_detr'),\
                13:('U1_Mean_STD',''),\
                14:('U2_Mean_STD',''),\
                15:('U1_STD_STD',''),\
                16:('U2_STD_STD','')}

    # loop through the segment IDs
    for n,(run_id,segment_id) in enumerate(sorted(segments)):
        t1=time.time()

        # from the table, get the segment info (including the filename)
        segment_info = segments[(run_id,segment_id)]

        # get the full path to the file with this segment data in
        filename = segment_info.filename(level2)
        print filename, bins
        # load the data from the file
        data = quiet.l2_read(filename)

        # down-sample the data into 1 sec bins, use only Q diodes and make:
        #    ds_array[(meanQ1, meanQ2, stdQ1, stdQ2),module,bin]
        #    where std has N-1 normalization, mean is not weighted by std
        l=len(data[0].tod.Q1)
        deltat=50
        ds=l/deltat
        if ds > 1:
            ds_array=zeros([8,nmodules,ds-1])
            for i in range(1,ds):
                for module in range(nmodules):
                    ds_array[0,module,i-1] = mean(data[module].tod.Q1[deltat*(i-1):deltat*i])
                    ds_array[1,module,i-1] = mean(data[module].tod.Q2[deltat*(i-1):deltat*i])
                    ds_array[2,module,i-1] = numpy.std(data[module].tod.Q1[deltat*(i-1):deltat*i],ddof=1)
                    ds_array[3,module,i-1] = numpy.std(data[module].tod.Q2[deltat*(i-1):deltat*i],ddof=1)
                    ds_array[4,module,i-1] = mean(data[module].tod.U1[deltat*(i-1):deltat*i])
                    ds_array[5,module,i-1] = mean(data[module].tod.U2[deltat*(i-1):deltat*i])
                    ds_array[6,module,i-1] = numpy.std(data[module].tod.U1[deltat*(i-1):deltat*i],ddof=1)
                    ds_array[7,module,i-1] = numpy.std(data[module].tod.U2[deltat*(i-1):deltat*i],ddof=1)  

    # bin the data into whatever bin size you input originally, make:
    #   binnddat[(mean(meanQ1),mean(mean(Q2),std(mean(Q1)),std(mean(Q2))),repeat detrend),module,bin]
        deltat = binsize # number of samples to rms over, 1 sammple = 1 sec
        binned=int(ds/deltat)
        if binned > 1:
            binneddat=zeros([12,nmodules,binned-1])
            for i in range(1,binned):
                for module in range(nmodules):
                    binneddat[0,module,i-1] = average(ds_array[0,module,deltat*(i-1):deltat*i]) # mean of demod Q1
                    binneddat[1,module,i-1] = average(ds_array[1,module,deltat*(i-1):deltat*i]) # mean of demod Q2
                    binneddat[2,module,i-1] = numpy.std(ds_array[0,module,deltat*(i-1):deltat*i],ddof=1) # std of demod Q1
                    binneddat[3,module,i-1] = numpy.std(ds_array[1,module,deltat*(i-1):deltat*i],ddof=1) # std of demod Q2
                    binneddat[4,module,i-1] = mean(matplotlib.mlab.detrend(ds_array[0,module,deltat*(i-1):deltat*i])) #mean of detrend Q1
                    binneddat[5,module,i-1] = mean(matplotlib.mlab.detrend(ds_array[1,module,deltat*(i-1):deltat*i])) #mean of detrend Q2
                    binneddat[6,module,i-1] = numpy.std(matplotlib.mlab.detrend_linear(ds_array[0,module,deltat*(i-1):deltat*i]),ddof=1) # std of det Q1
                    binneddat[7,module,i-1] = numpy.std(matplotlib.mlab.detrend_linear(ds_array[1,module,deltat*(i-1):deltat*i]),ddof=1) # std of det Q2
                    binneddat[8,module,i-1] = average(ds_array[4,module,deltat*(i-1):deltat*i]) # mean of demod Q1
                    binneddat[9,module,i-1] = average(ds_array[5,module,deltat*(i-1):deltat*i]) # mean of demod Q2
                    binneddat[10,module,i-1] = numpy.std(ds_array[4,module,deltat*(i-1):deltat*i],ddof=1) # std of demod Q1
                    binneddat[11,module,i-1] = numpy.std(ds_array[5,module,deltat*(i-1):deltat*i],ddof=1) # std of demod Q2
          

            # create the data output
            # make rmsarray[(mean(mean),mean(std),std(std),repeat for detrend),module]
            # use N-1 normalization for std dev.
            rmsarr=zeros([3,6,nmodules])
            for i in range(nmodules):
                rmsarr[0,0,i] = average(binneddat[0,i,0:binned-1])
                rmsarr[0,1,i] = average(binneddat[1,i,0:binned-1])
                rmsarr[0,2,i] = mean(binneddat[2,i,0:binned-1])
                rmsarr[0,3,i] = mean(binneddat[3,i,0:binned-1])
                rmsarr[0,4,i] = numpy.std(binneddat[2,i,0:binned-1],ddof=1)
                rmsarr[0,5,i] = numpy.std(binneddat[3,i,0:binned-1],ddof=1)
                rmsarr[1,0,i] = average(binneddat[4,i,0:binned-1])
                rmsarr[1,1,i] = average(binneddat[5,i,0:binned-1])
                rmsarr[1,2,i] = mean(binneddat[6,i,0:binned-1])
                rmsarr[1,3,i] = mean(binneddat[7,i,0:binned-1])
                rmsarr[1,4,i] = numpy.std(binneddat[6,i,0:binned-1],ddof=1)
                rmsarr[1,5,i] = numpy.std(binneddat[7,i,0:binned-1],ddof=1)
                rmsarr[2,0,i] = average(binneddat[8,i,0:binned-1])
                rmsarr[2,1,i] = average(binneddat[9,i,0:binned-1])
                rmsarr[2,2,i] = mean(binneddat[10,i,0:binned-1])
                rmsarr[2,3,i] = mean(binneddat[11,i,0:binned-1])
                rmsarr[2,4,i] = numpy.std(binneddat[10,i,0:binned-1],ddof=1)
                rmsarr[2,5,i] = numpy.std(binneddat[11,i,0:binned-1],ddof=1)

            # Open the output files
            outputfile={}

            for index in sorted(outputname.keys()):
                f='%s_%s%s_%s_allmods.txt' \
                   % (outputname[index][0],patch,outputname[index][1],bins)
                outputfile[index]=open(f,'a+')

            #outputfile[1] = open(string.join(['Q1_Mean_Mean',patch,bins,'allmods.txt'],'_'),"a+") # output mean(mean(ds_Q1))
            #outputfile[2] = open(string.join(['Q2_Mean_Mean',patch,bins,'allmods.txt'],'_'),"a+") # output mean(mean(ds_Q2))
            #outputfile[3] = open(string.join(['Q1_Mean_STD',patch,bins,'allmods.txt'],'_'),"a+") # output mean(mean(ds_Q1))
            #outputfile[4] = open(string.join(['Q2_Mean_STD',patch,bins,'allmods.txt'],'_'),"a+") # output mean(mean(ds_Q2))
            #outputfile[5] = open(string.join(['Q1_STD_STD',patch,bins,'allmods.txt'],'_'),"a+") # output mean(mean(dt_Q1))
            #outputfile[6] = open(string.join(['Q2_STD_STD',patch,bins,'allmods.txt'],'_'),"a+") # output mean(mean(dt_Q2))
            #outputfile[7] = open(string.join(['Q1_Mean_Mean',patch,'detr',bins,'allmods.txt'],'_'),"a+") # output mean(rms(dt_Q1))
            #outputfile[8] = open(string.join(['Q2_Mean_Mean',patch,'detr',bins,'allmods.txt'],'_'),"a+") # output mean(rms(dt_Q2))
            #outputfile[9] = open(string.join(['Q1_Mean_STD',patch,'detr',bins,'allmods.txt'],'_'),"a+") # output mean(rms(dt_Q1))
            #outputfile[10] = open(string.join(['Q2_Mean_STD',patch,'detr',bins,'allmods.txt'],'_'),"a+") # output mean(rms(dt_Q2))
            #outputfile[11] = open(string.join(['Q1_STD_STD',patch,'detr',bins,'allmods.txt'],'_'),"a+") # output mean(rms(dt_Q1))
            #outputfile[12] = open(string.join(['Q2_STD_STD',patch,'detr',bins,'allmods.txt'],'_'),"a+") # output mean(rms(dt_Q2))
            #outputfile[13] = open(string.join(['U1_Mean_STD',patch,bins,'allmods.txt'],'_'),"a+") # output mean(mean(ds_Q1))
            #outputfile[14] = open(string.join(['U2_Mean_STD',patch,bins,'allmods.txt'],'_'),"a+") # output mean(mean(ds_Q2))
            #outputfile[15] = open(string.join(['U1_STD_STD',patch,bins,'allmods.txt'],'_'),"a+") # output mean(mean(dt_Q1))
            #outputfile[16] = open(string.join(['U2_STD_STD',patch,bins,'allmods.txt'],'_'),"a+") # output mean(mean(dt_Q2))

#            outputfile2.write(' %d %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n' % (run_id, segment_id, data[0].time[0], rmsarr[0,1,0],rmsarr[0,1,1],rmsarr[0,1,2],rmsarr[0,1,3],rmsarr[0,1,4],rmsarr[0,1,5],rmsarr[0,1,6],rmsarr[0,1,7],rmsarr[0,1,8],rmsarr[0,1,9],rmsarr[0,1,10],rmsarr[0,1,11],rmsarr[0,1,12],rmsarr[0,1,13],rmsarr[0,1,14],rmsarr[0,1,15], rmsarr[0,1,16],rmsarr[0,1,17], rmsarr[0,1,18]))

            # Write to output files
            dtime=data[0].time[0]
            outputfile[1].write(constructOutputString(run_id,segment_id,dtime,nmodules,rmsarr,0,0))
            outputfile[2].write(constructOutputString(run_id,segment_id,dtime,nmodules,rmsarr,0,1))
            outputfile[3].write(constructOutputString(run_id,segment_id,dtime,nmodules,rmsarr,0,2))
            outputfile[4].write(constructOutputString(run_id,segment_id,dtime,nmodules,rmsarr,0,3))
            outputfile[5].write(constructOutputString(run_id,segment_id,dtime,nmodules,rmsarr,0,4))
            outputfile[6].write(constructOutputString(run_id,segment_id,dtime,nmodules,rmsarr,0,5))
            outputfile[7].write(constructOutputString(run_id,segment_id,dtime,nmodules,rmsarr,1,0))
            outputfile[8].write(constructOutputString(run_id,segment_id,dtime,nmodules,rmsarr,1,1))
            outputfile[9].write(constructOutputString(run_id,segment_id,dtime,nmodules,rmsarr,1,2))
            outputfile[10].write(constructOutputString(run_id,segment_id,dtime,nmodules,rmsarr,1,3))
            outputfile[11].write(constructOutputString(run_id,segment_id,dtime,nmodules,rmsarr,1,4))
            outputfile[12].write(constructOutputString(run_id,segment_id,dtime,nmodules,rmsarr,1,5))
            # 2,0 not written }- why?
            # 2,1 not written }
            outputfile[13].write(constructOutputString(run_id,segment_id,dtime,nmodules,rmsarr,2,2))
            outputfile[14].write(constructOutputString(run_id,segment_id,dtime,nmodules,rmsarr,2,3))
            outputfile[15].write(constructOutputString(run_id,segment_id,dtime,nmodules,rmsarr,2,4))
            outputfile[16].write(constructOutputString(run_id,segment_id,dtime,nmodules,rmsarr,2,5))


        # Report the timing info
        t2=time.time()
        dt=t2-t1
        rt=t2-t0
        pt=rt*nsegments/(60.0*(n+1))
        lt=rt*(nsegments-(n+1))/(60.0*(n+1))
        print 'That took %6.2f sec (%6.2f sec total, %6.2f min predict, %6.2f min left)' % (dt,rt,pt,lt)

    # Close file handles and exit
    for index in sorted(outputfile.keys()):
        outputfile[index].close()

    print '\nFinished!\n'

    sys.exit(0)

