#!/usr/bin/env python

"""
Read a W-band ces_stats.txt file
Select the object and cut type
Histogram each column for that object
Write the histograms to file
"""

#-------------------------------------------------------------------

    # Columns are

    #    CES id   Object   Mod  Di  Status   R     Sigma0     Alpha
    # Fknee    0-12.5Hz   Scan    0-0.2Hz    CMB     10-12.5Hz TOD_chisq
    # TOD_abs  chi_az   type-B   Wthr-10s  Wthr-30s

mapping = { 1:'CES id',\
            2:'Object',\
            3:'Mod',\
            4:'Di',\
            5:'Status',\
            6:'R',\
            7:'Sigma0',\
            8:'Alpha',\
            9:'Fknee',\
            10:'0-12.5Hz',\
            11:'Scan',\
            12:'0-0.2Hz',\
            13:'CMB',\
            14:'10-12.5Hz',\
            15:'TOD_chisq',\
            16:'TOD_abs',\
            17:'chi_az',\
            18:'type-B',\
            19:'Wthr-10s',\
            20:'Wthr-30s'}

panel = {1:[0,1,2,3],2:[4,5,6,7],3:[8,9,10,11],4:[12,13,14,14]}

#-------------------------------------------------------------------

def compare(x,y,op):
    from operator import eq,gt,lt
    return locals()[op](x,y)

def compare2(x,y,op):
    import operator
    return getattr(operator,op)(x,y)

#-------------------------------------------------------------------

def hist_slice(data,target,stat,op,nbins,col,percentile,dinum):

    """
    Take the data object
    Cut it in the requested way
    Histogram it for the required column
    Currently one target at a time
    """

#    data[(obj,id,mod,di,status)]=[r,sigma0,alpha,fknee,\
#                                  chisq_025,sscan,chisq_002,CMB,\
#                                  chisq_1025,TODchisq,TODabs,chi_az,\
#                                  typeb,wthr10,wther30]

    print 'Applying cuts (to %s)...' % stat
    dat=[]
    reportfreq=100000
    i=1
    for (obj,id,mod,di,status) in sorted(data.keys()):
        if obj == target and compare(status,stat,op):
            dat.append(data[(obj,id,mod,di,status)][col])
        if (i % reportfreq == 0): print 'Reached %i' % i
        i+=1

    if dat == [] :
        print 'Abort - dat is empty for target %s!' % target
        sys.exit(0)

    plot_solo(dat,nbins,col,target,op,stat,percentile,dinum)

    return

#-------------------------------------------------------------------

def hist_slices(data,target,stat,op,nbins):

    """
    Take the data object
    Cut it in the requested way
    Histogram it for the required column
    Currently one target at a time
    """

#    data[(obj,id,mod,di,status)]=[r,sigma0,alpha,fknee,\
#                                  chisq_025,sscan,chisq_002,CMB,\
#                                  chisq_1025,TODchisq,TODabs,chi_az,\
#                                  typeb,wthr10,wther30]

    dat={}

    for pkey in sorted(panel.keys()):
        for r in range(0,15):
            reportfreq=100000
            i=1
            dat[pkey]=[]
            for (obj,id,mod,di,status) in sorted(data.keys()):
                if obj == target and compare(status,stat,op):
                    dat[pkey].append(data[(obj,id,mod,di,status)][r])
                if (i % reportfreq == 0): print 'Reached %i' % i
                i+=1
        plot_multi(dat,nbins,target,op,stat,pkey)

    return

#-------------------------------------------------------------------

def plot_multi(dat,nbins,target,op,stat,pkey):

    """
    Plot four histograms to a single png
    """

    import matplotlib.pyplot as plt
    import numpy,sys,os
    from scipy import stats

    percentile=95
    #score=stats.scoreatpercentile(numpy.array(dat),percentile)
    #print '%s: Percentile %i -> upper lim %8.3f' \
    #      % (mapping[col+6], percentile, score)

    fign='hist_%s_st%s%i_panel%i_%i_%i.png' \
          % (target,op,stat,pkey,percentile,nbins)

    print panel[pkey]

    # Start the plot
    plt.figure(1)

    for index in panel[pkey]:
        print index
        plt.subplot(2,2,1+index%4)
        score=stats.scoreatpercentile(numpy.array(dat[index]),percentile)
        plt.yscale('log')
        plt.hist(dat[index],bins=nbins,\
                 range=(min(dat[index]),score),histtype='bar')
        plt.xlabel('%s' % mapping[index+6],fontsize=7)

    plt.title('%s - status %s%i (%i percent, %i bins)' \
              % (target,op,stat,percentile,nbins),fontsize=9)

    plt.savefig(fign)
    plt.close('all')

    return

#-------------------------------------------------------------------

def plot_solo(dat,nbins,col,target,op,stat,percentile,dinum):

    """
    Plot a single histogram to a single png
    """

    import matplotlib.pyplot as plt
    import numpy
    from scipy import stats

    score=stats.scoreatpercentile(numpy.array(dat),percentile)
    print '%s: Percentile %i -> upper lim %8.3f' \
          % (mapping[col+6], percentile, score)

    if dinum is not None:
        print 'dinum %i' % dinum
        ddi= 'dinum%i_' % dinum
    else:
        ddi=''

    fign='hist_%s_%sstatus_%s%i_%s_%i_%i.png' \
          % (target,ddi,op,stat,mapping[col+6],percentile,nbins)
    fign=os.path.join(object,fign)

    if dinum is not None:
        ddi= 'dinum%s' % dinum
    else:
        ddi='diall'

    plt.yscale('log')
    plt.hist(dat,bins=nbins,range=(min(dat),score),histtype='bar')
    plt.title('%s - %s status %s%i (%i percent, %i bins)' \
              % (target,ddi,op,stat,percentile,nbins))
    plt.xlabel('%s' % mapping[col+6])
    plt.savefig(fign)
    plt.close('all')

    return

#-------------------------------------------------------------------

if __name__ == '__main__':

    import sys,os,string
    from optparse import OptionParser
    #import scanf

    print
    # Set up options
    usage = "./histogram_ces_stats -h.py\n"
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose", action="store_true",\
                      dest="verbose",help="make lots of noise")

    parser.add_option("-i", "--infile", dest="infile",\
           help="specify input file", type="string",\
           metavar='ces_stats.txt')

    parser.add_option("-o", "--object", dest="object",\
           help="specify object", type="string",\
           metavar='(jupiter|patch_2a|all) etc')

    parser.add_option("-s", "--status", dest="status",\
           help="specify status cut", type="string",\
           metavar='(eq0|gt3)')

    parser.add_option("-n", "--nbins", dest="nbins",\
           help="specify nbins", type="int",\
           metavar='(10,20,50,100)')

    parser.add_option("-p", "--percentile", dest="percentile",\
           help="specify percentile", type="int",\
           metavar='(68,95,99)')

    parser.add_option("-d", "--dinum", dest="dinum",\
           help="specify dinum for histograms", type="int",\
           metavar='00')

    parser.add_option("-m", "--multiplot", action="store_true",\
                      dest="multiplot",help="attempt 4-up plots")

    # Assert defaults
    parser.set_defaults(\
        infile='input.txt',object='patch_2a',n=20,s='eq0',\
        verbose=False,multiplot=False,percentile=95,dinum=None)
    # Parse command line
    (options, args) = parser.parse_args()

    verbose=options.verbose
    infile=options.infile
    object=options.object
    status=options.status
    nbins=options.nbins
    percentile=options.percentile
    multiplot=options.multiplot
    dinum=options.dinum

    assert(infile is not None), 'ces_stats required!'

    indata=open(infile,'r')
    data={}
    reportfreq=100000
    i=1
    print 'Reading %s...' % infile

    for line in indata:
        if line.startswith('#'): continue
        try:
            id,obj,mod,di,stat,r,sigma0,alpha,fknee,\
              chisq_025,sscan,chisq_002,CMB,chisq_1025,\
              TODchisq,TODabs,chi_az,typeb,wthr10,wther30\
              =line.split()
        except:
            print 'Skipping %s' % line
            pass

        #print id,obj,mod,di,status,r,sigma0,alpha,fknee
        #print chisq_025,sscan,chisq_002,CMB,chisq_1025
        #print TODchisq,TODabs,chi_az,typeb,wthr10,wther30

        # key by obj (can always choose all)
        # And below it by status (for now, anyway)

        if object == obj or object == 'all': # 'all' not working yet..
            #if module is not None and module == int(mod):
            if dinum is None or \
                   (dinum is not None and dinum==int(di)+4*int(mod)):
                data[(str(obj),int(id),int(mod),int(di),int(stat))]=\
                   [float(r),float(sigma0),float(alpha),float(fknee),\
                    float(chisq_025),float(sscan),float(chisq_002),float(CMB),\
                    float(chisq_1025),float(TODchisq),float(TODabs),float(chi_az),\
                    float(typeb),float(wthr10),float(wther30)]

        if (i % reportfreq == 0): print 'Reached line %i' % i
        i+=1

    indata.close()
    print 'Read %s OK for object %s' % (infile,object)

    # Decode status word
    op,status=(status[0:2],int(status[2:]))
    print 'Status %s %i' % (op,status)
    #status=3; op = 'gt'

    if not multiplot:
        if not os.path.exists(object): os.mkdir(object)
        cols=range(0,15)
        for col in cols:
            hist_slice(data,object,status,op,nbins,col,percentile,dinum)
    else:
        hist_slices(data,object,status,op,nbins)

    sys.exit(0)

