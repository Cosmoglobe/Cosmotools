#!/usr/bin/env python

"""
"""

#-------------------------------------------------------------------

def exprToList(expr,dir='.'):

    """
    Glob a string into a list of files, usually for that dir
    """

    import glob,os

    ladder=glob.glob(os.path.join(dir,expr))

    return ladder

#-------------------------------------------------------------------

def plotCls(wrapped_data,directory,page):

    """
    """

    import matplotlib.pyplot as plt
    import pylab

    jks=wrapped_data.keys()

    modes=['EE','BB']

    plt.figure(1,figsize=(20,60),dpi=300)
    title='%s' % directory
    style='.'

    params={'xtick.labelsize': 12,'axes.labelsize':12}
    plt.rcParams.update(params)

    nwide=len(modes)
    njks=len(wrapped_data.keys())
    npanels=nwide*njks
    #print nwide,njks,npanels

    print 'Plotting...'

    limits=3#sigma
    lines =2#sigma

    #    print numpy.mean(wrapped_data[:][:][:,1])

    ipanel=1
    for ijk,jk in enumerate(jks):
        data=wrapped_data[jk]
        for imode,mode in enumerate(modes):
            uplim = max(data[mode][:,1])
            lolim = -min(data[mode][:,1])
            maxlim=max(uplim,lolim)

            if maxlim>3.5:
                plt.subplot(njks,nwide,ipanel,axisbg='red')
            elif maxlim>3:
                plt.subplot(njks,nwide,ipanel,axisbg='yellow')
            else:
                plt.subplot(njks,nwide,ipanel)

            #print 'Plotting panel %i (%i %i) ...' % (ipanel,ijk,imode)
            plt.axhline(y=0,color='black')
            plt.axhline(y=lines,color='black',linestyle=':')
            plt.axhline(y=-lines,color='black',linestyle=':')
            avg1 = numpy.mean(data['EE'][0:5,1])   
            avg2 = numpy.mean(data['BB'][0:5,1])
            avg = (avg1+avg2)/2.0
            avg1 = numpy.mean(data['EE'][10:22,1])   
            avg2 = numpy.mean(data['BB'][10:22,1])
            avg_n = (avg1+avg2)/2.0

            if mode=='EE':
                print avg, avg_n

            plt.axhline(y=avg,color='green',linestyle=':')
            plt.plot(data[mode][:,0],data[mode][:,1],style,color='blue')
            plt.xlim(0,1200)
            if uplim > limits:
                uplim=uplim+0.5

            if lolim > limits:
                lolim=lolim+0.5

            plt.ylim(-max(limits,lolim), max(limits,uplim))
            plt.xlabel('ell')
            if ipanel%nwide == 1:
                plt.ylabel('%s' % jk)
            if ipanel > (npanels-nwide):
                plt.xlabel('%s' % mode)
            else:
                plt.xlabel('')
            if ipanel == 2:
                plt.title(title)
            if ipanel < (npanels-nwide):
                params={'xtick.labelsize':12,'xtick.color':'black'}
                plt.rcParams.update(params)
            else:
                params={'xtick.labelsize':12,'xtick.color':'black'}
                plt.rcParams.update(params)
            ipanel+=1

    outf='knives.pdf'
    print 'Writing to file...'
#    plt.savefig(outf, facecolor=plt.facecolor, edgecolor='none')
    plt.savefig(outf)
    
    plt.close()
    print 'Look in %s' % outf

    return 0

#-------------------------------------------------------------------

def plotChis(data_all_modes,directory,ptes):

    """
    """

    import matplotlib.pyplot as plt
    import numpy

    #print wrapped_data.keys()

    plt.figure(1,figsize=(20,24),dpi=300)
    plt.suptitle(directory,fontsize='x-large')

    jks=data_all_modes.keys()
    njks=len(jks)
    nwide=2
    nbins=30
    limits=3#sigma
    #print njks,nwide

    if not ptes:
        range=(-limits,limits)
    else:
        range=(0.0,1.0)

    for ijk,jk in enumerate(jks):
        plt.subplot((njks+1)/2,nwide,ijk+1)
        data=data_all_modes[jk]
        plt.hist(data,bins=nbins,range=range,\
                 histtype='step',color='black',normed=True)
        if not ptes: plt.yscale('log')
        plt.ylabel('%s' % jk)
        if ijk > (njks-2*nwide):
            plt.xlabel('sigma')


    if not ptes:
        outf='chi_hist.pdf'
    else:
        outf='pte_hist.pdf'
    plt.savefig(outf)
    print 'Look in %s' % outf
    plt.close()

    return 0

#-------------------------------------------------------------------


def plotChi2(data_all_modes,directory,ptes):

    """
    """

    import matplotlib.pyplot as plt
    import numpy

    #print wrapped_data.keys()

    plt.figure(1,figsize=(20,48),dpi=300)
    plt.suptitle(directory,fontsize='x-large')

    jks=data_all_modes.keys()
    njks=len(jks)
    nwide=2
    nbins=40
    #limits=3#sigma
    #print njks,nwide

    if not ptes:
        range=(0.0,8.0)
    else:
        range=(0.0,1.0)

    for ijk,jk in enumerate(jks):
        plt.subplot((njks+1)/2,nwide,ijk+1)
        data=data_all_modes[jk]
        plt.hist(data,bins=nbins,range=range,\
                 histtype='step',color='black',normed=True)
        #if not ptes: plt.yscale('log')
        plt.ylabel('%s' % jk)
        if ijk > (njks-2*nwide):
            plt.xlabel('sigma')


    if not ptes:
        outf='chi2_hist.pdf'
    else:
        outf='pte2_hist.pdf'
    plt.savefig(outf)
    print 'Look in %s' % outf
    plt.close()

    return 0

#-------------------------------------------------------------------

def stripJKName(f):

    """
    Retrieve jackknife name from filename
    e.g. '/path/to/patch_6a/cls/norm_cl_patch_6a_typeb.dat' -> 'typeb'
    """

    import os

    f=os.path.basename(f)
    f=f.split('_')[-1]
    f=f.split('.')[0]

    return f

#-------------------------------------------------------------------

def fetchData(list_of_files,sims_only):

    """
    """

    import numpy,sys

    wrapped_data={}
    #print files_chi
    for f in list_of_files:
        #print f
        data=numpy.genfromtxt(f)
        # Remove data entries if required
        if sims_only:
            data=data[69:]

        numbins=numpy.shape(data)[0]
        fkey=stripJKName(f)
        wrapped_data[fkey]=data

        #print f,fkey,numbins
        #print fkey, numpy.shape(data_chi)

        # Cut out the modes
        # Not really necessary?
        #data={}
        #data['EE']=data_chi[:numbins/3]
        #data['EB']=data_chi[numbins/3:2*numbins/3]
        #data['BB']=data_chi[2*numbins/3:]

        # Take the chis and histogram them with data overplotted on sims

    return wrapped_data

#-------------------------------------------------------------------

if __name__ == '__main__':

    import os,sys,numpy

    try:
        indir=sys.argv[1]
    except IndexError:
        print '\nUsage: plotWBandPowSpecs.py [directory name]\n'
        sys.exit(0)

    indir=os.path.abspath(indir)
    print indir


    files_cl    =exprToList(expr='norm_cl_*.dat',dir=os.path.join(indir,'cls'))
    files_pte   =exprToList(expr='pte_*.dat',dir=os.path.join(indir,'cls'))
    files_chi   =exprToList(expr='chi_*.dat',dir=os.path.join(indir,'cls'))
    files_chisq =exprToList(expr='chisq_*.dat',dir=os.path.join(indir,'cls'))
    try:
        # If the run didn't complete cleanly..
        files_summ  =exprToList(\
        expr='jk_summary.dat',dir=os.path.join(indir,'cls'))[0]
    except:
        pass

    #print files_cl
    #print files_pte
    #print files_chi
    #print files_chisq
    #print files_summ

    # Plot the individual (data) power spectra

    # This is a little different to those below since it's split by mode
    wrapped_cls={}
    for f in files_cl:
        data_cl=numpy.genfromtxt(f)
        numbins=numpy.shape(data_cl)[0]
        fkey=stripJKName(f)
        #print f,fkey,numbins

        # Cut out the modes
        data={}
        data['EE']=data_cl[:numbins/3]
        data['EB']=data_cl[numbins/3:2*numbins/3]
        data['BB']=data_cl[2*numbins/3:]
        # BEWARE the ordering of these - at least they are tagged
        #for mode in ['EE','EB','BB']:
        #    print 'Fetched %s %s' % (fkey,mode), #,data[mode]
        wrapped_cls[fkey]=data

    directory=os.path.dirname(indir)

    # The aim here is to pass in dictionary slices for each page
    #for page,jk in enumerate(sorted(wrapped_data.keys())):
    #    status=plotCls(wrapped_data,directory,page)

    # Not used - would have been for paging the output
    # The idea was to slice the dictionary by passing plots in in, say, threes
    page=0
    status=plotCls(wrapped_cls,directory,page)

    # Now histogram the chis
    wrapped_chi=fetchData(files_chi,sims_only=True)
    ## TO BE DONE: Strip out data values and plot these separately over
    status=plotChis(wrapped_chi,directory,ptes=False)

    # Now we need a master histogram

    # Histogram the PTEs
    #wrapped_pte=fetchData(files_pte,sims_only=False)
    ## I think the next function call is wrong! Hist of PTEs..?
    ## ->> Omit it
    #status=plotChis(wrapped_pte,directory,ptes=True)

    wrapped_chisq=fetchData(files_chisq,sims_only=True)
    status=plotChi2(wrapped_chisq,directory,ptes=False)
    # I think the next line also does not make sense, so omit it for now
    #status=plotChi2(wrapped_pte,directory,ptes=True)

    sys.exit(0)


