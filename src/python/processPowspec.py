#!/usr/bin/python

"""
For a cls.fits file outputted by commander,
manipulate data
do KS
produce P
do plots
"""

#-------------------------------------------------------------------
# These two intrinsics are not defined until python 2.5

def all(iterable):
    for element in iterable:
        if not element:
            return False
    return True

def any(iterable):
    for element in iterable:
        if element:
            return True
    return False

# http://stackoverflow.com/questions/2892931/longest-common-substring-from-more-than-two-strings-python

def long_substr(data):
    substr = ''
    if len(data) > 1 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0])-i+1):
                if j > len(substr) and all(data[0][i:i+j] in x for x in data):
                    substr = data[0][i:i+j]
    return substr

def is_substr(find, data):
    if len(data) < 1 and len(find) < 1:
        return False
    for i in range(len(data)):
        if find not in data[i]:
            return False
        return True

#-------------------------------------------------------------------

#import exceptions
#class NumChainsError(exceptions.Exception):
#    def __init__(self):
#        return
#    def __str__(self):
#        print "","A nchains mismatch occured!"

#-------------------------------------------------------------------

def exprToList(expr,dir='.'):

    """
    Glob a string into a list of files, usually for that dir
    """

    import glob,os

    ladder=glob.glob(os.path.join(dir,expr))

    return ladder

#-------------------------------------------------------------------

def readChisqs(f):

    """
    Read a list of chisqs from an ascii file into a numpy array
    """

    import numpy

    data=sorted(numpy.genfromtxt(f))

    return data

#-------------------------------------------------------------------

def readChisToChisqs(f,both=False):

    """
    Read a list of chisqs from an ascii file into a numpy array
    """

    import numpy


    data=numpy.array(sorted(numpy.genfromtxt(f)))
    data_out=numpy.square(data)
    if both: return data_out,numpy.array(data)

    #length=numpy.shape(data)[0]
    #data_out=numpy.zeros(numpy.shape(data)[0])
    #for my_value in range(length):
    #    data_out[my_value] = data[my_value] * data[my_value]

    return data_out

#-------------------------------------------------------------------

def gaussianToPTE(array,mean=0.0,sigma=1.0):

    """
    For an array of chis, calculate their PTEs
    """

    from scipy import stats

    ptes=stats.norm.cdf(array,loc=mean,scale=sigma)

    return ptes

#-------------------------------------------------------------------

def ncx2ToPTE(array,ndof=1,mean=0.0):

    """
    For an array of chisqs, calculate their PTEs
    """

    from scipy import stats
    import numpy

    n=numpy.size(array)

    df=ndof; nc=n*mean*mean
    ptes=stats.ncx2.cdf(array,df,nc)

    return ptes

#-------------------------------------------------------------------

def chisqToPTE(array,ndof=1):

    """
    For an array of chisqs, calculate their PTEs
    """

    from scipy import stats

    ptes=stats.chi2.cdf(array,ndof)

    return ptes

#-------------------------------------------------------------------

def chisq_uniform(data,bins):

    """
    Do a chisq test against a uniform distribution
    Assume the errors are sqrt(number of points in each bin)
    """

    import numpy

    n=numpy.sum(data)
    nbins=len(bins)-1
    nb=float(n)/float(nbins) # Assumes uniform binning
    #print n,nb

    chisq=0.0
    for d in data:
        e=numpy.sqrt(d)
        m=float(nb)
        chi=(float(d)-e)/m
        #print '%f %f %f %f' % (d,e,m,chi)
        chisq+=chi*chi

    ndof=nbins-1
    chisq_red=chisq/float(ndof)
    pte=chisqToPTE(chisq,ndof=ndof)

    return chisq,ndof,chisq_red,pte

#-------------------------------------------------------------------

def reportStats(array,tag,file):

    """
    Do KS and chisq to uniform
    Print to stdout
    """

    import numpy

    patch,mode,realn=extractInfo(file)

    D,p=doKSTest(array)
    c,b=numpy.histogram(array,bins=10,range=(0,1.0))
    cq,ndof,rcq,pp=chisq_uniform(c,b)
    print 'stat %s %s %s %f %i %f %f %f %f %s' \
          % (patch,mode,realn,cq,ndof,rcq,pp,D,p,tag)
    
    return

#-------------------------------------------------------------------

def printPTEs(chi_dict,ptes_dict,ptes_shift_dict,\
              chisq_dict,ptes2_dict,ptes2_shift_dict):

    """
    """

    import numpy

    print '# val patch mode realn chi pte pte_shift chisq pte2 pte2_shift'
    print '# stat patch mode realn chisq ndof red_chisq PTE D_KS P_KS'
    for f in sorted(ptes_dict.keys()):
        patch,mode,realn=extractInfo(f)

        chi=chi_dict[f].tolist()
        ptes=ptes_dict[f].tolist()
        reportStats(ptes,'ptes',f)

        ptes_shift=ptes_shift_dict[f].tolist()
        reportStats(ptes_shift,'ptes_shift',f)

        chisq=chisq_dict[f].tolist()
        ptes2=ptes2_dict[f].tolist()
        reportStats(ptes2,'ptes2',f)

        ptes2_shift=ptes2_shift_dict[f].tolist()
        reportStats(ptes2_shift,'ptes2_shift',f)
#        print '# %s %s %s %f %f %f %f' % \
#              (patch,mode,realn,p,p_shift,p2,p2_shift)
        for n,z in enumerate(chi):
            print 'val %s %s %s %f %f %f %f %f %f' \
                  % (patch,mode,realn,chi[n],ptes[n],ptes_shift[n],\
                     chisq[n],ptes2[n],ptes2_shift[n])

    return

#-------------------------------------------------------------------

def extractInfo(chif):

    """
    Split a chi file name into patch, mode and whether data or sim
    """

    import os

    # Remove extension
    stem = '.'.join(os.path.basename(chif).split('.')[:-1])

    for m in ['EE','BB','EB']:
        if m in stem:
            mode=m
            break
        else:
            mode='all'

    for r in ['data','sim']:
        if r in stem:
            realisation=r
            break
        else:
            realisation='all'

    for p in ['2a','4a','6a','7b']:
        if p in stem:
            patch = p
            break
        else:
            patch = 'all'

    return patch,mode,realisation

#-------------------------------------------------------------------

def processPTEs(dir='.',expr='*_data.dat'):

    """
    """

    import numpy

    files=exprToList(expr=expr,dir=dir)

    data2={}; data={}; pte2={}; pte={}; pte_shift={}; pte2_shift={}

    ndof=1; shift=0.19; sigma=1.0
    for f in files:
        data2[f],data[f]=readChisToChisqs(f,both=True)
        pte2[f]=chisqToPTE(data2[f],ndof=ndof)
        pte[f]=gaussianToPTE(data[f],mean=0.0,sigma=sigma)
        pte_shift[f]=gaussianToPTE(data[f],mean=shift,sigma=sigma)
        pte2_shift[f]=ncx2ToPTE(data2[f],ndof=ndof,mean=shift)

    return data,pte,pte_shift,data2,pte2,pte2_shift

#-------------------------------------------------------------------

def getPTEs(f,g):

    """
    Read two lists of chi-squareds, one of the data (f) and the other of the 
    simulations (g)
    Algorithmically, I think the operation of getPTEs is, or should be,
    the same as that of fetchP
    """
    
    import numpy,math
    
    #print 'f file ', f
    #print 'g file ', g
    our_ptes=numpy.zeros(numpy.shape(f))
    
    #print our_ptes
    #print 'Shape f ', numpy.shape(f), ' Shape g ', numpy.shape(g)

    a = numpy.shape(f)[0]
    b = numpy.shape(g)[0]
    #print a, b
    
    for my_value in range(a):
        delta=10
        pte = 0.0
        for my_second_value in range(b):
           # pte = 0.0
            if (math.fabs(f[my_value] - g[my_second_value]) < delta):
                pte = 0.0
                delta=math.fabs(f[my_value] - g[my_second_value])
                our_ptes[my_value] = float(my_second_value)/float(b)

    return our_ptes

#-------------------------------------------------------------------

def applyPixelWindow(cls,window):

    """
    Apply binned pixel window to binned power spectrum
    """

    import numpy

    assert(numpy.shape(cls)==numpy.shape(window)),\
       'Shape mismatch for cls and attempted window'

    cls_windowed=cls*window*window

    return cls_windowed

#-------------------------------------------------------------------

def readPixelWindow(nside=128):

    """
    """

    import numpy,pyfits

    # Read bins
    b='bins_d50_n%i.dat' % nside
    bins=readBins(b)

    # Read pixel window (T,P)
    f='pixel_window_n%04i.fits' % nside
    d=pyfits.open(f)
    data=d['PIXEL WINDOW'].data

#    ps=readWMAPPowerSpectrum()
#    print 'p',numpy.shape(ps)
#    print 'd',numpy.shape(data[:,numpy.newaxis][1])
    index={'T':0,'P':1}

    mode='P'
#    print mode,numpy.shape(data[:,numpy.newaxis][index[mode]])

    window=integratePowerSpectrum(data[:,numpy.newaxis][index[mode]],\
                                  bins,ellmin=0,verbose=False)
    return window

#-------------------------------------------------------------------

def integratePowerSpectrum(ps,bins,ellmin=2,verbose=True,dat=False):

    """
    Bin a (fits-imported) WMAP power spectrum using the defined bins
    """

    import numpy

    (nell,ncols)=numpy.shape(ps)
    if verbose: print 'Power spectrum is %i x %i' % (nell,ncols)
    (nbins,dummy)=numpy.shape(bins)
    if verbose: print 'There are %i bins' % nbins

    # Set up the binned power spectrum and copy in the ell-bin values
    bps=numpy.zeros((nbins,ncols+2))
    bps[:,:2]=bins[:,:]

    for ibin in range(nbins):
        # Fetch the (copied) bin ranges
        ellstart,ellstop=bps[ibin,:2]
        #print ellstart,ellstop,ps[ellstart-ellmin:ellstop+1-ellmin,0]
        for icol in range(ncols):
            bps[ibin,icol+2]=numpy.mean(ps[ellstart-ellmin:ellstop+1-ellmin,icol-1])
        # This is an ell alignment issue
        if ellstop < ellmin:
            bps[ibin,2:]=0.0

    if verbose:
        for ibin in range(nbins):
            for icol in range(ncols+1):
                print bps[ibin,icol],
            print
        print 'The binned power spectrum is %i x %i' % numpy.shape(bps)

    return bps

#-------------------------------------------------------------------

def integratePowerSpectrumDat(ps,bins,dat=True):

    """
    Bin a (dat-imported) WMAP power spectrum using the defined bins
    Eventually to be merged with integratePowerSpectrum
    """

    import numpy
    
    (nell,ncols)=numpy.shape(ps)
    (nbins,dummy)=numpy.shape(bins)
    
    #print nell,nbins,ncols,bins
    
    bps=numpy.zeros((nbins,ncols+1))
    bps[:,:2]=bins[:,:] # Copy ell bins
    
    ellmin=2
    for ibin in range(nbins):
        ellstart,ellstop=bps[ibin,:2]
        print ellstart,ellstop,numpy.mean(ps[ellstart-ellmin:ellstop+1-ellmin,2])
        for icol in range(1,ncols):
            bps[ibin,icol+1]=numpy.mean(ps[ellstart-ellmin:ellstop+1-ellmin,icol])
            if ellstop < ellmin:
                bps[ibin,2:]=0.0

    for ibin in range(nbins):
        for icol in range(ncols+1):
            print bps[ibin,icol],
        print

    return bps

#-------------------------------------------------------------------

def readWMAPPowerSpectrum(f='/data4/quiet/jonathaz/map2cl_setup/wmap_lcdm_sz_lens_wmap5_cl_v3.dat'):

    """
    Columns are
    ell TT EE TE
    """

    import numpy

    data=numpy.genfromtxt(f)
    #(nell,ncols)=numpy.shape(data)

    return data

#-------------------------------------------------------------------

def readBins(d='bins_d50_n256.dat'):

    """
    Read bin definition file
    """

    import numpy

    b=open(d,'r')
    x=b.readlines()
    b.close()

    nbins=len(x)-2
    bins=numpy.zeros((nbins,2))

    ibin=0
    for i in x[:]:
        #if i.startswith('#'):
        #    x.remove(i)
        try:
            bins[ibin,:]=i.split()
            ibin+=1
        except:
            continue

    return bins

#-------------------------------------------------------------------

def chisqCalculator(data,model,ndiscardlow=3,ndiscardhigh=1,ndof=None):

    """
    For two arrays of values, and a list of errors,
    calculate the chisq between them
    """

    import numpy

    (nbins,dummy)=numpy.shape(data)
    print dummy

    chisq={'EE':0.0,'BB':0.0}
    n={'EE':1,'BB':3} # from zero
    m={'EE':-2,'BB':3} # from zero

    theory={'EE':model[ndiscardlow:-ndiscardhigh],\
            'EB':numpy.zeros(numpy.shape(model[ndiscardlow:-ndiscardhigh])),\
            'BB':numpy.zeros(numpy.shape(model[ndiscardlow:-ndiscardhigh]))}

    chisq_red={'EE':0.0,'BB':0.0}
    for ibin in range(nbins):
        for mode in chisq.keys():
            d=data[ibin,n[mode]]
            e=data[ibin,n[mode]+1]
            t=theory[mode][ibin,m[mode]]
            c=(d-t)/e
            chisq[mode]+=c*c
            chisq_red[mode]+=c*c/(ndof-1.0)

    return chisq,chisq_red

#-------------------------------------------------------------------

def calculateChisq(cls,iselect):

    """
    For map2cl spectra, return chisq-to-null
    """

    import numpy

    (numchains,numbins,ncols)=numpy.shape(cls)
    numspec=2
    chi=numpy.zeros((numchains,numbins,numspec))
    chisq=numpy.zeros((numchains,numbins,numspec))
    chisqtot=numpy.zeros((numchains,numspec))

    N=numbins-1

    pow={}; dpow={};
    print '# chain bin mode chisq'
    for ichain in range(numchains):
        for ibin in range(numbins):
            ell,pow['EE'],dpow['EE'],pow['BB'],dpow['BB']=cls[ichain,ibin,:]
            for mode in pow.keys():
                mind=ModesLookupMap2Cl[mode]
                chi[ichain,ibin,mind]=pow[mode]/dpow[mode]
                chisq[ichain,ibin,mind]\
                     =chi[ichain,ibin,mind]*chi[ichain,ibin,mind]
                chisqtot[ichain,mind]+=chisq[ichain,ibin,mind]
                if ichain==iselect:
                    if abs(chi[ichain,ibin,mind]) >= 3.0: print '***'
                    print 'chichi %i %i %i %s %f %f'\
                          %(ichain,ibin,mind,mode,chisq[ichain,ibin,mind],\
                            chi[ichain,ibin,mind])

        if ichain==iselect:
            print '#red %i %i %f %f' % (ichain,mind,chisqtot[ichain,mind],chisqtot[ichain,mind]/float(N))

    return chisq,chi

#-------------------------------------------------------------------

def calculatePr(cls,iselect):

    """
    For a map2cl spectrum, calculate P's
    """

    (numchain,numbin,ncols)=numpy.shape(cls)
    numspec=len(modes)
    pr=numpy.zeros((numchain,numspec,numbin))

    print numpy.shape(cls)
    # working here

    ColsLookup={'EE':1,'BB':3}

    for ibin in range(numbin):
        for chain in range(numchain):
            for mode in modes:
                mind=ModesLookupMap2Cl[mode]
                icol=ColsLookup[mode]
                # Mask where, for that mode and bin, realns > this one
                mask=(numpy.where(cls[:,ibin,icol]>cls[chain,ibin,icol]))
                pr[chain,mind,ibin]=float(numpy.shape(cls[mask])[0])/\
                                     float((numchain-1))
                if (pr[chain,mind,ibin] > 0.99): pr[chain,mind,ibin] = 0.99
                if (pr[chain,mind,ibin] < 0.01): pr[chain,mind,ibin] = 0.01

                if chain==iselect:
                    print chain,mind,mode,ibin,pr[chain,mind,ibin]

    return pr

#-------------------------------------------------------------------

def calculateSpectraStats(cls,iselect):

    """
    Wrapper for calculation of stats for map2cl output spectra
    """

    chisq,chi=calculateChisq(cls,iselect)
    pr=calculatePr(cls,iselect)
    dks,pks,pks_sims=calculateKS(pr,iselect)

    return chisq,chi,pr,dks,pks,pks_sims

#-------------------------------------------------------------------

def calculateKS(pr,iselect):

    """
    given P, calculate D and P(D)
    """

    import numpy,sys

    numchain,numspec,numbin=numpy.shape(pr)
    print numchain,numspec,numbin

    dks=numpy.zeros((numchain))
    pks=numpy.zeros((numchain))
    pks_sims=numpy.zeros((numchain))

    for ichain in range(numchain):
        # For both modes, and all bins
        EE=pr[ichain,ModesLookupMap2Cl['EE'],:].flatten().tolist()
        BB=pr[ichain,ModesLookupMap2Cl['BB'],:].flatten().tolist()
        dks[ichain],pks[ichain]=doKSTest(EE+BB)

    # New loop as need all the dks
    for ichain in range(numchain):
        pks_sims[ichain]=FetchP(dks.flatten().tolist(),dks[ichain])
        if ichain==iselect:
            print 'For chain %i,' % ichain
            print EE+BB
            print 'KS stat = %f' % dks[ichain]
            print 'KS prob = %f' % pks[ichain]
            print 'KS sims = %f' % pks_sims[ichain]

    return dks,pks,pks_sims

#-------------------------------------------------------------------

def plotChisqs(chisq,iselect,jackknife):

    """
    Plot chisqs to file for map2cl realisations
    """

    import numpy
    import matplotlib.pyplot as plt
    from scipy import stats

    # Plot bin histograms
    fig=plt.figure(figsize=(11,11),dpi=150) # Letter size
    fig.suptitle('')
    #xmax=numpy.max(chisq)
    xbins = numpy.arange(0,100,1.0)

    weights=1.0e-2*numpy.ones(numpy.shape(chisq))

    for imode,mode in enumerate(modes):
        mind=ModesLookupMap2Cl[mode]
        ax = fig.add_subplot(len(modes),1,imode+1)
        nsim, bins, dummy = ax.hist(chisq[iselect:,:,mind].flatten(),bins=xbins,color='blue',weights=weights[iselect:,:,mind].flatten(),alpha=0.5)
        n, bins, dummy = ax.hist(chisq[iselect,:,mind].flatten(),bins=xbins,color='red',alpha=0.5)
        #print sum(n), sum(nsim)
        ax.set_xlabel('chisq (%s)' % mode)
        plt.ylim(0,5)

    plotf='chisq_%s_%i_m2l.png' % (jackknife,iselect)
    plt.savefig(plotf)
    print 'Look at %s' % plotf

    return

#-------------------------------------------------------------------

def scanDirs(root,pattern='sims$',regexp='^map2cl',nside=None):

    """
    For a root directory(s), return a list of contained map2cl directories
    """

    import os,re,glob

    assert(nside is not None), 'nside required!'

    if nside==256:
        ddd=glob.glob(root)

    elif nside==128:
        if ',' in root:
            roots=root.split(',')
        else:
            roots=[root]

        pre=re.compile(pattern)
        qre=re.compile(regexp)

        # Would love to double [x for y in z] this
        dirs=[]
        for rr in roots:
            dirs.extend([os.path.join(rr,d) for d in os.listdir(rr) if \
                         os.path.isdir(os.path.join(rr,d)) and pre.search(d)])
        dirs.sort()

        ddd=[]
        for d in dirs:
            dd = [os.path.abspath(os.path.join(d,e)) for e in os.listdir(d) if \
                  os.path.isdir(os.path.join(d,e)) and qre.search(e)]
            if len(dd) > 1:
                print 'problem - too many subdirs', dd
                # print 'trying to continue...'
                # jkdict={0:'0',1:'1'}
                # return dd,jkdict
                return
            ddd.append(dd[0])

    # Now establish the jkdict
    jkdict={}
    pattern3='(?<=%s)\w+(?=_)'%long_substr([os.path.basename(d) for d in ddd])
    rre=re.compile(pattern3)
    for ijk,d in enumerate(ddd):
        try:
            jkdict[ijk]=rre.search(os.path.basename(d)).group(0)
        except:
            print 'Error in regexp search for jackknife name'

    return ddd,jkdict

#-------------------------------------------------------------------

def readPowspecs(dir,numbins=7,ncols=5,ndiscardlow=1,ndiscardhigh=1,\
                 pattern='cls_real\d{4,4}.dat',version=0,nside=None,\
                 override_chains=None):

    """
    Read map2cl-power spectra from a directory
    """

    import os,re,numpy

    pre=re.compile(pattern)

    files=os.listdir(dir)
    files=[f for f in files if pre.match(f)]
    files.sort()

    nchains=len(files)

    print 'The number of files found is %i' % nchains

    if nside==256:
        if version==2 or version==3:
            numbins=12
        else:
            numbins=10
    elif nside==128:
        numbins=7

    if override_chains is not None:
        nchains=override_chains
        print 'We are overriding the number of chains (%i)' % nchains

    #nchains=1
    #nchains = 10

    cls=numpy.zeros((nchains,numbins-ndiscardlow-ndiscardhigh,ncols))

    #print 'Reading pattern %s from %s' % (pattern,dir)

    for ichain,f in enumerate(files):
        #print 'THIS IS THE PATH %s' % os.path.join(dir,f)
        if ichain >= nchains: break # Exit after user-supplied nchains
        cls[ichain,:,:]=readPowspec(os.path.join(dir,f),version=version)

    return cls

#-------------------------------------------------------------------

def readPowspec(f,ndiscardlow=1,ndiscardhigh=1,version=0,\
                override_chains=None,verbose=False):

    """
    Read map2cl-style powspec file to numpy array cls, with columns:
    ell, EE, dEE, BB, dBB
    Top and bottom bins are stripped out if required
    optional version 1:
        ##  l_bin    l_1   l_2       EE      sigma_EE        BB       sigma_BB
        # 13.500     2    25    154.5631     65.9872     19.7837     15.1406
    option version 2:
        as v1, but with all spectra present as:
        TT TE TB EE EB BB
    option version 3:
        as v2, but now retain EB

    0  ell_bin
    1  ell_min
    2  ell_max
    3  TT
    5  TE
    7  TB
    9  EE
    11 EB
    13 BB
    """

    import numpy
    print 'PRINT THIS'
    data=numpy.genfromtxt(f)
    (nbins,ncols)=numpy.shape(data)
    verbose=True
    if version == 3:
        cols=[0,9,10,11,12,13,14]
        ndiscardlow=0
        ndiscardhigh=0
    if version == 2:
        cols=[0,9,10,13,14]
        ndiscardlow=0
        ndiscardhigh=0
    elif version == 1:
        cols=[0,3,4,5,6]
        ndiscardlow=1
        ndiscardhigh=1
    else:
        cols=numpy.arange(ncols) # i.e. all columns

    if verbose:
        for b in range(nbins): # was 5
            for mode in ['EE','EB','BB']:
                im=2*ModesLookup[mode]+3; id=im+1
                print '%ssigma %s %f %f' \
                      % (mode,f,data[b,0],data[b,im]/data[b,id])
                #print 'A %ssigma %s %f %f' \
                #      % (mode,f,data[b,0],data[b,im]/data[b,id])

    cls=data[ndiscardlow:nbins-ndiscardhigh,cols]

    return cls

#-------------------------------------------------------------------

def readCLs(f):

    """
    Read cls.fits into numpy data array
    """

    import pyfits

    # Read FITS
    fits=pyfits.open(f)
    header=fits[0].header
    data=fits[0].data

    # Fetch key parameters
    lmax=header['LMAX']
    numsamp=header['NUMSAMP']
    numchain=header['NUMCHAIN']
    numspec=header['NUMSPEC']

    print '\nFor %s, lmax = %i, numsamp = %i, numchain = %i, numspec = %i\n'\
          % (f,lmax,numsamp,numchain,numspec)

    fstem=f.split('.')[0]

    return data,lmax,numsamp,numchain,numspec,fstem

#-------------------------------------------------------------------

def readKSDS(ksdsf,numjks,numchain):

    """
    Read ks-ds.txt back in for histogramming
    """

    import numpy,scanf

    dks=numpy.zeros((numjks,numchain))

    data=open(ksdsf,'r')

    linef= '%i %s %i %f %f %f %f'

    jknames={}
    for line in data:
        if line.startswith('#'): continue
        ijk,jk,ichain,d,pks,dummy,psims=scanf.sscanf(line,linef)
        if ijk < 0: continue
        dks[ijk,ichain]=d
        jknames[str(ijk)]=jk

    jknames['all']='all'

    data.close()

    return dks,jknames

#-------------------------------------------------------------------

def probks(alam):

    """
    KS distribution (Q_KS)
    From Numerical Recipes
    """

    from math import exp

    eps1 = 0.001
    eps2 = 1.0e-8

    a2 = -2.0*alam*alam
    fac = 2.0
    probks = 0.0
    termbf = 0.0

    for j in range(1,100+1):
        term = fac * exp(a2*j*j)
        probks += term
        if (abs(term) <= eps1*termbf or abs(term) <= eps2*probks):
            return probks
        fac = -fac
        termbf = abs(term)
    # If failed to converge
    print '**WARNING - probks failed to converge'
    probks = 1.0
    return probks

#-------------------------------------------------------------------

def ks_uniform(x):

    """
    KS test for uniform distribution
    """

    from math import sqrt

    y=sorted(x)
    n=float(len(y))

    d=0.0
    fo=0.0

    for j in range(1,int(n)+1):
        fn = float(j) / n
        ff = y[j-1]
        dt = max([abs(fo-ff),abs(fn-ff)])
        if (dt > d): d = dt
        fo = fn

    arg = (sqrt(n) + 0.12 + 0.11/sqrt(n))*d
    prob = probks(arg)

    return (d, prob)


#-------------------------------------------------------------------

def p(d,n):

    """
    Given d, and n, return theoretical P
    P(D) for a theoretical distribution is usually calculated by the KS
    routines at the same time as D
    In this function, P(D) can be calculated post-hoc, given also N
    """

    from math import sqrt

    arg = (sqrt(n) + 0.12 + 0.11/sqrt(n))*d
    prob = probks(arg)
                    
    return prob

#-------------------------------------------------------------------

# Define EE and BB keys
ModesLookup={'TT':0,'TE':1,'TB':2,'EE':3,'EB':4,'BB':5}
ModesLookupMap2Cl={'EE':0,'BB':1}
modes=['EE','BB']

#-------------------------------------------------------------------

def reducePowerSpec(datain, numbins=7):

    """
    Take out low and high modes from power spectra data
    """
    # use like this:
    #dataell,dataoutEE,dataoutEEerror,dataoutBB,dataoutBBerror=reducePowerSpec(cls,7)

    import numpy,math


    (numchain,dummy,dummy)=numpy.shape(datain)
    numchain=11

    dataell=numpy.zeros((numbins-2,numchain+1))
    dataoutEE=numpy.zeros((numbins-2,numchain+1))
    dataoutEEerror=numpy.zeros((numbins-2,numchain+1))
    dataoutBB=numpy.zeros((numbins-2,numchain+1))
    dataoutBBerror=numpy.zeros((numbins-2,numchain+1))
    EEsigmaaway=numpy.zeros((numbins-2,numchain+1))
    BBsigmaaway=numpy.zeros((numbins-2,numchain+1))
    minmaxvalues=numpy.zeros((2,4,numchain+1))
    EEmin=numpy.zeros((numbins-2,numchain+1))
    EEmax=numpy.zeros((numbins-2,numchain+1))
    BBmin=numpy.zeros((numbins-2,numchain+1))
    BBmax=numpy.zeros((numbins-2,numchain+1))
    tempEEmin=numpy.zeros((numbins-2,numchain+1))
    tempEEmax=numpy.zeros((numbins-2,numchain+1))
    tempBBmin=numpy.zeros((numbins-2,numchain+1))
    tempBBmax=numpy.zeros((numbins-2,numchain+1))
    tempEEsigma=numpy.zeros((numbins-2,numchain+1))
    tempBBsigma=numpy.zeros((numbins-2,numchain+1))

    # Store means in a cube, one value per mode, bin and realisation
    for chain in range(numchain):
        a=0
        for ell_value in range(numbins-2):
            if ((float(datain[chain,ell_value,0]) > 45) and (float(datain[chain,ell_value,0]) < 345)):
                dataell[a,chain] = datain[chain,ell_value,0]
                dataoutEE[a,chain] = datain[chain,ell_value,1]
                dataoutEEerror[a,chain] = datain[chain,ell_value,2]
                dataoutBB[a,chain] = datain[chain,ell_value,3]
                dataoutBBerror[a,chain] = datain[chain,ell_value,4]
                EEsigmaaway[a,chain] = datain[chain,ell_value,1]/datain[chain,ell_value,2]
                BBsigmaaway[a,chain] = datain[chain,ell_value,3]/datain[chain,ell_value,4]
                EEmin[a,chain] = datain[chain,ell_value,1] - datain[chain,ell_value,2]
                EEmax[a,chain] = datain[chain,ell_value,1] + datain[chain,ell_value,2]
                BBmin[a,chain] = datain[chain,ell_value,3] - datain[chain,ell_value,4]
                BBmax[a,chain] = datain[chain,ell_value,3] + datain[chain,ell_value,4]
                if(chain > 0.5):
                    dataell[a,numchain]+=datain[chain,ell_value,0]/float(numchain-1)
                    dataoutEE[a,numchain]+=datain[chain,ell_value,1]/float(numchain-1)
                    dataoutEEerror[a,numchain]+=datain[chain,ell_value,2]/float(numchain-1)
                    dataoutBB[a,numchain]+=datain[chain,ell_value,3]/float(numchain-1)
                    dataoutBBerror[a,numchain]+=datain[chain,ell_value,4]/float(numchain-1)
                    #EEsigmaaway[a,numchain]+=(datain[chain,ell_value,1]/datain[chain,ell_value,2])/(float(numchain-1)
                    #BBsigmaaway[a,numchain]+=(datain[chain,ell_value,3]/datain[chain,ell_value,4])/float(numchain-1)
                a+=1

    



    for chain in range(numchain):
        for ell_value in range(numbins -2):
            tempEEmin[ell_value,chain] = datain[chain,ell_value,1] - datain[chain,ell_value,2]        
	    tempEEmax[ell_value,chain] = datain[chain,ell_value,1] + datain[chain,ell_value,2]        
	    tempBBmin[ell_value,chain] = datain[chain,ell_value,3] - datain[chain,ell_value,4]        
            tempBBmax[ell_value,chain] = datain[chain,ell_value,3] + datain[chain,ell_value,4]        
            tempEEsigma[ell_value,chain] = datain[chain,ell_value,1]/datain[chain,ell_value,2]
            tempBBsigma[ell_value,chain] = datain[chain,ell_value,3]/datain[chain,ell_value,4]

        minmaxvalues[0,0,chain] = math.floor(min(tempEEmin[:,chain]))
        minmaxvalues[1,0,chain] = math.ceil(max(tempEEmax[:,chain]))
        minmaxvalues[0,1,chain] = math.floor(min(tempBBmin[:,chain]))
        minmaxvalues[1,1,chain] = math.ceil(max(tempBBmax[:,chain]))
        minmaxvalues[0,2,chain] = math.floor(min(tempEEsigma[:,chain]))
        minmaxvalues[1,2,chain] = math.ceil(max(tempEEsigma[:,chain]))
        minmaxvalues[0,3,chain] = math.floor(min(tempBBsigma[:,chain]))
        minmaxvalues[1,3,chain] = math.ceil(max(tempBBsigma[:,chain]))


    for ell_value in range(numbins -2):
        print 'EEsigma ', EEsigmaaway[ell_value,0]
        print 'BBsigma ', BBsigmaaway[ell_value,0]
    for chain in range(1,10):
        for ell_value in range(numbins -2):
            print 'EEsigma',chain,' ', EEsigmaaway[ell_value,chain]
            print 'BBsigma',chain,' ', EEsigmaaway[ell_value,chain]
    #print 'EEsigmaavg ',  EEsigmaaway[ell_value,numchain]
    #print 'BBsigmaavg ',  BBsigmaaway[ell_value,numchain]


    return dataell,dataoutEE,dataoutEEerror,dataoutBB,dataoutBBerror,EEsigmaaway,BBsigmaaway, minmaxvalues


#-------------------------------------------------------------------

def calculateBurnin(data,lmin,lmax,dl,burnin,numchain,numbins,numspec,fstem):

    """
    For a cls array (from FITS),
    calculate stats by bin and
    write to file
    """

    import numpy
    from scipy import stats

    sumburn=numpy.zeros((numchain,numspec,numbins,numsamp))
    nosumburn=numpy.zeros((numchain,numspec,numbins,numsamp))
    # Store means in a cube, one value per mode, bin and realisation

    for ibin,ell in enumerate(range(lmin,lmax+dl,dl)):
        for chain in range(numchain):
            for mode in modes:
                mind=ModesLookup[mode]
                for sample in range(numsamp):
                    if (sample < 2):
                        if (sample == 0):
                            sumburn[chain,mind,ibin,sample]=data[(sample+1),chain,mind,ell]/(float(sample) +1.0)
                            nosumburn[chain,mind,ibin,sample]=data[(sample+1),chain,mind,ell]/1.0
                            #print 'This is it %f sample %d chain %d mind %d ell %d' % (data[(sample+1),chain,mind,ell], sample, chain, mind, ell)
                        if (sample == 1):
                            sumburn[chain,mind,ibin,sample]=data[sample,chain,mind,ell]
                            nosumburn[chain,mind,ibin,sample]=data[sample,chain,mind,ell]
                            #print 'This is it %f sample %d chain %d mind %d ell %d' % (data[sample,chain,mind,ell], sample, chain, mind, ell)
                        #sumburn[chain,mind,ibin,sample]=1 #test clause keep
			#print '%f %d' % (sumburn[chain,mind,ibin,sample], sample)
						#sumburn[chain,mind,ibin,sample]=data[sample,chain,mind,ell]/(float(sample) +1.0)
                    else:
                        #sumburn[chain,mind,ibin,sample]=((sumburn[chain,mind,ibin,(sample -1)]*float(sample)) + 0)/(float(sample)+1.0) # test clause keep
                        #print 'sample %d chain %d mind %d ell %d' % (sample, chain, mind, ell)
                        #print 'This is it %f sample %d chain %d mind %d ell %d' % (data[sample,chain,mind,ell], sample, chain, mind, ell)
                        sumburn[chain,mind,ibin,sample]=((sumburn[chain,mind,ibin,(sample -1)]*float(sample)) + data[sample,chain,mind,ell])/(float(sample)+1.0)
                        nosumburn[chain,mind,ibin,sample]=data[sample,chain,mind,ell]
			#print '%f %d' % (sumburn[chain,mind,ibin,sample], sample)
    #print 'Sumburn length'
    #print sumburn.shape
    #print sumburn[0,3,2,200]
    return sumburn,nosumburn

#-------------------------------------------------------------------

def mode_value(data):
    """
    From http://projects.scipy.org/scipy/ticket/905
    """
    counts = {}
    for x in data.flatten():
        counts[x] = counts.get(x,0) + 1
        maxcount = max(counts.values())
        modelist = []
        for x in counts:
            if counts[x] == maxcount:
                modelist.append(x)
    return modelist,maxcount

#-------------------------------------------------------------------

def writeMeanSpectrum(cls,idata,jname):

    """
    For one jackknife, write mean cl and error to file
    """

    import math,numpy

    print "This function isn't working - aborting"
    return None

    (numchains,numbins,numcols)=numpy.shape(cls)
    
    print 'numchains = %i' % numchains
    print 'numbins   = %i' % numbins
    print 'numcols   = %i' % numcols

    mu=numpy.zeros((numbins,2))
    std=numpy.zeros((numbins,2))

    outf='powspec.txt'
    out=open(outf,'w')

    header='# bin ell EE dEE BB dBB\n'
    out.write(header)
    for ibin in range(numbins):
        ellc=cls[0,ibin,0]
        for mode in modes:
                mind=ModesLookupMap2Cl[mode]
                col=2*(mind+1)
                mu[ibin,mind]=numpy.mean(cls[:,ibin,col])
                std[ibin,mind]=numpy.std(cls[:,ibin,col])
#                print mode,mind,mu[ibin,mind],std[ibin,mind]
        line= '%i %i %f %f %f %f' \
              % (ibin,ellc,mu[ibin,0],std[ibin,0],mu[ibin,1],std[ibin,1])
#        print line
        out.write('%s\n'%line)

    print 'Look in %s' % outf
    out.close()

    return

#-------------------------------------------------------------------


def calculateStats(data,lmin,lmax,dl,burnin,\
                   numchain,numbins,numspec,fstem,dosum):

    """
    For a cls array (from FITS),
    calculate stats by bin and
    write to file
    """

    import numpy
    from scipy import stats

    # Store means in a cube, one value per mode, bin and realisation
    mu=numpy.zeros((numchain,numspec,numbins))
    std=numpy.zeros((numchain,numspec,numbins))
    p=numpy.zeros((numchain,numspec,numbins))
    pr=numpy.zeros((numchain,numspec,numbins))
    cupper=numpy.zeros((numchain,numspec,numbins))
    clower=numpy.zeros((numchain,numspec,numbins))
    map=numpy.zeros((numchain,numspec,numbins))
    median=numpy.zeros((numchain,numspec,numbins))
    error=numpy.zeros((numchain,numspec,numbins))
    chisq=numpy.zeros((numchain,numspec,numbins))
    chisq2=numpy.zeros((numchain,numspec,numbins))
    chi=numpy.zeros((numchain,numspec,numbins))

    per_low = 16 # } 68% about median
    per_up  = 84 # } 68% about median

    histf='%s.hist' % fstem
    hist=open(histf,'w')

    # Calculate mean C_l's
    # NB For each ell in a bin, all the C_l values are the same
    #    -> so can ignore integrating across the bin
    hist.write('# bin lmin lmax realisation mu_EE rank_EE size_EE p_EE std_EE error_EE chisq_EE chisq2_EE chi_EE map_EE median_EE mu_BB rank_BB size_BB p_BB std_BB error_BB chisq_BB chisq2_BB chi_BB map_BB median_BB jk\n')
    for ibin,ell in enumerate(range(lmin,lmax+dl,dl)):
        #hist.write('# For bin %i (ell=%i-%i):\n' %(ibin,ell-dl/2,ell+dl/2))
        # For each realisation (data, and sims)
        for chain in range(numchain):
            # For each mode (EE, BB)
            if (dosum and chain==0): continue
            for mode in modes:
                mind=ModesLookup[mode]
                slice=data[burnin:numsamp,chain,mind,ell]
                mu[chain,mind,ibin]=numpy.mean(slice)
                std[chain,mind,ibin]=numpy.std(slice)
                #map[chain,mind,ibin]=float(stats.mode(slice)[0])
                #map[chain,mind,ibin]=mode_value(slice)[0][0]
                h=numpy.histogram(slice.flatten(),bins=100)
                #print h
                map[chain,mind,ibin]=h[1][h[0].argmax()]
                median[chain,mind,ibin]=numpy.median(slice,axis=0)
                clower[chain,mind,ibin]=stats.scoreatpercentile(slice,per_low)
                cupper[chain,mind,ibin]=stats.scoreatpercentile(slice,per_up)
                error[chain,mind,ibin]=(cupper[chain,mind,ibin]-clower[chain,mind,ibin])/2.0
                chisq[chain,mind,ibin]=(mu[chain,mind,ibin]/error[chain,mind,ibin])**2.0
                chisq2[chain,mind,ibin]=(mu[chain,mind,ibin]/std[chain,mind,ibin])**2.0
                chi[chain,mind,ibin]=(mu[chain,mind,ibin]/error[chain,mind,ibin])
        # Count number of mu's (for sims) smaller than mu_obs
        #  -- that number - call p(bin,spec,jk)
        # Can't combine these loops, because need all mu's available
        for chain in range(numchain):
            if (dosum and chain==0): continue
            hist.write('%i %i %i %i ' % (ibin,ell-dl/2,ell+dl/2,chain))
            for mode in modes:
                mind=ModesLookup[mode]
                # Mask where, for that mode and bin, realns > this one
                mask=(numpy.where(mu[:,mind,ibin]<mu[chain,mind,ibin]))
                p[chain,mind,ibin]=numpy.shape(mu[mask])[0]
                pr[chain,mind,ibin]=float(p[chain,mind,ibin])/float((numchain-1))
                if (pr[chain,mind,ibin] > 0.99):
                    pr[chain,mind,ibin] = 0.99
                if (pr[chain,mind,ibin] < 0.01):
                    pr[chain,mind,ibin] = 0.01
                hist.write('%f '    % mu[chain,mind,ibin])
                hist.write('%i '    % p[chain,mind,ibin])
                hist.write('%i '    % (numchain-1))
                hist.write('%3.2f ' % pr[chain,mind,ibin])
                hist.write('%f '    % std[chain,mind,ibin])
                hist.write('%f '    % error[chain,mind,ibin])
                hist.write('%f '    % chisq[chain,mind,ibin])
                hist.write('%f '    % chisq2[chain,mind,ibin])
                hist.write('%f '    % chi[chain,mind,ibin])
                hist.write('%f '    % map[chain,mind,ibin])
                hist.write('%f '    % median[chain,mind,ibin])
            hist.write('%s\n'   % f.split('.')[0])


    print 'Look in %s' % histf
    hist.close()

    return mu,p,pr,std,map,median,clower,cupper,error,chisq,chisq2,chi

#-------------------------------------------------------------------

def calculateCumStats(pr):

    """
    For a cls array (from FITS),
    calculate thecumstats
    write to file
    """

    import numpy
    import math
    from scipy import stats

    # Store means in a cube, one value per mode, bin and realisation

    cumdist=numpy.zeros((numchain,numspec,100))
    cumdistlength=numpy.zeros((numchain,numspec))
    cumdistposition=numpy.zeros((numchain,numspec))

    for chain in range(numchain):
        for mode in modes:
            mind=ModesLookup[mode]
            for j in range(0,100,1):
                for k in range(numbins):
                    #print pr[chain,mind,k], ' HERE IT IS chain ', chain, ' mode ', mind, ' j ', j, ' cumdist at this point ', cumdist[chain,mind,j]
                    if((pr[chain,mind,k] * 100) < (j + 1)):
                        cumdist[chain,mind,j] = cumdist[chain,mind,j] + 1.0/float(numbins);
                        #if (chain == 0):
                            #if (mind == 3):                        
                                #print j, ' ', cumdist[chain,mind,j]
    #print cumdist[0,3,:]

    for chain in range(numchain):
        for mode in modes:
            mind=ModesLookup[mode]
            largestdist = -1
            largestdistposition = -1
            thisdist = -1
            for j in range(0,100,1):
                thisdist = math.fabs(cumdist[chain,mind,j] - ((float(j)+0.0)/100.0))
                #if (chain == 0):
                    #if (mind == 3):
                        #print 'j ', j, ' this distance ', thisdist, ' cumdistposition  ', cumdistposition[chain,mind]
                if(thisdist > largestdist):
                    cumdistlength[chain,mind] = thisdist
                    largestdist = thisdist
                    thisdist = cumdistlength[chain,mind]
                    cumdistposition[chain,mind] = j
                    #print 'this distance ', thisdist, ' cumdistposition  ', cumdistposition[chain,mind]

    return cumdist,cumdistlength,cumdistposition

#-------------------------------------------------------------------

def plotPowerSpecMap2Cl(data,patch_name,filestem,iselect):

    """
    Plot power spectrum to file for selected chain
    """

    #to be called as such, in order to plot the power spectra for the data realisation, plus the first two simulations
    #plotPowerSpecMap2Cl(dataell,dataoutEE,dataoutEEerror,dataoutBB,dataoutBBerror,patch_name,filestem,0)
    #plotPowerSpecMap2Cl(dataell,dataoutEE,dataoutEEerror,dataoutBB,dataoutBBerror,patch_name,filestem,1)
    #plotPowerSpecMap2Cl(dataell,dataoutEE,dataoutEEerror,dataoutBB,dataoutBBerror,patch_name,filestem,2)

    import numpy
    import matplotlib.pyplot as plt

    # Plot Power Spectra
    fig=plt.figure(figsize=(8.5,11),dpi=150) # Letter size
    fig.suptitle('Patch %s, %s %i'%(patch_name,filestem,iselect), fontsize='x-large')

    ax = fig.add_subplot(2,1,1)
    ax.set_xlabel('x label', fontsize='x-large')
    ax.set_xlabel('ell')
    ax.set_ylabel('y label', fontsize='x-large')
    ax.set_ylabel('C_ell/uK^2 (EE)')
    ax.errorbar(data['ell'][:,iselect], data['EE'][:,iselect],xerr=0,yerr=[data['EEerror'][:,iselect],data['EEerror'][:,iselect]],marker='x',color='red',markersize=25,linestyle='None')
    ax.axhline(y=0,linestyle='-.')
    ax.set_xlim((0.0,350.0))
    ax.set_ylim((data['minmaxvalues'][0,0,iselect],data['minmaxvalues'][1,0,iselect]))
    ax = fig.add_subplot(2,1,2)
    ax.set_xlabel("x label", fontsize='x-large')
    ax.set_xlabel('ell')
    ax.set_ylabel("y label", fontsize='x-large')
    ax.set_ylabel('C_ell/uK^2 (BB)')
    ax.errorbar(data['ell'][:,iselect], data['BB'][:,iselect],xerr=0,yerr=[data['BBerror'][:,iselect],data['BBerror'][:,iselect]],marker='x',color='red',markersize=25,linestyle='None')
    ax.axhline(y=0,linestyle='-.')
    ax.set_xlim((0.0,350.0))
    ax.set_ylim((data['minmaxvalues'][0,1,iselect],data['minmaxvalues'][1,1,iselect]))
    #ax.set_ylim((-20,20.0))

    plotf='%s_%03i_map2cl_powspec.png' % (filestem,iselect)
    print 'look at %s_%03i_map2cl_powspec.png' % (filestem,iselect)
    plt.savefig(plotf)

    # Plot Power Spectra
    fig2=plt.figure(figsize=(8.5,11),dpi=150) # Letter size
    fig2.suptitle('Patch %s, %s %i'(patch_name,filestem,iselect), fontsize='x-large')

    ax = fig2.add_subplot(2,1,1)
    ax.set_xlabel('x label', fontsize='x-large')
    ax.set_xlabel('ell')
    ax.set_ylabel('y label', fontsize='x-large')
    ax.set_ylabel('sigma away (EE)')
    ax.errorbar(data['ell'][:,iselect], data['EEsigmaaway'][:,iselect],xerr=0,yerr=0,marker='x',color='red',markersize=25,linestyle='None')
    ax.axhline(y=0,linestyle='-.')
    ax.set_xlim((0.0,350.0))
    ax.set_ylim((data['minmaxvalues'][0,2,iselect],data['minmaxvalues'][1,2,iselect]))
    #ax.set_ylim((-20,20.0))
    ax = fig2.add_subplot(2,1,2)
    ax.set_xlabel('x label', fontsize='x-large')
    ax.set_xlabel('ell')
    ax.set_ylabel('y label', fontsize='x-large')
    ax.set_ylabel('sigma away (BB)')
    ax.errorbar(data['ell'][:,iselect], data['BBsigmaaway'][:,iselect],xerr=0,yerr=0,marker='x',color='red',markersize=25,linestyle='None')
    ax.axhline(y=0,linestyle='-.')
    ax.set_xlim((0.0,350.0))
    ax.set_ylim((data['minmaxvalues'][0,3,iselect],data['minmaxvalues'][1,3,iselect]))
    #ax.set_ylim((-20,20.0))

    plotf='%s_%03i_map2cl_powspec_sigma.png' % (filestem,iselect)
    print 'look at %s_%03i_map2cl_powspec_sigma.png' % (filestem,iselect)
    plt.savefig(plotf)


    return

#-------------------------------------------------------------------

def plotBurns(noburns,burns,lmin,lmax,dl,numchain,numbins,numspec,\
                  patch_name,filestem,iselect):

    """
    Plot burns of chains to file for selected chain
    """

    import numpy
    import matplotlib.pyplot as plt
    # Plot burns
    fig=plt.figure(figsize=(8.5,11),dpi=150) # Letter size
    fig.suptitle('Patch %s, %s'%(patch_name,filestem))
    #fig.suptitle(filestem)
    #xbins = numpy.arange(0,20,0.2)
    x=numpy.arange(numsamp)
    #x=numpy.array(range(400),dtype=int)
    for ibin,ell in enumerate(range(lmin,lmax+dl,dl)):
        for imode,mode in enumerate(modes):
            mind=ModesLookup[mode]

            thisburn=numpy.zeros(numsamp)
            for sample in range(numsamp):
                thisburn[sample]=burns[iselect,mind,ibin,sample]
#            print 'specific burn shape'        
#            print thisburn.shape[0]


            ax = fig.add_subplot(numbins,len(modes),imode+1+len(modes)*ibin)
            ax.set_ylabel('%i-%i'%(ell-dl/2,ell+dl/2))
            ax.plot(x,thisburn,'r-',linewidth=1)            
            ax.set_xlabel('Burns (%s)'%mode)

    plotf='%s_%03i_burns.png' % (filestem,iselect)
    plt.savefig(plotf)
    print 'Look at %s' % plotf

    fig2=plt.figure(figsize=(8.5,11),dpi=150) # Letter size
    fig2.suptitle('Patch %s, %s'%(patch_name,filestem))
    #fig.suptitle(filestem)
    #xbins = numpy.arange(0,20,0.2)
    x=numpy.arange(numsamp)
    #x=numpy.array(range(400),dtype=int)
    for ibin,ell in enumerate(range(lmin,lmax+dl,dl)):
        for imode,mode in enumerate(modes):
            mind=ModesLookup[mode]

            thisdata=numpy.zeros(numsamp)
            for sample in range(numsamp):
                thisdata[sample]=noburns[iselect,mind,ibin,sample]

            ax = fig2.add_subplot(numbins,len(modes),imode+1+len(modes)*ibin)
            ax.set_ylabel('%i-%i'%(ell-dl/2,ell+dl/2))
            ax.plot(x,thisdata,'r-',linewidth=1)
            ax.set_xlabel('Samples (%s)'%mode)

    plotf='%s_%03i_samples.png' % (filestem,iselect)
    plt.savefig(plotf)
    print 'Look at %s' % plotf

    return

#-------------------------------------------------------------------

def plotHistogram(data,lmin,lmax,dl,burnin,numchain,numbins,numspec,\
                  mu,median,map,clower,cupper,filestem,iselect):

    """
    Plot histograms of chains to file for selected chain
    """

    import numpy
    import matplotlib.pyplot as plt

    # Plot bin histograms
    fig=plt.figure(figsize=(8.5,11),dpi=150) # Letter size
    fig.suptitle(filestem)
    xbins = numpy.arange(0,20,0.2)
    #chisq=stats.chi2.pdf(xbins,3)
    #ax.plot(xbins,50.0*chisq,'r-',linewidth=3)

    for ibin,ell in enumerate(range(lmin,lmax+dl,dl)):
        for imode,mode in enumerate(modes):
            mind=ModesLookup[mode]
            #print imode+1+len(modes)*ibin,ell,mode
            #print 'Drawing bin %i (ell=%i-%i)' %(ibin,ell-dl/2,ell+dl/2)
            ax = fig.add_subplot(numbins,len(modes),imode+1+len(modes)*ibin)
            ax.set_ylabel('%i-%i'%(ell-dl/2,ell+dl/2))
            n, bins, dummy = ax.hist(data[burnin:numsamp,iselect,mind,ell], bins=xbins)
            ax.set_xlabel('C_ell/uK^2 (%s)'%mode)
            ax.axvline(mu[iselect,mind,ibin],color='red')
            ax.axvline(median[iselect,mind,ibin],color='blue')
            ax.axvline(map[iselect,mind,ibin],color='cyan')
            ax.axvline(clower[iselect,mind,ibin],color='yellow')
            ax.axvline(cupper[iselect,mind,ibin],color='yellow')

    plotf='%s_%03i.png' % (filestem,iselect)
    plt.savefig(plotf)
    print 'Look at %s' % plotf

    return

#-------------------------------------------------------------------

def plotPowerSpec(data,lmin,lmax,dl,burnin,numchain,numbins,numspec,\
                  mu,median,clower,cupper,patch_name,filestem,iselect):

    """
    Plot power spectrum to file for selected chain
    """

    import numpy
    import matplotlib.pyplot as plt

    # Plot Power Spectra
    fig=plt.figure(figsize=(8.5,11),dpi=150) # Letter size
    fig.suptitle('Patch %s, %s'%(patch_name,filestem))

    x=numpy.array([50,100,150,200,250,300])
    for imode,mode in enumerate(modes):
        mind=ModesLookup[mode]
        ax = fig.add_subplot(len(modes),1,imode+1)
        ax.set_xlabel('ell')
        ax.set_ylabel('C_ell/uK^2 (%s)'%mode)
        ax.errorbar(x, mu[iselect,mind,:],xerr=0,yerr=[mu[iselect,mind,:]-clower[iselect,mind,:],cupper[iselect,mind,:]-mu[iselect,mind,:]],marker='x',color='red',markersize=25,linestyle='None')
        ax.errorbar(x, median[iselect,mind,:],xerr=0,yerr=[mu[iselect,mind,:]-clower[iselect,mind,:],cupper[iselect,mind,:]-mu[iselect,mind,:]],marker='x',color='blue',markersize=26,linestyle='None')
        ax.axhline(y=0,linestyle='-.')
        ax.set_xlim((0.0,350.0))

    plotf='%s_%03i_powspec.png' % (filestem,iselect)
    plt.savefig(plotf)
    print 'Look at %s' % plotf

    return

#-------------------------------------------------------------------

def plotSummedPowspec(data,lmin,lmax,dl,burnin,numchain,numbins,\
                  numspec,map,patch_name,filestem,iselect):

    """
    Plot power spectrum to file for simulated data
    Ignore the real data!
    """

    import numpy
    import matplotlib.pyplot as plt

    # Plot Power Spectra
    fig=plt.figure(figsize=(8.5,11),dpi=150) # Letter size
    fig.suptitle('Patch %s, %s'%(patch_name,filestem))

    x=numpy.array([50,100,150,200,250,300])
    for imode,mode in enumerate(modes):
        mind=ModesLookup[mode]
        ax = fig.add_subplot(len(modes),1,imode+1)
        ax.set_xlabel('ell')
        ax.set_ylabel('C_ell/uK^2 (%s)'%mode)
        ax.plot(x,map[iselect+1:,mind,:],'r-',linewidth=1)
#        ax.errorbar(x, mu[iselect,mind,:],xerr=0,yerr=[mu[iselect,mind,:]-clower[iselect,mind,:],cupper[iselect,mind,:]-mu[iselect,mind,:]],marker='x',color='red',markersize=25,linestyle='None')
#        ax.errorbar(x, median[iselect,mind,:],xerr=0,yerr=[mu[iselect,mind,:]-clower[iselect,mind,:],cupper[iselect,mind,:]-mu[iselect,mind,:]],marker='x',color='blue',markersize=26,linestyle='None')
        ax.axhline(y=0,linestyle='-.')
        ax.set_xlim((0.0,350.0))

    plotf='%s_%03i_powspec.png' % (filestem,iselect)
    plt.savefig(plotf)
    print 'Look at %s' % plotf

    return

#-------------------------------------------------------------------

def plotProb(pr,cumdist,cumdistlength,cumdistposition,patch_name,filestem,iselect):

    """
    Plot probs to file for selected chain
    """

    import numpy
    import matplotlib.pyplot as plt
    import math

    # Plot bin histograms
    fig=plt.figure(figsize=(8.5,11),dpi=150) # Letter size
    fig.suptitle('Patch %s, %s'%(patch_name,filestem))
    #fig.suptitle(filestem)

    for imode,mode in enumerate(modes):
        mind=ModesLookup[mode]
        ax = fig.add_subplot(len(modes),1,imode+1)
        ax.set_xlabel('Prob(mu_obs > mu_sims) (%s)' % mode)
        x=pr[iselect,mind,:]
        y=numpy.zeros(numpy.shape(x))
        ax.plot(x,y,linestyle='None', marker='x', markersize=10)
        ax.set_xlim((0,1.0))
        ax.set_ylim((-0.05,0.05))

    plotf='%s_%03i_probs.png' % (filestem,iselect)
    plt.savefig(plotf)
    print 'Look at %s' % plotf

    fig2=plt.figure(figsize=(8.5,11),dpi=150) # Letter size
    fig2.suptitle('Patch %s, %s'%(patch_name,filestem))

    for imode,mode in enumerate(modes):
        mind=ModesLookup[mode]
        ax = fig2.add_subplot(len(modes),1,imode+1)
        ax.set_xlabel('Cum Prob(mu_obs > mu_sims) (%s)' % mode)
        xbins = numpy.arange(0,1,0.01)
	xd=x.shape[0]
	xs=xd +1
        xl=xd/100.0
        mypos = cumdistposition[iselect,mind]
        myposd = float(mypos)/100.0
        thisheight = cumdist[iselect,mind,mypos] 
        if ((thisheight - myposd) > 0):
            arrow = -1
        else:
            arrow = 1
        myheight = thisheight * 6
        print 'mypos ', mypos, ' myposd ', myposd, ' myheight ', myheight
        fakey = numpy.arange(0,xd,xl)
        x=pr[iselect,mind,:]
        n, bins, patches = ax.hist(x, bins=xbins, cumulative=True, histtype='step', facecolor='red', alpha=0.5)
        l = plt.plot(bins, fakey, 'r--', linewidth=1)
        if (arrow == 1):
            ax.annotate('D', xy=(myposd, myheight),  xycoords='data',
                        xytext=(myposd, 0), textcoords='axes fraction',
                        arrowprops=dict(facecolor='black', shrink=0.00),
                        horizontalalignment='center', verticalalignment='bottom',)
        if (arrow == -1):
            ax.annotate('D', xy=(myposd, myheight),  xycoords='data',
                        xytext=(myposd, 0), textcoords='axes fraction',
                        arrowprops=dict(facecolor='black', shrink=0.00),
                        horizontalalignment='center', verticalalignment='bottom',)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, xs)


                        #cumdist[chain,mind,j] = cumdist[chain,mind,j] + 1.0/float(numbins);
                 #   cumdistlength[chain,mind] = thisdist
                #    cumdistposition[chain,mind] = j


    plotf='%s_%03i_cum_probs.png' % (filestem,iselect)
    plt.savefig(plotf)
    print 'Look at %s' % plotf


    return

#-------------------------------------------------------------------

def plotAllProb(totalprEE,totalprBB,patch_name,filestem,addthis):

    """
    Plot probs to file for selected chain
    """

    import numpy
    import matplotlib.pyplot as plt

    # Plot bin histograms
    fig=plt.figure(figsize=(8.5,11),dpi=150) # Letter size
    fig.suptitle('Patch %s, all_jks'%patch_name)

    x1=totalprEE[:]
    x2=totalprBB[:]
    y1=numpy.zeros(numpy.shape(x1))
    y2=numpy.zeros(numpy.shape(x2))

    ax = fig.add_subplot(2,1,1)
    ax.set_xlabel('Prob(mu_obs > mu_sims) EE')
    ax.plot(x1,y1,linestyle='None', marker='x', markersize=10)
    ax.set_xlim((0,1))
    ax.set_ylim((-0.05,0.05))
    ax = fig.add_subplot(2,1,2)
    ax.set_ylim((-0.05,0.05))
    ax.set_xlim(0, 1)
    ax.set_xlabel('Prob(mu_obs > mu_sims) EE')
    ax.plot(x2,y2,linestyle='None', marker='x', markersize=10)

    plotf='%s_%03i_all_probs.png' % ('alljks',(0+addthis))
    plt.savefig(plotf)
    print 'Look at %s' % plotf


    fig2=plt.figure(figsize=(8.5,11),dpi=150) # Letter size
    fig2.suptitle('Patch %s, all_jks'%patch_name)

    ax = fig2.add_subplot(2,1,1)
    ax.set_xlabel('Cum Prob(mu_obs > mu_sims) EE')
    xbins = numpy.arange(0,1,0.01)
    xd=x1.shape[0]
    xs=xd +1
    xl=xd/100.0
    fakey = numpy.arange(0,xd,xl)
    n, bins, patches = ax.hist(x1, bins=xbins, cumulative=True, histtype='step', facecolor='red', alpha=0.5)
    l = plt.plot(bins, fakey, 'r--', linewidth=1)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, xs)
    ax = fig2.add_subplot(2,1,2)
    ax.set_xlabel('Cum Prob(mu_obs > mu_sims) BB')
    n, bins, patches = ax.hist(x2, bins=xbins, cumulative=True, histtype='step', facecolor='red', alpha=0.5)
    l = plt.plot(bins, fakey, 'r--', linewidth=1)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, xs)

    plotf='%s_%03i_all_cum_probs.png' % ('alljks',(0+addthis))
    plt.savefig(plotf)
    print 'Look at %s' % plotf

    fig3=plt.figure(figsize=(8.5,11),dpi=150) # Letter size
    fig3.suptitle('Patch %s, all_jks'%patch_name)

    ax = fig3.add_subplot(2,1,1)
    ax.set_xlabel('Cum Prob(mu_obs > mu_sims) EE')
    xbins = numpy.arange(0,1,0.01)
    xd=x1.shape[0]
    xs=xd +1
    xl=xd/100.0
    fakey = numpy.arange(0,xd,xl)
    n, bins, patches = ax.hist(x1, bins=xbins, cumulative=False, histtype='step', facecolor='red', alpha=0.5)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 10)
    ax = fig3.add_subplot(2,1,2)
    ax.set_xlabel('Cum Prob(mu_obs > mu_sims) BB')
    n, bins, patches = ax.hist(x2, bins=xbins, cumulative=False, histtype='step', facecolor='red', alpha=0.5)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 10)

    plotf='%s_%03i_all_non_cum_probs.png' % ('alljks',(0+addthis))
    plt.savefig(plotf)
    print 'Look at %s' % plotf

    return

#-------------------------------------------------------------------


def plotChisq(chisq,idata):

    """
    Plot histogram of chisq to file
    """

    import numpy
    import matplotlib.pyplot as plt
    from scipy import stats

    # Plot bin histograms
    fig=plt.figure(figsize=(11,11),dpi=150) # Letter size
    fig.suptitle('')
    xbins = numpy.arange(0,12,0.2)

    weights=1.0e-2*numpy.ones(numpy.shape(chisq))

    for imode,mode in enumerate(modes):
        mind=ModesLookup[mode]
        #print imode+1+len(modes)*ibin,ell,mode
        #print 'Drawing bin %i (ell=%i-%i)' %(ibin,ell-dl/2,ell+dl/2)
        ax = fig.add_subplot(len(modes),1,imode+1)
        nsim, bins, dummy = ax.hist(chisq[:,idata+1:,mind,:].flatten(),bins=xbins,color='blue',alpha=0.5,weights=weights[:,idata+1:,mind,:].flatten(),normed=False)
        #n, bins, dummy = ax.hist(chisq[idata:,mind,:], bins=xbins)
        n, bins, dummy = ax.hist(chisq[:,idata,mind,:].flatten(),bins=xbins,color='red',alpha=0.5,normed=False)
        #print nsim,sum(nsim)
        #print n,sum(n)
        print sum(n), sum(nsim)
        #assert(sum(n)==chisq[:,idata,mind,:].flatten().size), 'size mismatch!'
        ax.set_xlabel('chisq (%s)' % mode)
        plt.ylim(0,14)

        #for k in range(10):
        #    curves=stats.chi2.pdf(xbins,k)
        #    ax.plot(xbins,curves,'r-',linewidth=1,color='black')

    plotf='chisq.png'
    plt.savefig(plotf)
    print 'Look at %s' % plotf

    return

#-------------------------------------------------------------------

def xfrange(start, stop=None, step=None):

    """
    Like range(), but returns list of floats instead
    All numbers are generated on-demand using generators
    from http://code.activestate.com/recipes/66472-frange-a-range-function-with-float-increments
    """

    if stop is None:
        stop = float(start)
        start = 0.0
        
    if step is None:
        step = 1.0
            
    cur = float(start)
            
    while cur < stop:
        yield cur
        cur += step

    return

#-------------------------------------------------------------------

def dumpKS(dmin=0,dmax=1.0,dd=0.001,Nmin=1,Nmax=100,dn=1,f='ks.grid'):

    """
    Dump the KS distribution into a file for plotting
    """

    out=open(f,'w')

    out.write('# d n PKS\n')

    for n in range(Nmin,Nmax+dn,dn):
        for d in xfrange(dmin,dmax+dd,dd):
            out.write('%f %i %f\n'%(d,n,p(d,n)))

    out.close()

    print 'Look in %s' % f

    return

#-------------------------------------------------------------------

def plotKSHistogram(dks,jknames,ijk=''):

    """
    Plot histogram of KS D statistics
    """

    import numpy
    import matplotlib.pyplot as plt
    from scipy import stats

    # Plot bin histograms
    fig=plt.figure(figsize=(11,11),dpi=150) # Letter size
    fig.suptitle('jackknife %s - %s'%(str(ijk),jknames[ijk]))
    xbins = numpy.arange(0,0.7,0.005)
    #print numpy.shape(dks)

    n=numpy.size(dks[:,idata+1:])
    #print n
    weights=numpy.ones(numpy.shape(dks))/float(n)
    #weights=numpy.ones(numpy.shape(dks))
    #eps=numpy.finfo(numpy.float).eps
    assert(weights[:,idata+1:].sum()-1.0<1.0e-6)

    ax = fig.add_subplot(1,1,1)
    nsim, bins, dummy = ax.hist(dks[:,idata+1:].flatten(),bins=xbins,\
                                weights=weights[:,idata+1:].flatten(),\
                                normed=False,color='blue',alpha=1.0,\
                                cumulative=-1,histtype='step',label='sims')
    ax.set_xlabel('$D$')
    ax.set_ylabel('$P_{KS}(D)$')
    ax.set_yticks(numpy.arange(0.0,1.1,0.1))
    ax.set_xticks(numpy.arange(0.0,0.8,0.1))
    ax.grid(color='k',linestyle='-', linewidth=1)

    P_theory=numpy.zeros(numpy.shape(xbins))
    nbins=12
    for i,x in enumerate(xbins):
        P_theory[i]=P(x,nbins)
    ax.plot(xbins,P_theory,color='r',linewidth=2,label='theory ($n=%i$)'%nbins)

    plt.ylim(0.0,1.01)
    plt.legend()

    plotf='ks-ds-%s.png' % ijk
    plt.savefig(plotf)
    print 'Look at %s' % plotf

    return

#-------------------------------------------------------------------

def plotScaledMean(chi,idata):

    """
    Plot histogram of (mean/error) (i.e. ~ sqrt(chisq)) to file
    """

    import numpy
    import matplotlib.pyplot as plt
    from scipy import stats

    # Plot bin histograms
    fig=plt.figure(figsize=(11,11),dpi=150) # Letter size
    fig.suptitle('')
    xbins = numpy.arange(0,4,0.2)

    weights=1.0e-2*numpy.ones(numpy.shape(chi))

    for imode,mode in enumerate(modes):
        mind=ModesLookup[mode]
        ax = fig.add_subplot(len(modes),1,imode+1)
        nsim, bins, dummy = ax.hist(chi[:,idata+1:,mind,:].flatten(),bins=xbins,color='blue',alpha=0.5,weights=weights[:,idata+1:,mind,:].flatten(),normed=False)
        n, bins, dummy = ax.hist(chi[:,idata,mind,:].flatten(),bins=xbins,color='red',alpha=0.5,normed=False)
        #print sum(n), sum(nsim)
        ax.set_xlabel('mean/error (%s)' % mode)
        plt.ylim(0,30)

    plotf='mean-scaled.png'
    plt.savefig(plotf)
    print 'Look at %s' % plotf

    return

#-------------------------------------------------------------------

def doKSTest(values):

    """
    Carry out the KS test for a list of values
    Return d and p_KS
    """

    from scipy import stats

    (d,p)=ks_uniform(values)
    (d_scipy,p_scipy)=stats.kstest(values,'uniform')

    #print 'KS stat = %f (%f scipy)' % (d,d_scipy)
    ## print 'KS prob = %f' % (1.0-p)
    #print 'KS prob = %f (%f scipy)' % (p,p_scipy)

    return d,p

#-------------------------------------------------------------------

def fetchP(ds,d):

    """
    Given a list of D's from simulations (and data), calculate a P value
    as the position of a data-D in that distribution
    This function is designed to return P(D) given a D distribution for sims
    """

    ds_sorted=sorted(ds)

    index=ds_sorted.index(d)
                
    num_sims=len(ds_sorted)-1

    P=float(index)/float(num_sims)

    return 1.0-P

#-------------------------------------------------------------------

    # Read the bins
    #binsf='/data4/quiet/commander_setup/multi_n128/bins_d50.dat'
    #b=open(binsf,'r')
    #b.readline()
    #x,y,z=b.readline().split()
    #print x
    #for i in range(int(x)):
    #    print b.readline()
    #b.close()

#-------------------------------------------------------------------

def pChisq(chisq,n):

    """
    Calculate the analytic PTE for a chisq value with n degrees of freedom
    """

    from scipy import stats

    return 1-stats.chisqprob(chisq,n)

#-------------------------------------------------------------------

def writeChisq(chisq,idata,jkdict,map2cl=True):

    """
    Write chisq stats to file
    """

    reportf='chisq.txt'

    (numjks,numchain)=numpy.shape(chisq)

    # ks-ds.txt should contain P's and D's for sims and for all jacks
    # NB prob_theory should be the same as ks
    k=open(reportf,'w')
    k.write('# ijk jk ichain chisq\n')
    for ichain in range(numchain):
        for ijk in range(numjks):
            if ichain==idata:
                print '%i %s %f' % (ijk,jkdict[ijk],chisq[ijk,ichain])
            line='%i %s %i %5.3f\n' % (ijk,jkdict[ijk],ichain,chisq[ijk,ichain])
            k.write(line)

    print 'Look in %s' % reportf
    k.close()

    return

#-------------------------------------------------------------------

def writeKS(probs,dks,ks,numbins,idata,jkdict,map2cl=False,\
            reportf='ks-ds.txt',texf='ks-ds.tex'):

    """
    Report/write the KS results for each jackknife, for the data chain
    """

    import numpy

    if map2cl:
        dict=ModesLookupMap2Cl
    else:
        dict=ModesLookup
   
    reportf='ks-ds.txt'; texf='ks-ds.tex'
    prob_theory=numpy.zeros(numpy.shape(ks))
    prob_sims=numpy.zeros(numpy.shape(ks))
    (numjks,numchain)=numpy.shape(ks)

    # ks-ds.txt should contain P's and D's for sims and for all jacks
    # NB prob_theory should be the same as ks
    k=open(reportf,'w')
    t=open(texf,'w')
    k.write('# ijk jk ichain d P_KS(d) P_KS(d) P_sims(d)\n')
    t.write('\\begin{table}[h]\n\\begin{tabular}{llll}\nJackknife & $D$ & $P_{\mathrm{theory}}$ & $P_{\mathrm{sims}}$ \\\\\n\hline\n')
    N=numbins*len(modes)
    for ichain in range(numchain):
        for ijk in range(numjks):
            # Use the D against dks for sims
            prob_theory[ijk,ichain]=p(dks[ijk,ichain],N)
            prob_sims[ijk,ichain]=fetchP(dks[ijk,:].flatten().tolist(),dks[ijk,ichain])
            if ichain==idata:
                print '%i %s %i %f %f %f %f' % (ijk,jkdict[ijk],ichain,dks[ijk,ichain],ks[ijk,ichain],prob_theory[ijk,ichain],prob_sims[ijk,ichain])
                tex='%s & %5.3f & %4.2f & %4.2f \\\\\n' % (jkdict[ijk],dks[ijk,ichain],ks[ijk,ichain],prob_sims[ijk,ichain])
                t.write(tex)
            line='%i %s %i %5.3f %4.2f %4.2f %4.2f\n' % (ijk,jkdict[ijk],ichain,dks[ijk,ichain],ks[ijk,ichain],prob_theory[ijk,ichain],prob_sims[ijk,ichain])
            k.write(line)

    # Report KS results for all jackknives taken together
    for ichain in range(numchain):
        EE=probs[:,ichain,dict['EE'],:].flatten().tolist()
        BB=probs[:,ichain,dict['BB'],:].flatten().tolist()
        D_all,p_all=doKSTest(EE+BB)
        if ichain==idata:
            print 'For all jks, D = %f, p = %f' % (D_all,p_all)
            tex='%s & %5.3f & %4.2f & %4.2f \\\\\n' % ('\%all',D_all,p_all,-1.0)
            t.write(tex)
        line='%i %s %i %5.3f %4.2f %4.2f %4.2f\n' % (-1,'all',ichain,D_all,p_all,-1.0,-1.0)
        k.write(line)

    print 'Look in %s (and %s)' % (reportf, texf)
    k.close()
    t.write('\\end{tabular}\n\\caption{K--S test results.\\label{table:xx}}\n\\end{table}\n')
    t.close()

    return

#-------------------------------------------------------------------

if __name__ == '__main__':

    import os,sys,numpy,pyfits
    from optparse import OptionParser

    # Set up options
    usage = '%s -h\n' % sys.argv[0]
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="make lots of noise")

    parser.add_option("-p", "--patch", dest="patch",\
                      help="specify QUIET CMB patch", type="string",\
                      metavar='(2a|4a|6a|7b)')

    parser.add_option("-i", "--idata", dest="idata",\
                      help="specify idata", type="int",\
                      metavar='0')

    parser.add_option("-b", "--burnin", dest="burnin",\
                      help="specify burnin", type="int",\
                      metavar='100')

    parser.add_option("-s", "--summed", action="store_true", dest="dosum",
                      help="analyze summed instead of diff data")

    parser.add_option("-k", "--ks", action="store_true", dest="doks",
                      help="post-analyze K-S D distributions")

    parser.add_option("-d", "--dir", dest="directory",\
                      help="specify directory of spectra", type="string",\
                      metavar='.')

    parser.add_option("-m", "--m", dest="multi",\
                      help="specify base path for multiple jackknives",\
                      type="string",metavar='/path1,/path2')

    parser.add_option("-j", "--jackknife", dest="jackknife",\
                      help="specify qualifying stem for jackknife",\
                      type="string",metavar='.')

    parser.add_option("-n", "--nside", dest="nside",
                      help="specify nside", type="int")

    parser.add_option("-c", "--nchains", dest="numchains_in",
                      help="specify nchains", type="int")

    parser.add_option("-x", "--version", dest="version",
                      help="cls.dat version to use", type="int")

    # Assert defaults
    parser.set_defaults(verbose=True,patch=None,idata=0,burnin=100,\
                        dosum=False,doks=False,directory=None,\
                        jackknife='jk',multi=None,version=0,nside=None,
                        numchains_in=None)

    # Parse command line
    (options, args) = parser.parse_args()

    verbose=options.verbose
    patch=options.patch
    # Index of data (as opposed to sims)
    idata=options.idata
    # Set burnin estimate
    burnin=options.burnin
    dosum=options.dosum
    doks=options.doks
    directory=options.directory
    jackknife=options.jackknife
    multi=options.multi
    version=options.version
    nside=options.nside
    numchains_in=options.numchains_in

    assert(nside==128 or nside==256), 'disallowed nside (%i)' % nside

    assert(patch is not None), 'patch required!'

    # Avoid looking at the data - to fix
    #assert(dosum!=(idata==0)), "You can't do that.."

    if multi is not None:
        # Fetch the list of jackknives and establish parameters
        dirs,jkdict=scanDirs(multi,nside=nside)
        numjks=len(dirs)
        numspec=len(modes)
        cls[ijk,:]=readPowspecs(directory,version=version,nside=nside,\
                                override_chains=1)
        #(numchains,numbins,numcols)=numpy.shape(readPowspecs(dirs[0],\
        #                                        version=version,\
        #                                        nside=nside))
        #cls=readPowspecs(dirs[0],version=version,nside=nside)
        if numchains_in is not None: numchains=numchains_in
        print 'numjks    = %i' % numjks
        print 'numchains = %i' % numchains
        print 'numbins   = %i' % numbins
        print 'numcols   = %i' % numcols
        print 'numspec   = %i' % numspec

        # Initialize arrays
        cls=numpy.zeros((numjks,numchains,numbins,numcols))
        chi=numpy.zeros((numjks,numchains,numbins,numspec))
        chisq=numpy.zeros((numjks,numchains,numbins,numspec))
        pr=numpy.zeros((numjks,numchains,numspec,numbins))
        dks=numpy.zeros((numjks,numchains))
        pks=numpy.zeros((numjks,numchains))
        pks_sims=numpy.zeros((numjks,numchains))
        # Now load the data
        print 'Reading power spectra...'
        # For each jackknife
        for ijk,directory in enumerate(dirs):
            print '\nReading %i %s' % (ijk,jkdict[ijk])
            cls[ijk,:]=readPowspecs(directory,version=version,nside=nside,\
                                    override_chains=numchains)
            #if dosum:
            #    print '\nCalculating and writing mean power spectrum...'
            #    writeMeanSpectrum(cls[ijk,:],idata,jkdict[ijk])
            #    sys.exit(0)
            print '\nDeriving some statistics from spectra for %s...' \
                  % jkdict[ijk]
            chisq[ijk,:],chi[ijk,:],pr[ijk,:],dks[ijk,:],\
                 pks[ijk,:],pks_sims[ijk,:]=\
                 calculateSpectraStats(cls[ijk,:],idata)
            print '\nPlotting chisq distribution for %s...' % jkdict[ijk]
            plotChisqs(chisq[ijk,:],idata,jkdict[ijk])

            print '\nConverting and plotting power spectra...'
            data={}
            data['ell'],data['EE'],data['EEerror'],data['BB'],data['BBerror'],\
                data['EEsigmaaway'],data['BBsigmaaway'], data['minmaxvalues']\
                =reducePowerSpec(cls[ijk,:],7)
            #for ichain in [idata,1,2,numchains]:
            for ichain in [idata]:
                plotPowerSpecMap2Cl(data,patch,jkdict[ijk],ichain)

        # Now we have all the jackknives,
        print '\nWriting KS stats to file for all jks...'
        writeKS(pr,dks,pks,numbins,idata,jkdict,map2cl=True)

        # Now we have all the jackknives,
        print '\nWriting chisq stats to file for all jks...'
        writeChisq(chisq,idata,jkdict,map2cl=True)

        sys.exit(0)


    if directory is not None:
        print '\nReading power spectra...'
        cls=readPowspecs(directory,version=version)
        (numchains,numbins,numcols)=numpy.shape(cls)
        print 'numchains = %i' % numchains
        print 'numbins   = %i' % numbins
        print 'numcols   = %i' % numcols

        numchains=10

        print '\nDeriving some statistics from spectra...'
        chisq,chi,pr,dks,pks,pks_sims=calculateSpectraStats(cls,idata)

        print '\nWriting KS stats to file...'
        writeKS(pr[numpy.newaxis,:],dks[numpy.newaxis,:],\
                pks[numpy.newaxis,:],numbins,idata,jkdict={0:jackknife},\
                map2cl=True)

        print '\nPlotting chisq distribution...'
        plotChisqs(chisq,idata,jackknife)

        print '\nConverting and plotting power spectra...'
        data={}
        #data['ell'],data['EE'],data['EEerror'],data['BB'],data['BBerror'],\
        #        data['EEsigmaaway'],data['BBsigmaaway'], data['minmaxvalues']\
        #        =reducePowerSpec(cls,7)
        #for ichain in [idata,1,2,numchains]:
        #    plotPowerSpecMap2Cl(data,patch,jackknife,ichain)

        sys.exit(0)

    # A postprocessing step for KS stats
    if doks:
        numjks=12
        numchain=101
        ksf='ks-ds.txt'
        dks,jknames=readKSDS(ksf,numjks,numchain)
        print numpy.shape(dks)
        # Plot P(D) histograms
        plotKSHistogram(dks,jknames,ijk='all')
        for ijk in range(numjks):
            # newaxis is needed to replace the slice-lost jk axis
            plotKSHistogram(dks[numpy.newaxis,ijk,:],jknames,ijk=str(ijk))
        sys.exit(0)

    # Specify bins
    numbins=6
    lmin=50; lmax=300; dl=50

    # Fetch these parameters from first FITS file (read header only)
    first=pyfits.open(args[0])
    header=first[0].header
    numchain=header['NUMCHAIN']
    numspec=header['NUMSPEC']
    numsamp=header['NUMSAMP']
    del first,header

    patch_name=patch
    numjks=len(args)

    if not dosum:
        probs=numpy.zeros((numjks,numchain,numspec,numbins))
        chi=numpy.zeros((numjks,numchain,numspec,numbins))
        poorchi=numpy.zeros((numjks,numchain,numspec,numbins))
        ks=numpy.zeros((numjks,numchain))
        dks=numpy.zeros((numjks,numchain))

        # Rob to clean these up
        # (Suggest dict - totalpr['EE2'], etc, or just index)
        totalprEE=numpy.zeros(numjks*numbins)
        totalprBB=numpy.zeros(numjks*numbins)
        totalprEE2=numpy.zeros(numjks*numbins)
        totalprBB2=numpy.zeros(numjks*numbins)
        totalprEE3=numpy.zeros(numjks*numbins)
        totalprBB3=numpy.zeros(numjks*numbins)
        burns=numpy.zeros((numjks,numchain,numspec,numbins,numsamp))

    jkdict={}
    for ijk,f in enumerate(args):
        data,lmax_dummy,numsamp,numchain,numspec,fstem=readCLs(f)
        jkdict[ijk]=fstem
        mu,p,pr,std,map,median,clower,cupper,error,chisq,chisq2,meanscal=\
            calculateStats(data,lmin,lmax,dl,burnin,numchain,\
                           numspec,numbins,fstem,dosum)
        # For summed power spectra:
        if dosum:
            for bin in range(numbins):
                mEE=numpy.mean(map[idata+1:,ModesLookup['EE'],bin])
                sEE=numpy.std(map[idata+1:,ModesLookup['EE'],bin])
                mBB=numpy.std(map[idata+1:,ModesLookup['BB'],bin])
                sBB=numpy.std(map[idata+1:,ModesLookup['BB'],bin])
                print '%i %f %f %f %f'%(bin,mEE,sEE,mBB,sBB)
            for ichain in range(numchain):
                if ichain==idata: continue
                plotHistogram(data,lmin,lmax,dl,burnin,numchain,numbins,\
                              numspec,mu,median,map,clower,cupper,fstem,ichain)
            plotSummedPowspec(data,lmin,lmax,dl,burnin,numchain,numbins,\
                              numspec,map,patch_name,fstem,idata)
            # Get out of here...
            sys.exit(0)

        # Calculate cumulative probabilities
        cumdist,cumdistlength,cumdistposition=calculateCumStats(pr)
        # ??
        burns,noburns=calculateBurnin(data,lmin,lmax,dl,burnin,numchain,\
                                      numspec,numbins,fstem)
        # Plot probs for data and one sim
        for ipp in range(idata+2):
            plotProb(pr,cumdist,cumdistlength,cumdistposition,patch_name,\
                     fstem,ipp)
        # Plot burns for data and two sims
        for ipb in range(idata+3):
            plotBurns(noburns,burns,lmin,lmax,dl,numchain,numbins,numspec,\
                      patch_name,fstem,ipb)

        probs[ijk,:,:,:]=pr[numpy.newaxis,:,:,:]
        chi[ijk,:,:,:]=chisq[numpy.newaxis,:,:,:]
        poorchi[ijk,:,:,:]=meanscal[numpy.newaxis,:,:,:]

        for ichain in range(numchain):
            # For both modes, and all bins
            EE=probs[ijk,ichain,ModesLookup['EE'],:].flatten().tolist()
            BB=probs[ijk,ichain,ModesLookup['BB'],:].flatten().tolist()
            dks[ijk,ichain],ks[ijk,ichain]=doKSTest(EE+BB)
            if ichain==idata:
                print 'For chain %i (jk %s),' % (ichain,jkdict[ijk])
                print EE+BB
                print 'KS stat = %f' % dks[ijk,ichain]
                print 'KS prob = %f' % ks[ijk,ichain]

        for idd in [idata,idata+1]:
            # Only plot samples histogram for data and first simulation
            plotHistogram(data,lmin,lmax,dl,burnin,numchain,numbins,numspec,\
                          mu,median,map,clower,cupper,fstem,idd)
            # Plot power spectrum for data and first simulation
            plotPowerSpec(data,lmin,lmax,dl,burnin,numchain,numbins,numspec,\
                          mu,median,clower,cupper,patch_name,fstem,idd)

    # Not sure what this loop is, but Rob is going to tidy...
    for ijk in range(numjks):
        for addto in range(numbins):
            totalprEE[((ijk * numbins) + addto)]=probs[ijk,idata,ModesLookup['EE'],addto]
            if (totalprEE[((ijk * numbins) + addto)] > 0.99):
                totalprEE[((ijk * numbins) + addto)]=0.99
            if (totalprEE[((ijk * numbins) + addto)] < 0.01):
                totalprEE[((ijk * numbins) + addto)]=0.01
            totalprBB[((ijk * numbins) + addto)]=probs[ijk,idata,ModesLookup['BB'],addto]
            if (totalprBB[((ijk * numbins) + addto)] > 0.99):
                totalprBB[((ijk * numbins) + addto)]=0.99
            if (totalprBB[((ijk * numbins) + addto)] < 0.01):
                totalprBB[((ijk * numbins) + addto)]=0.01
            totalprEE2[((ijk * numbins) + addto)]=probs[ijk,idata+1,ModesLookup['EE'],addto]
            if (totalprEE2[((ijk * numbins) + addto)] > 0.99):
                totalprEE2[((ijk * numbins) + addto)]=0.99
            if (totalprEE2[((ijk * numbins) + addto)] < 0.01):
                totalprEE2[((ijk * numbins) + addto)]=0.01
            totalprBB2[((ijk * numbins) + addto)]=probs[ijk,idata+1,ModesLookup['BB'],addto]
            if (totalprBB2[((ijk * numbins) + addto)] > 0.99):
                totalprBB2[((ijk * numbins) + addto)]=0.99
            if (totalprBB2[((ijk * numbins) + addto)] < 0.01):
                totalprBB2[((ijk * numbins) + addto)]=0.01
            totalprEE3[((ijk * numbins) + addto)]=probs[ijk,idata+2,ModesLookup['EE'],addto]
            if (totalprEE3[((ijk * numbins) + addto)] > 0.99):
                totalprEE3[((ijk * numbins) + addto)]=0.99
            if (totalprEE3[((ijk * numbins) + addto)] < 0.01):
                totalprEE3[((ijk * numbins) + addto)]=0.01
            totalprBB3[((ijk * numbins) + addto)]=probs[ijk,idata+2,ModesLookup['BB'],addto]
            if (totalprBB3[((ijk * numbins) + addto)] > 0.99):
                totalprBB3[((ijk * numbins) + addto)]=0.99
            if (totalprBB3[((ijk * numbins) + addto)] < 0.01):
                totalprBB3[((ijk * numbins) + addto)]=0.01

    # Now analyze/report all jackknives together
    plotAllProb(totalprEE,totalprBB,patch_name,fstem,idata)
    plotAllProb(totalprEE2,totalprBB2,patch_name,fstem,idata+1)
    plotAllProb(totalprEE3,totalprBB3,patch_name,fstem,idata+2)

    #plotChisq(chi,idata)
    # Plot the poor man's chisq distribution
    plotScaledMean(poorchi,idata)

    # Calculate and write out the KS stats to file
    writeKS(probs,dks,ks,numbins,idata,jkdict)

    sys.exit(0)
