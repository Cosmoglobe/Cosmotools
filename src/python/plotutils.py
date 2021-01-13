#!/usr/bin/python

"""
Collection of utilities for making null tests, including directory structure \n
$Author: laura $
$Date: 2009-11-19 22:27:42 +0100 (Thu, 19 Nov 2009) $
Usage:
      developed on python 2.4-2.5
      ./plotutils.py ?? ?? \n
Footnotes: \n
"""

import os,sys,re,string,platform,time,datetime
from numpy import *

# -------------------------------------------------------------

def ReadClErrFile(filename='trialcl.dat'):

    """
    Read PCL power spectrum into memory
    Returns a dictionary of lines keyed by field
    """

    import scanf
    from math import pi

    file=open(filename)

    powerspec={'EE':[],'BB':[]}
    corrs=False
    correlations=[]
    lineformat=' %i %i %f %f # %s bandpower %f'
    for line in file:
        if not line.startswith('#') and not corrs: # Ignore comment lines
            lmin,lmax,mean,rms,mode,bias=scanf.sscanf(line,lineformat)
            #lmin,lmax,mean,rms,mark,vec,power,dummy=[el for el in line.split()]
            lc=(lmin+lmax)/2.0
            cl=2.0*pi*mean/(lc*(lc+1.0))
            clrms=2.0*pi*rms/(lc*(lc+1.0))
            powerspec[mode].append([lc,lmin,lmax,cl,clrms,mean,rms,bias])
        if line.startswith('# correlations:'):
            corrs=True
            n=0
            # These are (n_EE+n_BB) * ?(n_EE+n_BB)
            for mode in powerspec:
                n+=len(powerspec[mode])
        if corrs:
            correlations.append([el for el in line.split()])

    file.close()

    return powerspec, correlations

# -------------------------------------------------------------

def ReadClFile(filename='trialcl.dat'):

    """
    Read PCL power spectrum into memory
    Returns a dictionary of lines keyed by field
    """

    import scanf
    from math import pi

    file=open(filename)

    powerspec={'EE':[],'BB':[]}
    lineformat=' %i %i %f # %s bandpower'
    for line in file:
        if not line.startswith('#'): # Ignore comment lines
            lmin,lmax,power,mode=scanf.sscanf(line,lineformat)
            #lmin,lmax,power,mark,vec,dummy=[el for el in line.split()]
            lc=(lmin+lmax)/2.0
            cl=(2.0*pi*power)/(lc*(lc+1.0))
            powerspec[mode].append([lc,lmin,lmax,cl,power])

    file.close()

    correlations=None

    return powerspec, correlations

# -------------------------------------------------------------


def CalculateNullChiSquared(powerspectrum,ellsq=False,reduce=True):

    """
    Calculate the Chi^2 for a power spectrum compared to a null model
    Optionally reduce the ChiSquared and return it too
    """

    # The definition of Chi^2 for points with errors when E=0 is based on
    # Numerical Recipes (edn 3) 15.1.6, p778.

    ChiSquared=0.0
    BinChiSquareds=[]
    for [lc,lmin,lmax,cl,clrms,mean,rms] in powerspectrum:
        if ellsq:
            BinChiSquared = (mean/rms)**2
        else:
            BinChiSquared = (cl/clrms)**2
        ChiSquared += BinChiSquared
        BinChiSquareds.append(BinChiSquared)

    if reduce:
        # What about those horrible degrees of freedom..?
        M=0
        N=len(powerspectrum)
        dof=N-M
        reducedChiSquared = ChiSquared / float(dof)
        return ChiSquared,BinChiSquareds,dof,reducedChiSquared
    else:
        return ChiSquared,BinChiSquareds

# -------------------------------------------------------------

def WriteChiSquared(date,id,patch,flavour,mode,chisq,chisqred, \
                    ellsq=False,ndof=20,filename='chisq.txt'):

    """
    Open Chi-squared file for null suite, if not there already
    """

    import os,secateur

    # Create column header info if creating new file
    if os.path.exists(filename):
        header = None
    else:
        header = '# date runid patch flavour mode ellsq chisq ndof chisqred\n'

    #file=open(filename,'a') # Append or create
    print filename
    file = secateur.OpenFile(filename,'a')

    if header is not None:
        file.write(header)

    line = '%s-%s %s %s %s %s %f %0d %f\n' % (date, id, patch, flavour, mode, \
            str(ellsq), chisq[mode], ndof, chisqred[mode])

    file.write(line)

    file.close()

    return line

# -------------------------------------------------------------

def WriteBinsChiSquared(date,id,patch,flavour,mode, \
                    binschisq,ellsq=False,filename='binsq.txt'):

    """
    Open Bins Chi-squared file for null suite, if not there already
    """

    import secateur

    # Create column header info if creating new file
    if os.path.exists(filename):
        header = None
    else:
        header = '# date runid patch flavour mode ellsq binnum chisq\n'

    #file=open(filename,'a') # Append or create
    file = secateur.OpenFile(filename,'a')

    if header is not None:
        file.write(header)

    for n,bin in enumerate(binschisq[mode]):
        line = '%s-%s %s %s %s %s %d %f\n' % (date, id, patch, flavour,\
                       mode, str(ellsq), n+1, bin)
        file.write(line)

    file.close()

    return line

# -------------------------------------------------------------

def ConstructClFilename(date,id,patch,flavour,extension='diff_pixels.txt.cls.dat'):

    """
    Build Cl filename (including subdir) from constituent parts
    """

    filename = 'pcl-%s-%s/powspec/pcl-%s-%s-patch_%s_%s.%s' % (date,id,date,id,patch,flavour,extension)

    return filename

# -------------------------------------------------------------

def DeconstructClFilename(filename):

    """
    Split C_ell filename into constituent information
    """

    import re


    els = re.split('_',filename)
    dummy,date,id,dummy=re.split('-',els[0])
    patch=els[1]
    flavour=els[2]
    #season=re.split('.',els[2])[0]

    return date,id,patch,flavour

# -------------------------------------------------------------

def PlotPowerSpectrum(powerspec,errpowerspec,chisq,redchisq,\
                      outfile='x.png',path='/data4/quiet/nulltests',\
                      micro=True,mode=None,ellsq=False,ndof=20):

    """
    Plot power spectrum from PCL C_ell file to png file
    You choose whether you want l(l+1)/2pi x C_ell, or not
    """

    import os, re, secateur
    from matplotlib import pyplot as plt
    from matplotlib import ticker as ticker

    if mode is None: mode = 'EE'

    plotdata=[]
    for n,bin in enumerate(powerspec[mode]):
        [lc,lmin,lmax,cl,power]=bin
   	bin2=errpowerspec[mode][n]
        [lc2,lmin2,lmax2,cl2,clrms2,mean2,rms2,bias2]=bin2
        # XX Force noise bias to be zero for now
        mean = 0
	value = power - mean
	errorb = rms2
        if ellsq:
            plotdata.append([lc,value,errorb])
        else:
	    # Currently broken - to be fixed
            print 'THIS DOESNT WORK, CANT UNCORRECT L(L+1)'
            plotdata.append([lc,cl,clrms])

    pdata=array(plotdata)

    ellstart=0.0
    ellstop=400.0
    cllim=1.0e-6
    stretch=5.0e8
    cllim=cllim*stretch

    x=pdata[:,0]
    y=pdata[:,1]
    dy=pdata[:,2]
    xlabel='$\ell$'
    if ellsq:
        elllab='L(L+1)/2\pi' # Can't get the LaTeX to work!
    else:
        elllab='' ####### NB The dollar WAS provided by elllab
    if micro and not ellsq: # Optionally convert K to \mu K
        y*=1.0e12
        dy*=1.0e12
        ylabel='%s$C_{\ell,\rm{%s}}$/$\mu$K$^2$' % (elllab,mode)
        cllim*=1.0e12
    elif ellsq:
        ylabel='%s$C_{\ell,\rm{%s}}$/K$^2$' % (elllab,mode)
        # These two lines are micro hack - to be JZ Revisited
        elllab=''
        ylabel='%s$\ell(\ell+1)/2\pi.C_{\ell,%s}/\mu$K$^2$' % (elllab,mode)
    else:
        ylabel='%s$C_{\ell,\rm{%s}}$/K$^2$' % (elllab,mode)

    # Don't overwrite the plot if it's already there
    if os.path.exists(outfile):
        os.rename(outfile,outfile+'.old')        

    els = re.split('_',outfile)
    dummy,date,run,dummy=re.split('-',els[0])
    patch,season=els[1:3]
    #season=re.split('.',els[2])[0]
    title='patch %s-%s: %s-%s %s' % (patch, mode, date, run, season)
    marker='+'
    ls=''
    if mode=='EE': color='red'
    elif mode=='BB': color='blue'

    def millions(x, pos):
        'The two args are the value and tick position'
        return '$%1.1fM' % (x*1e-6)

    def format1(x, pos):
        'The two args are the value and tick position'
        return '%2.0e' % (x)

    plt.errorbar(x,y,yerr=dy,ls=ls,marker=marker,color=color)
    plt.subplot(111).yaxis.set_major_formatter(formatter=ticker.FuncFormatter(format1)) # Set the y tics explicitly...
    plt.xlim(ellstart,ellstop)
    plt.ylim(-cllim,cllim)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    #plt.axes().set_aspect(aspect='equal')
    #plt.show()

    chisqinfo='%s %s Chi^2: %f / %d = %f' \
               % (mode, str(ellsq), chisq[mode], ndof, redchisq[mode])
    plt.text(0.05*ellstart,cllim,chisqinfo)
    plt.savefig(outfile)
    plt.clf()
    plt.close('all')

    secateur.FixPerms(outfile)
                            
    return plotdata,outfile

# -------------------------------------------------------------

def plot_single_Null(filename1, filename2):

    import pylab, numpy
    import matplotlib
    from matplotlib import pyplot  as plt
    import string

    """
    figure out EEs, put in array
    figure out BBs, put in array
    """

    EE = numpy.loadtxt(filename1)
    #BB = pylab.load(filename2)

    npts = len(EE[:,0])
    print 'read in array'
    els = zeros([npts])
    for i in range (0,npts):
       els[i] = (EE[i,0] + EE[i,1])/2.0
    plt.errorbar(els,EE[:,2],yerr=EE[:,3])
    plt.xlim(0,400)
    plt.ylim(-20,20)
    plt.ylabel('Cl EE, (microK sq)')
    plt.savefig(string.join([filename1,'_plot.png'],''))
    plt.clf()
    plt.close('all')
    print 'closed all, exiting'

    print 'calculating chisq statistic now'
    chisq = 0
    for  i in range (0,npts):
       chisq = chisq + pow((EE[i,2]/EE[i,3]),2)
   
    print chisq, 'Ndof: ', npts


###########################################################################
######### MAIN PART IS HERE ###############################################
###########################################################################

if __name__ =="__main__":
    
    plot_single_Null(sys.argv[1],sys.argv[2]) 


############################################################################
############ END ###########################################################
############################################################################

