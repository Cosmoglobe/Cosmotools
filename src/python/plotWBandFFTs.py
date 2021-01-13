#!/usr/bin/env python

"""
Plot FFT power spectrum for all modules of a CES.
pngs are colour-coded by status: 0 = blue; contains 32 = red; others = black
"""

import matplotlib
matplotlib.use("Agg")

#-------------------------------------------------------------------

# Define ordering of diodes
Diodes=['Q1','U1','U2','Q2']
DiodesLookup={'Q1':0,'U1':1,'U2':2,'Q2':3}
Q1='Q1'; U1='U1'; U2='U2'; Q2='Q2'
NDIODES=len(Diodes)
NMODULES=81

#-------------------------------------------------------------------

def _cut_data(chain,n=50):
    m=len(chain)/n*n
    return chain[:m]

#-------------------------------------------------------------------

def calcPowerSpectrum(tod,fsamp,fit_length=5000,cut_multiple=50):

    """
    Fit a 1/f + white noise power spectrum to the power spectrum data,
    of the form
    P(f) = sigma^2 (1+(fknee/f)^alpha )
    The arguments are the TOD itself and the sampling frequency.
    Optional keywords are fit_length, a first iteration of the range of
    f to fit to (default 5000)
    and cut_multiple, which specifies how the chain is to be cut to
    improve speed (default 50)

    If this runs too slowly you could pass tod[::n],fsamp/n as the
    arguments for some thinning fraction n.

    The fit uses the algorithm described in astro-ph/0405462
    In particular you should not just fit least squares to the power
    spectrum in real space - 
    you have to do it in log space and include a constant offset.
    Just read the paper.
    """

    from numpy import array, pi, arange, log, where, exp
    from numpy.fft import fft, fftfreq

    #Cut the chain down so it does not a large prime length (speeds things up)
    #If this runs too slowly it might be worth considering passing
    #a thinned chain.
    chain=_cut_data(tod,cut_multiple)

    #Compute the spacing and sizes off the Fourier space array
    f=fftfreq(len(chain),1.0/fsamp)
    df = f[1]-f[0]
    n = len(chain) / 2 - 1

    #Construct the x-axis of the power spectrum
    j = slice(1,n+1)
    #k = f[j]

    #Compute the power spectrum.
    Pj = (abs(fft(chain))**2)[j]

    # Rescale (a la Newburgh)
    #Pj /= numpy.sqrt(float(n))**2
    Pj /= float(n)

    #We fit to the log of the power spectrum,
    #and only fit to the initial part of it
    #so as not to weight the white noise tail to much. 
    log_Pj = log(Pj[:fit_length])

    # Calculate max frequency
    max_freq=df*n

    return log_Pj,max_freq,df

#-------------------------------------------------------------------

def loadCESStats(ces_stats='./ces_statistics.txt'):

    """
    Read a CES statsfile into memory
    """

    indata=open(ces_stats,'r')
    accept={}; powspec={}; chisq={}
    print 'Reading %s (for the statuses etc.)...' % ces_stats
    for line in indata:
        if line.startswith('#'): continue
        try:
            cesid,obj,mod,di,status,r,sigma0,alpha,fknee,\
               chisq_025,sscan,chisq_002,CMB,chisq_1025,\
               TODchisq,TODabs,chi_az,typeb,wthr10,wther30,jump\
               =line.split()
        except:
            print 'Skipping %s' % line
            pass

        #accept[(int(cesid),int(mod),int(di))]=int(status)
        powspec[(int(cesid),int(mod),int(di))]=(sigma0,alpha,fknee)
        #chisq[(int(cesid),int(mod),int(di))]=(chisq_025,sscan,chisq_002,\
        #                                      CMB,chisq_1025)

    indata.close()

    stats=(accept,powspec,chisq)

    print '%s read OK' % ces_stats

    return stats

#-------------------------------------------------------------------

def queryStatus(cut,status):

    """
    Ask whether cut has been applied by querying the status word
    """

    return cut == cut & status

#------------------------------------------------------------------- 

def plotSingleModuleFFT(powspecs,samples,module,ces,status,curve,mjdrange):

    """
    For a single module, plot the FFTs
    Later add best-fit noise models
    """

    import matplotlib.pyplot as plt
    from math import log
    import numpy

    plt.figure(1)

    # It would be nice to have the MJD in the title
    #but these are less available
    title='ML%i RW%02i %f-%f' % (ces,module,mjdrange[0],mjdrange[-1])

    #title='*%i* - %s ML%i RW%02i di%i %03i %s - %f-%f' \
    #       % (accept[(module,di%4)],target[-2:],ces,module,di%4,di,\
    #          Diodes[di%4],mjdstart,mjdstop)

    for idiode in range(NDIODES):
        plt.subplot(2,2,idiode+1)

        #print curve[(mod,idiode)]

        if status[(mod,idiode)] == 0:
            colour='blue'
        else:
            colour='black'
            if queryStatus(32,status[(mod,idiode)]):
                colour='red'

        plt.plot(samples[idiode],powspecs[idiode],'.',\
                 markersize=0.6,color=colour)

        #print mod,idiode
        (sigma0,alpha,fknee)=curve[(mod,idiode)]
        fitted_curve=[]
        flag=None
        for sample,n in enumerate(samples[idiode]):
            try:
                # I'm not sure about a factor of (n/)2 here,
                # but it seems to work empirically
                fitted_curve.append(log((sigma0**2)*(1.0+(fknee/float(n))**-alpha)))
            except ValueError:
                print sigma0,alpha,fknee
                #fitted_curve.append(log(sigma0**2))
                #print 'Ignoring ValueError in hockey stick'
                flag='%i %i: %f %f %f %i' % (mod,idiode,sigma0,fknee,-alpha,n)
                break
        if flag is not None: print flag
        #assert(len(samples[idiode])==len(fitted_curve)), \
        #  'array lengths: %i != %i' % (len(samples[idiode]),len(fitted_curve))
        if not (sigma0 < numpy.finfo(numpy.float).eps):
            plt.plot(samples[idiode],fitted_curve,color='green')
        else:
            print 'sigma0 = %f -- ignoring curve' % sigma0

        plt.xscale('log')
        plt.xlim(1e-4,max(samples[0]))

        if idiode >=2:
            plt.xlabel('Frequency/Hz')
        if idiode%2 == 0:
            plt.ylabel('log10(Power)')
        if idiode == 0:
            plt.title(title)

        [xmin,xmax,ymin,ymax]=plt.axis()
        xpos=5.0e-4
        ypos=ymin+0.1*(ymax-ymin)
        plt.text(xpos,ypos,'%s - di%i *%i*' % \
                 (Diodes[idiode],idiode,status[module,idiode]))

    #title='*%i* - %s ML%i RW%02i di%i %03i %s - %f-%f' \
    #       % (accept[(module,di%4)],target[-2:],ces,module,di%4,di,\
    #          Diodes[di%4],mjdstart,mjdstop)



    fign='ML%04i_RW%02i_fft.png' % (ces,mod)

    # Output plot and tidy up
    plt.savefig(fign)
    plt.clf()
    print '-> Look in %s' % fign
    plt.close('all')

    return

#-------------------------------------------------------------------

if __name__ == '__main__':

    import sys,numpy
    import h5py,quietW

    from optparse import OptionParser

    print
    # Set up options
    usage = "./plotWBandFFTs.py -h\n"
    parser = OptionParser(usage=usage)

    parser.add_option("-c", "--ces", dest="ces",\
                      help="specify ces", type="int",\
                      metavar='412')

    parser.add_option("-i", "--infile", dest="statsfile",\
                      help="specify input ces_stats file", type="string",\
                      metavar='ces_statistics.txt')

    parser.set_defaults(ces=None,statsfile=None)

    # Parse command line
    (options, args) = parser.parse_args()

    ces=options.ces
    statsfile=options.statsfile

    dostats=True
    if dostats:
        if statsfile is None:
            statsfile='./ces_statistics.txt'

        indata=open(statsfile,'r')
        accept={}; curve={}
        print 'Reading %s (for the statuses)...' % statsfile
        for line in indata:
            if line.startswith('#'): continue
            try:
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
                curve[(int(mod),int(di))]=\
                  (float(sigma0),float(alpha),float(fknee))
        indata.close()

        #stats=loadCESStats(statsfile)
        # accept and chisq are currently empty
        #(accept,powspec,chisq)=stats

    #fields=['time','tod','tp']
    fields=['time','tod']
    data=quietW.l3_read(ces=ces,fields=fields)

    tod=data.tod
    mjdrange=(data.time[0],data.time[-1])

    #print numpy.shape(tod)

    # Set sampling frequency/Hz
    fsamp=50

    for mod in range(NMODULES):
        powspecs=[]; samples=[];
        for di in range(NDIODES):
            idiode=di+4*mod
            powspec,numax,dnu=calcPowerSpectrum(tod[idiode,:],fsamp,\
                                                fit_length=-1,cut_multiple=50)
            #print len(powspec)
            powspecs.append(powspec)
            sample=[(isamp+1)*dnu for isamp in range(len(powspec))]
            samples.append(sample)

        plotFFT=True
        if plotFFT:
            plotSingleModuleFFT(powspecs,samples,mod,\
                                ces,accept,curve,mjdrange)

    dumpfft=False
    if dumpfft:
        for isamp,sample in enumerate(powspec):
            print isamp+1,(isamp+1)*dnu,sample


    # These seem to be all the same, so only print once
    print '\nFrequency/Hz step %f -> %f ' % (dnu,numax)
    print '\nFinished!\n'

    sys.exit(0)
