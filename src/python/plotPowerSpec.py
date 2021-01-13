#!/usr/bin/python

"""
Script to plot single power spectrum from one PCL output C_ell file
Also writes the Chi^2 and bin Chi^2's to file
$Author: laura $
$Date: 2009-11-19 22:27:42 +0100 (Thu, 19 Nov 2009) $
Usage:
      developed on python 2.4-2.5
      ./plotPowerSpec.py clfile\n
      OR ./plotPowerSpec.py args\n
"""

# -------------------------------------------------------------

if __name__ == '__main__':

    import sys

    narg=2
    narg2=5
    if len(sys.argv) != narg+1 and len(sys.argv) != narg2+1:
        print >> sys.stderr, '%s <cl.dat> <clerrb.dat> [mode]' % sys.argv[0]
        print '    output to CWD for now'
        print '    (optional) mode = [EE|BB|both]'
        print 'OR ARGS: <date> <id> <patchnum> <flavour> <mode>'
        sys.exit(2)

    import os,plotutils,secateur

    micro=True # Plot power spectra in \mu K rather than K
    #ellsq=False
    ellsq=True # Whether to calculate/plot ell(ell+1)/2\pi x C_ell, or not
    root='/data4/quiet/nulltests/pcl' # Input rootdir for power spectra
    outdir='.' # Output dir for plots and chisq files
    powerspecerrorsuffix='diff_pixels.txt.cls.dat'
    powerspecdatasuffix='diff_pixels.txt.clest.dat'

    mode=None # Initialize: process EE, BB or both (None -> both)

    # Set up files etc
    if len(sys.argv)>narg2+1: # KLUGE TO MAKE IT RUN!!!!
        date,id,patch,flavour,mode=sys.argv[1:]
        clerrfile=os.path.join(root,'pcl-'+date+'-'+id+'/powspec/pcl-'+date+'-'+id+'-patch_'+patch+'_'+flavour+'.'+powspecerrorsuffix)
        clfile=os.path.join(root,'pcl-'+date+'-'+id+'/powspec/pcl-'+date+'-'+id+'-patch_'+patch+'_'+flavour+'.'+powspecdatasuffix)
    else:
        clfile=sys.argv[1]
        clerrfile=sys.argv[2]
        date,id,patch,flavour=plotutils.DeconstructClFilename(clfile)
        subdir=os.path.join('pcl-'+date+'-'+id,'powspec')
        subdir=os.path.join(root,subdir)
        clfile=os.path.join(subdir,clfile)
        chisqfile='%s-%s-patch_%s.chisq' % (date,id,patch)
        binschisqfile='%s-%s-patch_%s.binsq' % (date,id,patch)
        chisqfile=os.path.join(subdir,chisqfile)
        binschisqfile=os.path.join(subdir,binschisqfile)

    modes=[]
    if mode=='all' or mode=='both' or mode is None:
        modes.extend(['EE','BB'])
    else:
        modes.append(mode)

    # Read power spectrum but not correlation matrix from file
    ps,corr=plotutils.ReadClFile(clfile)
    # Read noise power spectrum and correlation matrix from file
    noiseps,noisecorr=plotutils.ReadClErrFile(clerrfile)

    chi2={}
    chi2red={}
    binschi2={}
    binschi2red={}
    ndof={}

    # Calculate Chi^2
    chisq = False # Laura protection switch
    for mode in ps.keys():
        if chisq:
            chi2[mode],binschi2[mode],ndof[mode],chi2red[mode]= \
                plotutils.CalculateNullChiSquared(ps[mode],ellsq)
            print 'Chi^2 %s %s %f / %d = %f' \
                % (flavour, mode, chi2[mode], ndof[mode], chi2red[mode])
        else:
            chi2[mode]=-99.0
            binschi2[mode]=-99.0
            ndof[mode]=-99.0
            chi2red[mode]=-99.0
            

    # Generate the plot and write Chi^2 to file
    for mode in modes:
        plotfile=secateur.StripExtension(os.path.basename(clfile))[0]\
                  +'-'+mode+'.png'
        data,outfile=plotutils.PlotPowerSpectrum(ps,noiseps,chi2,chi2red,\
                  plotfile,outdir,micro,mode,ellsq,ndof[mode])
        print 'Now look at', outfile
        if chisq:
            plotutils.WriteChiSquared(date,id,patch,flavour,mode,\
                  chi2,chi2red,ellsq,ndof[mode],chisqfile)
            plotutils.WriteBinsChiSquared(date,id,patch,flavour,mode,\
                  binschi2,ellsq,binschisqfile)
        print 'Updated %s.chisq .binsq' \
              % (secateur.StripExtension(os.path.basename(chisqfile))[0])

    sys.exit(0)

# -------------------------------------------------------------
