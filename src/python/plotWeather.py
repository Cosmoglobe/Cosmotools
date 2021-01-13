#!/usr/bin/python


import os,sys,pylab,numpy,matplotlib,string,math
import datetime,numpy,gaussfit,file_read_utils
from quiet import N_MODULE_MAX_Q,N_MODULE_MAX_W
from matplotlib import pyplot as plt
from matplotlib import ticker as ticker
from numpy import *

user=os.getenv('USER')
if user == 'newburgh': user='%s/repository'%user

#path_to_stats = '/usit/titan/u1/%s/quiet_svn/oslo/src/python/weather_cutting/statisticsfiles/' % user
#path_to_cuts = '/usit/titan/u1/%s/quiet_svn/oslo/src/python/weather_cutting/cuts/' % user
#path_to_plots = '/usit/titan/u1/%s/quiet_svn/oslo/src/python/weather_cutting/plots/' % user
nbins=200

###-----------------------------------------------------------------

def get_avg_in_range(data,module,min,max):

   """
   """

   n=len(data[:,0])
   summer=0.0
   counter=0.0
   output=0.0
   offset=3
   for i in range(n):
      if data[i,module+offset] > min and data[i,module+offset] < max:
         summer=summer+data[i,module+offset]
         counter=counter+1
   if counter !=0:
      output = summer*1.0/counter*1.0
   else:
      output = 0.0

   return output

###-----------------------------------------------------------------

def get_rms_in_range(data,module,min,max):

   """
   """

   n=len(data[:,0])
   summer=0.0
   counter=0.0
   output=0.0
   offset=3
   mean=get_avg_in_range(data,module,min,max)
   for i in range(n):
      if data[i,module+offset] > min and data[i,module+offset] < max:
         summer=summer+(data[i,module+offset]-mean)**2
         counter=counter+1
   if counter !=0:
      output = numpy.sqrt(summer*1.0/counter*1.0)
   else:
      output = 0.0

   return output

####-------------------------------------------------------------

def find_hist_max(data,module,min,max):

   """
   """

   # What is offset?
   offset=3

   # Calculate the histogram
   #print 'nbins %i' % nbins
   histthing=numpy.histogram(data[:,module+offset],nbins,(min,max)) 
   # histthing is a tuple containing
   #    1. the number of histogram entries
   #    2. the bin edges
   n=len(histthing[0])
   print 'm%i (%i histogram entries)' % (module,numpy.sum(histthing[0]))

   # maxval is the peak value of the histogram
   maxval=numpy.max(histthing[0][2:n])
   binthing=0
   counter=0

   # Find (I think) the position, binthing, of the peak
   for i in range(2,n):
      if histthing[0][i] == maxval:
         binthing=histthing[1][i]

   # Count the number of entries to the left of the peak
   for i in range(0,n):
      if histthing[1][i] <= binthing:
	 counter += 1
      histthing[1][i] = histthing[1][i] + 0.5*(histthing[1][1]-histthing[1][0])

   binthing=binthing + 0.5*(histthing[1][1]-histthing[1][0])

   print binthing,counter

   return binthing,counter,histthing,maxval

####-------------------------------------------------------------

def get_histogram(hist):

   """
   Copy histogram hist (which is a tuple of 1-D numpy arrays)
   into a 2-D numpy array newhist
   """

   n=len(hist[0][:])
   newhist=zeros([2,n])
   for i in range(n):
      newhist[0][i] = hist[0][i]
      newhist[1][i] = hist[1][i]

   return newhist

####-------------------------------------------------------------

def get_reflected_histogram(hist,counter,binthing):

   """
   """

   newhist=zeros([2,2*counter+1])
   for i in range(counter):
      newhist[0][i] = hist[0][i]
      newhist[1][i] = hist[1][i]
      newhist[0][2*counter - i] = hist[0][i]
      newhist[1][2*counter - i] = 2*binthing - hist[1][i]
   newhist[0][counter] = hist[0][counter]
   newhist[1][counter] = hist[1][counter]

   return newhist

####-------------------------------------------------------------

def get_reflected_data(peak_loc,count,data,module):

   """
   """

   offset=3
   j=0
   n=len(data[:,module+offset])
   newdata=zeros([n])
   for i in range(n):
      if data[i,module+offset] <= peak_loc:
         newdata[j] = data[i,module+offset]
	 newdata[j*2] = 2*peak_loc - data[i,module+offset]
 	 j=j+1	

   return newdata

###------------------------------------------------------------------

def cutcriteria_result(data,mod,diode,cutcriteria,peak_loc,rms,binsize,patch):

   """
   """

   npts=len(data[:,0])
   offset=3
   print mod, diode, peak_loc, rms
   if cutcriteria == -1.0:
      newcutcriteriah = peak_loc + 5.0*rms
      newcutcriterial = peak_loc - 5.0*rms
   elif cutcriteria == -3.0:
      newcutcriteriah = peak_loc + 3.0*rms
      newcutcriterial = peak_loc - 3.0*rms
   else:
      newcutcriterial = cutcriteria
      newcutcriteriah = 100

   cutf=string.join([path_to_cuts,'weather_sigma_cutlist_',\
                     str(binsize),'_patch_',patch,'.txt'],'')
   acceptf=string.join([path_to_cuts,'weather_sigma_acceptlist_',\
                        str(binsize),'_patch_',patch,'.txt'],'')
   outf=string.join([path_to_cuts,'howmanysigma_',str(binsize),'_patch_',\
                     patch,'.txt'],'')

   for f in [cutf,acceptf,outf]:
      if os.path.exists(f): os.rename(f,'%s.old'%f)

   cutfile = open(cutf,'w')
   acceptfile = open(acceptf,'w')
   outfile = open(outf,'w')

   cutfile.write('#Weather cutting on criteria: accept RMS smaller than %e \n' % (cutcriteria))
   for i in range(npts):
      if data[i,mod+offset] > newcutcriteriah or data[i,mod+offset] < newcutcriterial:
         cutfile.write(' %d %d %d %d %e \n' % (data[i,0],data[i,1],mod,diode,data[i,mod+offset]))
      else:
         acceptfile.write(' %d %d %d %d %e \n' % (data[i,0],data[i,1],mod,diode,data[i,mod+offset]))
      outfile.write(' %d %d %d %d %e \n ' % (data[i,0],data[i,1],mod,diode,(data[i,mod+offset]-peak_loc)*1.0/rms*1.0))
   acceptfile.close()   
   cutfile.close()
   outfile.close()

   #return newcutcriteriah,newcutcriterial
   return newcutcriterial,newcutcriteriah

####-----------------------------------------------------------------

def histogramthisbiatch(data,module,label,max,low,high):

   """
   """

   #params={'axes.labelsize':10,'ytick.labelsize':8,'xtick.labelsize':8,'legend.font size':6}
   #pylab.rcParams.update(params)

   offset=3
   # Plot the histogram
   plt.hist(data[:,module+offset],bins=nbins,histtype='bar',range=(0,max))
   plt.ylabel(label)
   plt.axvline(low,color='r')
   plt.axvline(high, color='r')

   return

####-------------------------------------------------------------

def histogramanybiatch(data,module,label,max,vline):

   """
   """

   params={'axes.labelsize':8,'ytick.labelsize':6,'xtick.labelsize':6,'legend.font size':6}
   pylab.rcParams.update(params)

   plt.hist(data,bins=nbins,histtype='bar',range=(0,max))
   plt.axvline(vline,color='r')

   return

####-------------------------------------------------------------

def plotbydate(data,module,label,max,low,high,fixed=None):

   """
   """

   offset = 3
   counter = 0
   lentd = len(data[:,0])
   leny = len(data[0,:])
   newar=zeros([lentd,leny])
   for i in range(lentd):
      if data[i,module+offset] != -100.0:
         newar[counter,:] = data[i,:]
         counter = counter + 1
   datear=getdate(newar[0:counter,2])
   plt.plot_date(datear,newar[0:counter,module+offset],'b+')
   plt.axhline(low,color='r')
   plt.axhline(high,color='r')
   plt.ylim(low*0.25,high*2.5)
   if low < 0.0: plt.ylim(low*2.5,high*2.5)
      #print 'low negative %20.10f' % low, high>low
   if (low > high):
      print 'low > high (%20.10f > %20.10f) --- exiting' % (low,high)
      sys.exit(1)
   if fixed is not None: plt.ylim(-1e-6,fixed)
   plt.ylabel(label)

   return

####-------------------------------------------------------------

def getdate(array):

   """
   """

   n=len(array)
   datearr=[]
   for i in range(n):
        jd = array[i] + 2400000.5
        unixtime = (jd - 2440587.5)*86400
        datething=datetime.datetime.utcfromtimestamp(unixtime)
        datearr.append(datething)

   return datearr

####-------------------------------------------------------------

def format1(x, pos):

   """
   This allows the plots to use scientific notation
   The two args are the value and tick position
   """

   return '%2.0e' % (x)

####-------------------------------------------------------------

def exponent(x):

   """
   Return the base-10 exponent of a float
   Tested somewhat
   """

   from math import floor,log10

   return int(floor(log10(abs(x))))


###########################################################################
######### MAIN PART IS HERE ###############################################
###########################################################################

if __name__ == '__main__':

   from optparse import OptionParser

   usage = 'Looks for files:Qx_STD_STD_6a_binsz_.txt.txt\n -1 = cut on 1 sigma\n -2 = cut on gaussian FWHM criteria: 3*sigma\n other: low criteria is your value, high value is 100\n creates files: weather_cutlist_basoffilename.txt\n            and weather_acceptlist_baseoffilename.txt\n'

   parser = OptionParser(usage=usage)

   parser.add_option('-p', '--patch', dest='patch',\
                     help='specify QUIET CMB patch', type='string',\
                     metavar='(2a|4a|6a|7b)')

   parser.add_option('-t', '--timescale', dest='timescale',\
                     help='specify time scale for binning',\
                     type='int',metavar='(10|30)')

   parser.add_option('-c', '--cut', dest='cut',\
                     help='specify cut criterion',\
                     type='int',metavar='(-1|-2)')

   parser.add_option('-f', '--fixed', dest='fixed',\
                     help='specify fixed axis range',\
                     type='float',metavar='(20e-7 works well)')

   parser.add_option('-b', '--band', dest='band',\
                     help='specify band', type='string',\
                     metavar='(Q|W)')

   parser.add_option('-o', '--outdir', dest='outdir',\
                     help='specify output directory', type='string',\
                     metavar='.')

   parser.set_defaults(patch=None,timescale=None,cut=None,band=None,\
                       outdir=None,fixed=None)

   # Parse command line
   (options, args) = parser.parse_args()

   patch=options.patch
   timescale=options.timescale
   cut=options.cut
   band=options.band
   outdir=options.outdir
   fixed=options.fixed

   binsize=timescale
   cutcriteria=cut
   patch=patch

   outdir='%s/' % os.path.abspath(outdir)

   path_to_stats=outdir
   path_to_cuts=outdir
   path_to_plots=outdir

   if band == 'Q':
      nmodules=N_MODULE_MAX_Q-2 # i.e. ignore TT
   elif band == 'W':
      nmodules=N_MODULE_MAX_W-6 # i.e. ignore TT
   else:
      nmodules=None

   #binsize=int(sys.argv[1])
   #cutcriteria=int(sys.argv[2])
   #patch=sys.argv[3]

   #Q1,Q2,U1,U2 = file_read_utils.read_STDSTD_statsfile(binsize,patch)
#   Q1,Q2,U1,U2 = file_read_utils.read_STDSTD_nobad_statsfile(binsize,patch)
   Q1,Q2,U1,U2 = file_read_utils.read_STDSTD_noise_statsfile(binsize,patch,\
                                                             outdir)

   # Set up data for the loop over diodes
   diodes=[(1,Q1,'Q1',0),\
           (2,Q2,'Q2',3),
           (3,U1,'U1',1),
           (4,U2,'U2',2)]

   # These are known-to-be-dead diodes/modules
   dead=set([(7,'Q1'),(7,'Q2'),(7,'U1'),(7,'U2'),\
             (8,'Q1'),(8,'Q2'),(8,'U1'),(8,'U2'),\
             (28,'Q1'),(28,'Q2'),(28,'U1'),(28,'U2'),\
             (42,'Q1'),(42,'Q2'),(42,'U1'),(42,'U2'),\
             (81,'Q1'),(81,'Q2'),(81,'U1'),(81,'U2'),\
             (4,'U2'),(38,'Q1'),(40,'Q1')])


   # Set max

   if binsize == 300 or binsize == 480:
      mmax=2E-7
   else:
      mmax=8E-7

   if band=='W':
      #max=80e-7
      mmax=fixed

   # Set up fitting parameters
   init_amp=8E8   # Peak amplitude first guess
   init_sig=2E-7  # Gaussian's sigma first guess
   maxfev=800 # Otherwise use the default fitting niter of 800


   # Now loop over modules
   print 'Analyzing %s-band with %i poln modules' % (band,nmodules)
   print 'entered: ', binsize, cutcriteria, mmax
   
   tol = 0 # Tolerance for missing CESs out of the histogram
   for mod in range(nmodules):
      # Set names of output files
      title=string.join(['Patch ', patch, ', Module ',str(mod),' timescale:',str(binsize),'s'],'')
      figname=string.join([path_to_plots,'Patch_',patch,'_Module',str(mod),'_',str(binsize),'_hist.png'],'')
      title_p=title
      figname_p=string.join([path_to_plots,'Patch_',patch,'_Module',str(mod),'_',str(binsize),'_date.png'],'')

      ### PLOT HISTOGRAMS

      # Set up plotting parameters
      colour='white'
      params={'axes.labelsize':10,'ytick.labelsize':8,\
              'xtick.labelsize':8,'legend.font size':6}
      pylab.rcParams.update(params)

      for (panel,diode,dname,posn) in diodes:
         # Count the number of CESs
         (nces,ndiodeplustwo)=numpy.shape(diode)
         #ndiode = ndiodeplustwo - 2

         # Set up the subpanels
         if (mod,dname) in dead: # Mask dead modules
            colour='red'
         else:
            colour='white'

         plt.subplot(4,1,panel)#,axisbg=colour) -- broken by next line
         # Set scientific notation for the xtics
         plt.subplot(4,1,panel).xaxis.set_major_formatter(\
           formatter=ticker.FuncFormatter(format1))
         plt.subplot(4,1,panel).set_axis_bgcolor(colour)

         # Calculate the histogram
         #     and find the its peak
         # Points are necessarily > 0
         # Only search up to mmax
         # Data are inside the array called diode
         peak_loc,counter,hist,maxval=find_hist_max(diode,mod,0,mmax)

         print '%s %i' % (dname,counter)

         # Reshape the histogram
         newhist=get_histogram(hist)

         #if (mod,dname) == (0,'U1'):
         #   print newhist,'n',numpy.shape(newhist)
         #   print diode,numpy.shape(diode)
         #   print hist
         #   print peak_loc,counter,maxval
         #   #sys.exit(1)

         #print newhist

         # Keep track of the number of missing CESs
         if int(numpy.sum(newhist)) != nces:
            otol=tol
            tol = max(nces-int(numpy.sum(newhist)),tol)
            print '*** %i %s missing (%i ->) %i CESs' % (mod,dname,otol,tol)

         #assert(nces-1<=int(numpy.sum(newhist))<=nces),\
         #   'Some of our CESs are missing (%i < %i) for %i %s!' \
         #   % (int(numpy.sum(newhist)),nces,mod,dname)

         # Fit a Gaussian to the histogrammed data
         #     and return the 3 parameters of the best-fit gaussian
         results=gaussfit.get_fitparams(newhist,peak_loc,init_amp,init_sig,\
                                        maxfev=maxfev)

         # Evaluate the gaussian at each histogram bin location
         func = gaussfit.gaussvalue(results,hist[1][:])

         # Calculate the positions of the low and high bars
         #     and write files
         low,high=cutcriteria_result(diode,mod,posn,cutcriteria,peak_loc,\
                                     results[1],binsize,patch)
         #print 'th',peak_loc,counter,maxval,init_amp,init_sig,mmax,low,high

         # Plot the histogram and bars, but only up to mmax
         histogramthisbiatch(diode,mod,dname,mmax,low,high)

         # Set the xrange
         if fixed is not None: plt.xlim(-1e-6,fixed)

         # Plot the curve
         plt.plot(hist[1][:],func,'m-')
         if panel == 1: plt.title(title)

      # Output plot and tidy up
      plt.savefig(figname)
      plt.clf()
      plt.close('all')

      ### PLOT BY TIME
      for (panel,diode,dname,posn) in diodes:
         if (mod,dname) in dead: # Mask dead modules
            colour='red'
         else:
            colour='white'
         # Set up and format the subplot
         plt.subplot(4,1,panel)
         plt.subplot(4,1,panel).yaxis.set_major_formatter(\
           formatter=ticker.FuncFormatter(format1))
         plt.subplot(4,1,panel).set_axis_bgcolor(colour)

         peak_loc,counter,hist,maxval=find_hist_max(diode,mod,0,mmax)
         newhist=get_histogram(hist)
         results=gaussfit.get_fitparams(newhist,peak_loc,init_amp,init_sig,\
                                        maxfev=maxfev)
         #func1 = gaussfit.gaussvalue(results,hist[1][:])
         low,high=cutcriteria_result(diode,mod,posn,cutcriteria,peak_loc,\
                                     results[1],binsize,patch)
         #print 'tt',peak_loc,counter,maxval,init_amp,init_sig,mmax,low,high
         plotbydate(diode,mod,dname,mmax,low,high,fixed=fixed)
         if panel == 1: plt.title(title_p)

      # Output plot and tidy up
      plt.savefig(figname_p)
      plt.clf()
      plt.close('all')

   print '\nMaximum loss of CESs was %i\n' % tol

   print 'Finished!\n'
   sys.exit(0)

############################################################################
############ END ###########################################################
############################################################################
