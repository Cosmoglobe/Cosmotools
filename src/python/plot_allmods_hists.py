#!/usr/bin/python

import os,sys,pylab,numpy,matplotlib,string
import datetime,numpy,gaussfit,file_read_utils
from quiet import N_MODULE_MAX_Q,N_MODULE_MAX_W
from matplotlib import pyplot as plt
from numpy import zeros
from numpy import *

user=os.getenv('USER')
if user == 'newburgh': user='%s/repository'%user

#path_to_stats = '/usit/titan/u1/%s/quiet_svn/oslo/src/python/weather_cutting/statisticsfiles/' % user
#path_to_cuts = '/usit/titan/u1/%s/quiet_svn/oslo/src/python/weather_cutting/cuts/' % user
#path_to_plots = '/usit/titan/u1/%s/quiet_svn/oslo/src/python/weather_cutting/plots/' % user
nbins=200

###-----------------------------------------------------------------

def get_avg_in_range(data,module,min,max):
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
   offset=3
   print 'm%i'%module
   histthing=numpy.histogram(data[:,module+offset],nbins,(min,max)) 
   n=len(histthing[0])
   maxval=numpy.max(histthing[0][2:n])
   binthing=0
   counter=0 
   for i in range(2,n):
      if histthing[0][i] == maxval:
         binthing=histthing[1][i]
   for i in range(0,n):
      if histthing[1][i] <= binthing:
	 counter=counter+1
      histthing[1][i] = histthing[1][i] + 0.5*(histthing[1][1]-histthing[1][0])

   binthing=binthing + 0.5*(histthing[1][1]-histthing[1][0])
   return binthing,counter,histthing,maxval

####-------------------------------------------------------------

def get_histogram(hist,counter,binthing):
   n=len(hist[0][:])
   newhist=zeros([2,n])
   for i in range(n):
      newhist[0][i] = hist[0][i]
      newhist[1][i] = hist[1][i]
   return newhist

####-------------------------------------------------------------

def get_reflected_histogram(hist,counter,binthing):
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

def get_reflected_data(peak,count,data,module):
   offset=3
   j=0
   n=len(data[:,module+offset])
   newdata=zeros([n])
   for i in range(n):
      if data[i,module+offset] <= peak:
         newdata[j] = data[i,module+offset]
	 newdata[j*2] = 2*peak - data[i,module+offset]
 	 j=j+1	
   return newdata
###------------------------------------------------------------------

def cutcriteria_result(data,mod,diode,cutcriteria,peak,rms,binsize,patch):
   npts=len(data[:,0])
   offset=3
   print mod, diode, peak, rms
   if cutcriteria == -1.0:
      newcutcriteriah = peak + 5.0*rms
      newcutcriterial = peak - 5.0*rms
   elif cutcriteria == -3.0:
      newcutcriteriah = peak + 3.0*rms
      newcutcriterial = peak - 3.0*rms
   else:
      newcutcriterial = cutcriteria
      newcutcriteriah = 100
   
   cutfile = open(string.join([path_to_cuts,'weather_sigma_cutlist_',str(binsize),'_patch_',patch,'.txt'],''),"a+")
   acceptfile = open(string.join([path_to_cuts,'weather_sigma_acceptlist_',str(binsize),'_patch_',patch,'.txt'],''),"a+")
   outfile = open(string.join([path_to_cuts,'howmanysigma_',str(binsize),'_patch_',patch,'.txt'],''),"a+")

   cutfile.write('#Weather cutting on criteria: accept RMS smaller than %e \n' % (cutcriteria))
   for i in range(npts):
      if data[i,mod+offset] > newcutcriteriah or data[i,mod+offset] < newcutcriterial:
         cutfile.write(' %d %d %d %d %e \n' % (data[i,0],data[i,1],mod,diode,data[i,mod+offset]))
      else:
         acceptfile.write(' %d %d %d %d %e \n' % (data[i,0],data[i,1],mod,diode,data[i,mod+offset]))
      outfile.write(' %d %d %d %d %e \n ' % (data[i,0],data[i,1],mod,diode,(data[i,mod+offset]-peak)*1.0/rms*1.0))
   acceptfile.close()   
   cutfile.close()
   outfile.close()
   return newcutcriteriah,newcutcriterial

####-----------------------------------------------------------------

def histogramthisbiatch(data,module,label,max,low,high):

   params={'axes.labelsize':10,'ytick.labelsize':8,'xtick.labelsize':8,'legend.font size':6}
   pylab.rcParams.update(params)

   offset=3
   plt.hist(data[:,module+offset],bins=nbins,histtype='bar',range=(0,max))
   plt.ylabel(label)
   plt.axvline(low,color='r')
   plt.axvline(high, color='r')

   return

####-------------------------------------------------------------

def histogramanybiatch(data,module,label,max,vline):

   params={'axes.labelsize':8,'ytick.labelsize':6,'xtick.labelsize':6,'legend.font size':6}
   pylab.rcParams.update(params)

   plt.hist(data,bins=nbins,histtype='bar',range=(0,max))
   plt.axvline(vline,color='r')

####-------------------------------------------------------------

def plotbydate(data,module,label,max,low,high):
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
   plt.ylabel(label)

####-------------------------------------------------------------

def getdate(array):
   n=len(array)
   datearr=[]
   for i in range(n):
        jd = array[i] + 2400000.5
        unixtime = (jd - 2440587.5)*86400
        datething=datetime.datetime.utcfromtimestamp(unixtime)
        datearr.append(datething)
   return datearr


###########################################################################
######### MAIN PART IS HERE ###############################################
###########################################################################

if __name__ == '__main__':

   from optparse import OptionParser

   usage = 'Looks for files:Qx_Mean_STD_6a_binsz_.txt and Qx_STD_STD_6a_binsz_.txt.txt\n -1 = cut on 1 sigma\n -2 = cut on gaussian FWHM criteria: 3*sigma\n other: low criteria is your value, high value is 100\n creates files: weather_cutlist_basoffilename.txt\n            and weather_acceptlist_baseoffilename.txt\n'

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

   parser.add_option('-b', '--band', dest='band',\
                     help='specify band', type='string',\
                     metavar='(Q|W)')

   parser.add_option('-o', '--outdir', dest='outdir',\
                     help='specify output directory', type='string',\
                     metavar='.')

   parser.set_defaults(patch=None,timescale=None,cut=None,band=None,\
                       outdir=None)

   # Parse command line
   (options, args) = parser.parse_args()

   patch=options.patch
   timescale=options.timescale
   cut=options.cut
   band=options.band
   outdir=options.outdir

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

   print 'Analyzing %s-band with %i poln modules' % (band,nmodules)

   if binsize == 300 or binsize == 480:
      max=2E-7
   else:
      max=8E-7

   print 'entered: ', binsize, cutcriteria, max
   #Q1,Q2,U1,U2 = file_read_utils.read_STDSTD_statsfile(binsize,patch)
#   Q1,Q2,U1,U2 = file_read_utils.read_STDSTD_nobad_statsfile(binsize,patch)
   Q1,Q2,U1,U2 = file_read_utils.read_STDSTD_noise_statsfile(binsize,patch,\
                                                             outdir)

   for i in range(nmodules):
      iter=str(i)
      title=string.join(['Patch ', patch, ', Module ',iter,' binsize:',str(binsize)],'')
      figname=string.join([path_to_plots,'Patch_',patch,'_Module',iter,'_',str(binsize),'_hist.png'],'')
      title_p=title
      figname_p=string.join([path_to_plots,'Patch_',patch,'_Module',iter,'_',str(binsize),'_date.png'],'')


      ### PLOT HISTOGRAMS

      # Set histogram ranges
      rangemax=8E8
      range2max=2E-7

      # Set up data for the loop over diodes
      diodes=[(1,Q1,'Q1',0),\
              (2,Q2,'Q2',3),
              (3,U1,'U1',1),
              (4,U2,'U2',2)]

      # Set up plotting
      for (panel,diode,dname,posn) in diodes:
         plt.subplot(4,1,panel)
         peak,counter,hist,maxval=find_hist_max(diode,i,0,max)
         newdata=get_histogram(hist,counter,peak)
         results=gaussfit.get_fitparams(newdata,peak,rangemax,range2max)
         func = gaussfit.gaussvalue(results,hist[1][:])
         low,high=cutcriteria_result(diode,i,posn,cutcriteria,peak,results[1],binsize,patch)
         histogramthisbiatch(diode,i,dname,max,low,high)
         plt.plot(hist[1][:],func,'m-')
         if panel == 1: plt.title(title) 
      # Output plot and tidy up
      plt.savefig(figname)
      plt.clf()
      plt.close('all')

      ### PLOT BY TIME
      for (panel,diode,dname,posn) in diodes:
         plt.subplot(4,1,panel)
         peak,counter,hist,maxval=find_hist_max(diode,i,0,max)
         newdata=get_histogram(hist,counter,peak)
         results=gaussfit.get_fitparams(newdata,peak,rangemax,range2max)
         func = gaussfit.gaussvalue(results,hist[1][:])
         low,high=cutcriteria_result(diode,i,posn,cutcriteria,peak,results[1],binsize,patch)
         plotbydate(diode,i,dname,max,low,high)
         if panel == 1: plt.title(title_p)
      # Output plot and tidy up
      plt.savefig(figname_p)
      plt.clf()
      plt.close('all')

############################################################################
############ END ###########################################################
############################################################################
