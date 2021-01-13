import sys, quiet, math
import quietPyExt
import numpy
from numpy import zeros
from numpy import fft,average,log10,sqrt
import fitPowerSpectrum
eulerm = 0.57721566490153286060

# get nearest multiple of 50 -------------------------------------
def get_nffts(nsamp):
    cut=int(nsamp/50)
    intlength = (cut-1)*50
    return intlength

# create a frequency array ----------------------------------------
def get_freqarr(l):
   from numpy import zeros 
   freqarr=zeros([l])
   for i in range(0,l):
      freqarr[i] = (i*1.0/l*1.0)*25.0
   return freqarr

# estimate #1 of the noise floor  ---------------------------------
def get_floor(array,lowfreq,highfreq):
   x=len(array[:,0])
   y=len(array[0,:])

   lowit=int(lowfreq*y/25)
   higit=int(highfreq*x/25)
   noisearr=zeros([x])

   for i in range(0,x):
      noisearr[i] = average(array[i,lowit:higit])

   return noisearr
       
# integrate power in the bump ---------------------------------------
def get_bump_power(array, noisearr, scanfreq, lowfreq,highfreq):
   x=len(array[:,0])
   y=len(array[0,:])

   scanit=int(scanfreq*y/25)
   lowit=int(lowfreq*y/25)
   highit=int(highfreq*y/25)

   power=zeros([x])
   for i in range(0,x):
      for j in range(lowit, highit):
         power[i] = power[i] + (array[i,j]-noisearr[i])
      power[i] = power[i]/(highit-lowit)

   return power

# create array with 1/f spectrum --------
def get_spectrum(spectrum,freqarr,l):
    func = zeros([l])
    for i in range(0,l):
       #func[i] = float(i*1.0)**(-1*spectrum[2]) + spectrum[1]
	func[i] = (spectrum[1])*(1+(spectrum[0]/freqarr[i])**spectrum[2])
    return func

# calculate the chi2 to 1/f spectrum to data, log -------------------
def get_chi_sq(data,spectrum,freqarr,scanfreq,toggle=0):
    n=len(data)
    pearson = zeros([2,n])
    chiarray=zeros([2,9])
    #second=zeros([n])
    #second[:]=abs(spectrum[1] + spectrum[1]*(spectrum[0]/freqarr[:])**spectrum[2])
    #pearson[0,:] = (data[:] - log10(abs(second[:])) - eulerm)**2
    #pearson[1,:] = ((data[:] - log10(second[:]) - eulerm)**2)/abs(log10(second[:]))
    pearson[0,:] = (fitPowerSpectrum._residuals(spectrum,freqarr,data))**2
    pearson[1,:] = ((fitPowerSpectrum._residuals(spectrum,freqarr,data))**2)

    int_subsets=zeros([9])
    freq_subsets=[0.02,(scanfreq - 0.5*scanfreq),scanfreq,(scanfreq+0.5*scanfreq),\
	(2.0*scanfreq-0.2*scanfreq),2.0*scanfreq, (2.0*scanfreq + 0.2*scanfreq),7.0, 24.0]
    for i in range(0,9):
       int_subsets[i] = int(freq_subsets[i] * (n/25))
    chiarray[0,0]=average(pearson[0,1:n])
    chiarray[1,0]=average(pearson[1,1:n])
    if toggle==0:
       for i in range(1,9):
	 chiarray[0,i] = average(pearson[0,int_subsets[i-1]:int_subsets[i]])
	 chiarray[1,i] = average(pearson[1,int_subsets[i-1]:int_subsets[i]])
    else:
      chiarray[:,1:8] = 0.0
    return chiarray[0,:],chiarray[1,:]
          
# plot the fft and overlay the 1/f spectrum ---------------------
def plot_ffts(data, spectrum, freqarr, run_id, scan_id,scanfreq,j):
      from matplotlib import pyplot as plt
      import string
      n=len(freqarr)
      start=1
      chiarray=zeros([2,4,9])
      plt.subplot(2,2,1)
      holder=get_spectrum(spectrum[0,:],freqarr,n)
      chiarray[0,0,:],chiarray[1,0,:] = get_chi_sq(log10(data[0,:]),spectrum[0,:],freqarr,scanfreq,1)
      plt.loglog(freqarr[start:n],data[0,start:n],'r+')
      plt.loglog(freqarr[start:n],holder[start:n],'k')
      plt.axvline(x=0.02,color='b')
      plt.axvline(x=(scanfreq-0.5*scanfreq),color='g')
      plt.axvline(x=scanfreq,color='r')
      plt.axvline(x=scanfreq+0.5*scanfreq,color='c')
      plt.axvline(x=2.0*scanfreq-0.2*scanfreq,color='m')
      plt.axvline(x=2.0*scanfreq,color='y')
      plt.axvline(x=2.0*scanfreq+0.2*scanfreq,color='b')
      plt.axvline(x=7.0,color='g')
      plt.ylabel(string.join(['DD1, CHI2:',str(chiarray[0,0,0])],''))
  
      plt.subplot(2,2,2)
      holder=get_spectrum(spectrum[1,:],freqarr,n)
      plt.loglog(freqarr[start:n],data[1,start:n],'r+')
      plt.loglog(freqarr[start:n],holder[start:n],'k')
      chiarray[0,1,:],chiarray[1,1,:] = get_chi_sq(log10(data[1,:]),spectrum[1,:],freqarr,scanfreq,1)
      plt.ylabel(string.join(['DD2, CHI2:',str(chiarray[0,1,0])],''))

      plt.subplot(2,2,3)
      holder=get_spectrum(spectrum[2,:],freqarr,n)
      plt.loglog(freqarr[start:n],data[2,start:n],'r+')
      plt.loglog(freqarr[start:n],holder[start:n],'k')
      chiarray[0,2,:],chiarray[1,2,:] = get_chi_sq(log10(data[2,:]),spectrum[2,:],freqarr,scanfreq,1)
      plt.ylabel(string.join(['DD3, CHI2:',str(chiarray[0,2,0])],''))
      
      plt.subplot(2,2,4)
      holder=get_spectrum(spectrum[3,:],freqarr,n)
      plt.loglog(freqarr[start:n],data[3,start:n],'r+')
      plt.loglog(freqarr[start:n],holder[start:n],'k')
      chiarray[0,3,:],chiarray[1,3,:] = get_chi_sq(log10(data[3,:]),spectrum[3,:],freqarr,scanfreq,1)
      plt.ylabel(string.join(['DD4, CHI2:',str(chiarray[0,3,0])],''))
      plt.savefig(string.join(['FFT_',str(run_id),'_',str(scan_id),'_module',str(j),'.png'],''))
      plt.clf()
      plt.close('all')

###########################################################################
######### MAIN PART IS HERE ###############################################
###########################################################################

if __name__ =="__main__":

  # main part of the code starts
  level2 = "/data4/quiet/level2/"
  #runlist = quiet.L2Runlist("/data4/quiet/runlist_l2.txt")
  #runlist = quiet.L2Runlist("/data4/quiet/runlist_l2_tworuns.txt")
  runlist = quiet.L2Runlist(sys.argv[1])
  segments = runlist.segments()
  Statsfile=open(sys.argv[2],"w+")
  Statsfile.write("run seg mod DD  Fknee  sigma     alpha     full    lowest lower    SS    middle    harmon   high    highest ...same normalized ....\n")
  
  #loop through the segment IDs
  for (run_id,segment_id) in segments:

	segment_info = segments[(run_id,segment_id)]
	filename = segment_info.filename(level2)
	print filename
        #SSfile=open('BumpContrib.txt',"w+")
        
	data = quiet.l2_read(filename)

        # set up the FFT
        nffts = get_nffts(data[0].nsamp)
        nfreqs = nffts/2
        freqarr=get_freqarr(nfreqs)
	# get scan frequency, print it out
	scanfreq = quietPyExt.get_scan_freq(data[0])

	for j in range(0,19):
	   fftarr=zeros([4,nffts])
	   spectrum=zeros([4,3])
	   spectrum[:,:] = 1.0
           chi=zeros([2,4,9])
           # calculate the ffts
           fftarr[0,:]=(abs(numpy.fft.fft(data[j].tod.Q1[0:nffts])))**2
           fftarr[1,:]=(abs(numpy.fft.fft(data[j].tod.U1[0:nffts])))**2
           fftarr[2,:]=(abs(numpy.fft.fft(data[j].tod.U2[0:nffts])))**2
           fftarr[3,:]=(abs(numpy.fft.fft(data[j].tod.Q2[0:nffts])))**2
   	   # get white noise floor, knee, slope
	   try:
              spectrum[0,:] = fitPowerSpectrum.fitPowerSpectrum(data[j].tod.Q1[0:nffts],50)
	      chi[0,0,:],chi[1,0,:] = get_chi_sq(log10(fftarr[0,0:nfreqs]),spectrum[0,:],freqarr,scanfreq,0)
           except:
              pass
           try:
              spectrum[1,:] = fitPowerSpectrum.fitPowerSpectrum(data[j].tod.U1[0:nffts],50)
              chi[0,1,:],chi[1,1,:] = get_chi_sq(log10(fftarr[1,0:nfreqs]),spectrum[1,:],freqarr,scanfreq,0)
           except:
              pass
           try:
              spectrum[2,:] = fitPowerSpectrum.fitPowerSpectrum(data[j].tod.U2[0:nffts],50)
	      chi[0,2,:],chi[1,2,:] = get_chi_sq(log10(fftarr[2,0:nfreqs]),spectrum[2,:],freqarr,scanfreq,0)
           except:
              pass
           try:
              spectrum[3,:] = fitPowerSpectrum.fitPowerSpectrum(data[j].tod.Q2[0:nffts],50)
	      chi[0,3,:],chi[1,3,:] = get_chi_sq(log10(fftarr[3,0:nfreqs]),spectrum[3,:],freqarr,scanfreq,0)
           except:
              pass
	   # plot the FFT, overlay 1/f spectrum
           toggle=int(sys.argv[3])
           if toggle==1:
              try:
                 plot_ffts(fftarr[:,0:nfreqs],spectrum,freqarr, run_id,segment_id,scanfreq,j)
	      except:
 	         pass
	   # get various Chi2's
 	   for k in range(0,4):
	      Statsfile.write('%d %d %d %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n' \
		% (run_id, segment_id, j, k, \
		spectrum[k,0], spectrum[k,1],spectrum[k,2], \
		chi[0,k,0],chi[0,k,1],chi[0,k,2],chi[0,k,3],chi[0,k,4],\
		chi[0,k,5],chi[0,k,6],chi[0,k,7], chi[0,k,8], \
		chi[1,k,0],chi[1,k,1],chi[1,k,2],chi[1,k,3],chi[1,k,4],\
		chi[1,k,5],chi[1,k,6],chi[1,k,7],chi[1,k,8]))
 	   # estimate power in bump
           #lowfreq = scanfreq - (scanfreq/10.0)
           #highfreq = scanfreq + (scanfreq/10.0)
	   #bump = get_bump_power(fftarr,spectrum,scanfreq,lowfreq,highfreq)
 	   #SSfile.write('Scan: %d Seg: %d Module: %d %f %f %f %f \n' % (run_id, segment_id,j,bump[0],bump[1],bump[2],bump[3]))
  	

  Statsfile.close()

############################################################################
############ END ###########################################################
############################################################################

