import sys,numpy,pylab,string
from matplotlib import pyplot as plt
from numpy import zeros
import binning
import quiet

# --------------------------------------------------
def get_binned_tod(data,binsize,j):
#   print 'in getbinnedtod'
   n=len(data[j].tod[0])
   nbins=n/(binsize+1)
   binneddata = zeros([4,nbins])
   for i in range(0,4):
      binneddata[i]=binning.bin_given_binsize(data[j].tod[i],binsize)

   return binneddata

# --------------------------------------------------
def get_binned_TP(data,binsize,j):
#   print 'in getbinnedtod'
   n=len(data[j].TP[0])
   nbins=n/(binsize+1)
   binneddata = zeros([4,nbins])
   for i in range(0,4):
      binneddata[i]=binning.bin_given_binsize(data[j].TP[i],binsize)

   return binneddata

# --------------------------------------------------
def get_binned_tod_std(data,binsize,j):
#   print 'in getbinnedtod'
   n=len(data[j].tod[0])
   nbins=n/(binsize+1)
   binneddata = zeros([4,nbins])
   for i in range(0,4):
      binneddata[i]=binning.bin_given_binsize_std(data[j].tod[i],binsize)

   return binneddata


# --------------------------------------------------
def plot_tods(todarr,run_id,segment_id,mod,label): 
 #     print 'in plotting'
      params={'axes.labelsize':12,'axes.titlesize':12,'ytick.labelsize':8,'xtick.labelsize':8,'legend.fontsize':10,'markerscale':0.5,'figure.subplot.wspace':0.4,'figure.subplot.hspace':0.3}
      plt.rcParams.update(params)

      try:
        plt.subplot(2,2,1)
        plt.plot(todarr[0,:])
        plt.ylabel('DD1 (V)')
        plt.title(string.join([label, ' for Run ',str(int(run_id)),' , Seg ',str(int(segment_id)), ', RQ: ',str(int(mod))],''))
      except:
	print 'nope'
        pass
      try:
        plt.subplot(2,2,2)
  	plt.plot(todarr[1,:])
        plt.ylabel('DD2 (V)')	
      except:
	print 'nope'
        pass
      try:
        plt.subplot(2,2,3)
	plt.plot(todarr[2,:])
        plt.ylabel('DD3 (V)')
      except:
	print 'nope'
        pass
      try:
        plt.subplot(2,2,4)
	plt.plot(todarr[3,:])
        plt.ylabel('DD4 (V)')
      except:
	print 'nope'
        pass
      plt.savefig(string.join([label,'_TOD_',str(run_id),'_',str(segment_id),'_module',str(mod),'.png'],''))
      plt.clf()
      plt.close('all')


###########################################################################
######### MAIN PART IS HERE ###############################################
###########################################################################

if __name__ =="__main__":

  # main part of the code starts
        level2 = "/data4/quiet/level2/"
        runlist = quiet.L2Runlist(sys.argv[1])
        binsize = int(sys.argv[2])
        segments = runlist.segments()
  #loop through the segment IDs
#  for (run_id,segment_id) in segments:
#  thing=zeros([4,2])
#  thing[0,:] = [691,2]
#  thing[1,:] = [764,1]
#  thing[2,:] = [766,1]
#  thing[3,:] = [1486,3]
#  for i in range(4): 
        run_id = 508
        segment_id = 1   
      	segment_info = segments[(run_id,segment_id)]
        filename = segment_info.filename(level2)
        print filename

        data = quiet.l2_read(filename)#,TP=True)

	# plot the TOD
        for j in range(0,1):
	   
           print 'plotting module: ', j
	   todarr = get_binned_tod(data,binsize,j)
	   #todarr = get_binned_tod_std(data,binsize,j)
	   #todarr = get_binned_TP(data,binsize,j)
           try:
           	#plot_tods(todarr,run_id,segment_id,j,'STD')
   	        plot_tods(todarr,run_id,segment_id,j,'TOD')

           except:
		print 'couldnt plot'
              	pass

############################################################################
############ END ###########################################################
############################################################################
  
