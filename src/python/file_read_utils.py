#!/usr/bin/python

import os,sys,pylab,numpy,matplotlib,string,datetime,numpy,gaussfit
from matplotlib import pyplot as plt
from numpy import zeros
from numpy import *

user=os.getenv('USER')
if user == 'newburgh' : user='%s/repository'%user

sys.path.append('/usit/titan/u1/%s/quiet_svn/oslo/src/python/laura_programs'%user)
import plot_cut_statistics_various


path_to_stats = '/usit/titan/u1/%s/quiet_svn/oslo/src/python/weather_cutting/statisticsfiles/' % user
path_to_cuts = '/usit/titan/u1/%s/quiet_svn/oslo/src/python/weather_cutting/cuts/' % user
path_to_plots = '/usit/titan/u1/%s/quiet_svn/oslo/src/python/weather_cutting/plots/' % user
nbins=200
path_to_noise='/usit/titan/u1/%s/quiet_svn/oslo/src/python/weather_cutting/noisetest/' % user

typeb_reject='/usit/titan/u1/newburgh/repository/quiet_svn/oslo/src/python/weather_cutting/auxil/typeb_badness_forcutting.txt'

####------------------------------------------------------------------
def getbadtypeb(typeb_reject):
    A=pylab.load(typeb_reject)
    n=len(A[:,0])
    ntypeb = 1800
    plotarr=zeros([19,4,ntypeb,2])
    plotarr[:,:,:,:] = -100
    counter=zeros([19,4])
    countert = 0
    maxc = 0
    for i in range(0,n):
       for mod in range(0,19):
          for diode in range(0,4):
             if int(A[i,2]) == mod and int(A[i,3]) == diode:
                plotarr[mod,diode,counter[mod,diode],:] = [A[i,0], A[i,1]]
                counter[mod,diode] = counter[mod,diode] + 1
                if counter[mod,diode] > maxc:
                    maxc = counter[mod,diode]

    return plotarr, counter, maxc

####-------------------------------------------------------------------
def makebad():
    listobad=[329,330,331,332,333,334,335,336,337,338,339,340,341,563,\
        564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,\
        580,581,582,583,584,585,586,587,588,589,590,591,592,593]
    setobad=set(listobad)
    #setobad.add(639)
    return setobad

def read_STDSTD_nobad_statsfile(binsize,patch):

   print 'in nobad'
   filename=string.join([path_to_stats,'Q1_STD_STD_',patch,'_',str(binsize),'_allmods.txt'],'')
   filename2=string.join([path_to_stats,'Q2_STD_STD_',patch,'_',str(binsize),'_allmods.txt'],'')
   filename3=string.join([path_to_stats,'U1_STD_STD_',patch,'_',str(binsize),'_allmods.txt'],'')
   filename4=string.join([path_to_stats,'U2_STD_STD_',patch,'_',str(binsize),'_allmods.txt'],'')
   A=matplotlib.mlab.load(filename)
   B=matplotlib.mlab.load(filename2)
   C=matplotlib.mlab.load(filename3)
   D=matplotlib.mlab.load(filename4)
   plotarr,counter,maxc =getbadtypeb(typeb_reject)
   badrq11 = matplotlib.mlab.load('/usit/titan/u1/newburgh/repository/quiet_svn/oslo/src/python/pswitchoff/rq11_ps21_badscans.txt')
   badrq12 = matplotlib.mlab.load('/usit/titan/u1/newburgh/repository/quiet_svn/oslo/src/python/pswitchoff/rq12_ps21_badscans.txt')
   # set values of all to -100 if scan is in the list of bad scans
   badlist=makebad()
   lena = len(A[:,0])
   for i in range(0,lena):
      if A[i,0] not in badlist:
         pass
      else:
         A[i,2:-1] = -100
         A[i,-1] = -100
	 B[i,2:-1] = -100
         B[i,-1] = -100
	 C[i,2:-1] = -100
         C[i,-1] = -100
         D[i,2:-1] = -100
         D[i,-1] = -100

   # set RQ12, RQ11 values to -100 if the phase switches were off, or if its in the typeb list
   offset=3
   for data in range(0,len(A[:,0])):
      for module in range(0,2):
         diode=0
         for typeb in range(0,int(counter[module,diode])):
            if plotarr[module,diode,typeb,0] == A[data,0] and plotarr[module,diode,typeb,1] == A[data,1]:
                   A[data,module+offset] = -100
         diode=3
         for typeb in range(0,int(counter[module,diode])):
            if plotarr[module,diode,typeb,0] == B[data,0] and plotarr[module,diode,typeb,1] == B[data,1]:
                   B[data,module+offset] = -100
         diode=1
         for typeb in range(0,int(counter[module,diode])):
            if plotarr[module,diode,typeb,0] == C[data,0] and plotarr[module,diode,typeb,1] == C[data,1]:
                   C[data,module+offset] = -100
         diode=2
         for typeb in range(0,int(counter[module,diode])):
            if plotarr[module,diode,typeb,0] == D[data,0] and plotarr[module,diode,typeb,1] == D[data,1]:
                   D[data,module+offset] = -100
      for rq11 in range(0,len(badrq11[:,0])):
         if  A[data, 0] == badrq11[rq11,0] and A[data,1] == badrq11[rq11,1]:
            print A[data,0], A[data,1]
            A[data,11+offset] = -100
            B[data,11+offset] = -100
            C[data,11+offset] = -100
            D[data,11+offset] = -100
      for rq12 in range(0,len(badrq12[:,0])):
         if  A[data, 0] == badrq12[rq12,0] and A[data, 1] == badrq12[rq12,1]:
            A[data,12+offset] = -100
            B[data,12+offset] = -100
            C[data,12+offset] = -100
            D[data,12+offset] = -100

   return A, B, C, D

####-------------------------------------------------------------------
def read_STDSTD_noise_statsfile(binsize,patch):

   filename=string.join([path_to_noise,'Q1_STD_STD_',patch,'_',str(binsize),'.txt'],'')
   filename2=string.join([path_to_noise,'Q2_STD_STD_',patch,'_',str(binsize),'.txt'],'')
   filename3=string.join([path_to_noise,'U1_STD_STD_',patch,'_',str(binsize),'.txt'],'')
   filename4=string.join([path_to_noise,'U2_STD_STD_',patch,'_',str(binsize),'.txt'],'')
   A=matplotlib.mlab.load(filename)
   B=matplotlib.mlab.load(filename2)
   C=matplotlib.mlab.load(filename3)
   D=matplotlib.mlab.load(filename4)

   return A, B, C, D



####-------------------------------------------------------------------
def read_STDSTD_statsfile(binsize,patch):

   filename=string.join([path_to_stats,'Q1_STD_STD_',patch,'_',str(binsize),'_allmods.txt'],'')
   filename2=string.join([path_to_stats,'Q2_STD_STD_',patch,'_',str(binsize),'_allmods.txt'],'')
   filename3=string.join([path_to_stats,'U1_STD_STD_',patch,'_',str(binsize),'_allmods.txt'],'')
   filename4=string.join([path_to_stats,'U2_STD_STD_',patch,'_',str(binsize),'_allmods.txt'],'')
   A=matplotlib.mlab.load(filename)
   B=matplotlib.mlab.load(filename2)
   C=matplotlib.mlab.load(filename3)
   D=matplotlib.mlab.load(filename4)

   return A, B, C, D
####-------------------------------------------------------------------

def read_MeanSTD_statsfile(binsize,patch):

   filename=string.join([path_to_stats,'Q1_Mean_STD_',patch,'_',str(binsize),'_allmods.txt'],'')
   filename2=string.join([path_to_stats,'Q2_Mean_STD_',patch,'_',str(binsize),'_allmods.txt'],'')
   filename3=string.join([path_to_stats,'U1_Mean_STD_',patch,'_',str(binsize),'_allmods.txt'],'')
   filename4=string.join([path_to_stats,'U2_Mean_STD_',patch,'_',str(binsize),'_allmods.txt'],'')
   A=matplotlib.mlab.load(filename)
   B=matplotlib.mlab.load(filename2)
   C=matplotlib.mlab.load(filename3)
   D=matplotlib.mlab.load(filename4)

   return A, B, C, D

