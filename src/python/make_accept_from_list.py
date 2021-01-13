#!/usr/bin/python

import sys
import pylab
import numpy
import matplotlib
import string
from numpy import zeros
from numpy import *
#import quartile_finder
#import Laurausefulfunctions



###########################################################################
######### MAIN PART IS HERE ###############################################
###########################################################################

if __name__ =="__main__":

   filename=sys.argv[1]
   A=matplotlib.mlab.load(filename)

   n=len(A[:,1])
   counter=0
   writefile=open('holder.txt','a+')   
   writefile.write(' %d %d %d %d ' % (A[0,0],A[0,1],A[0,2],A[0,3]))
   for i in range(1,n):
      if A[i,0] == A[i-1,0] and A[i,1] == A[i-1,1]: # equal to one before, add line
         writefile.write(' %d %d ' % (A[i,2],A[i,3]))  
         counter= counter+1
      else: # new run, write to file accordingly 
         writefile.write('%d \n' % (counter))
         counter=1
         writefile.write('%d %d %d %d ' % (A[i,0],A[i,1],A[i,2],A[i,3]))
   writefile.close()
############################################################################
############ END ###########################################################
############################################################################


