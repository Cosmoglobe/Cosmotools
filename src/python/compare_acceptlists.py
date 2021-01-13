import sys, quiet, math,string
from numpy import zeros
from numpy import fft,average,log10,sqrt

# ---------------------------
def accept_to_list(acceptlist):
	n=acceptlist.acount()
	listlist=zeros([n,19,4,5])
	for i in range(0,n):
	    for module in range(0,19):
		for diode in range(0,4):
		   if acceptlist.scans.items()[i][1][module][diode]:
			thing=0
			listlist[i,module,diode,:] = [acceptlist.scans.items()[i][0][0],acceptlist.scans.items()[i][0][1],module,diode,thing]
		   else:
			thing=1
			listlist[i,module,diode,:] = [acceptlist.scans.items()[i][0][0],acceptlist.scans.items()[i][0][1],module,diode,thing]
	return listlist

# --------------------------
def fast_iterate_accept_accept(A,B,accepts_common,rejects_common,accepts_A,accepts_B):
   Alist=accept_to_list(A)
  # Blist=accept_to_list(B)

   nA=len(Alist[:,0,0,0])
   for i in range(0,nA):
      print 'on run: ', A.scans.items()[i][0]
      for module in range(0,19):
         for diode in range(0,4):
            if Alist[i,module,diode,4] == 0:
                try:
                   if B.isIncluded(Alist[i,module,diode,0],Alist[i,module,diode,1],Alist[i,module,diode,2],Alist[i,module,diode,3]):
                      accepts_common.write('%d %d %d %d \n' % (Alist[i,module,diode,0],Alist[i,module,diode,1],Alist[i,module,diode,2],Alist[i,module,diode,3]))
                   else:
                      accepts_A.write('%d %d %d %d \n' % (Alist[i,module,diode,0],Alist[i,module,diode,1],Alist[i,module,diode,2],Alist[i,module,diode,3]))
                except:
		   accepts_A.write('%d %d %d %d \n' % (Alist[i,module,diode,0],Alist[i,module,diode,1],Alist[i,module,diode,2],Alist[i,module,diode,3]))
            else:
		try:
                   if B.isIncluded(Alist[i,module,diode,0],Alist[i,module,diode,1],Alist[i,module,diode,2],Alist[i,module,diode,3]):
                      accepts_B.write('%d %d %d %d \n' % (Alist[i,module,diode,0],Alist[i,module,diode,1],Alist[i,module,diode,2],Alist[i,module,diode,3]))
                   else:
                      rejects_common.write('%d %d %d %d \n' % (Alist[i,module,diode,0],Alist[i,module,diode,1],Alist[i,module,diode,2],Alist[i,module,diode,3]))
		except:
		      rejects_common.write('%d %d %d %d \n' % (Alist[i,module,diode,0],Alist[i,module,diode,1],Alist[i,module,diode,2],Alist[i,module,diode,3]))


# --------------------------
def fast_iterate_accept_reject(A,B,accepts_common,rejects_common,accepts_A,accepts_B):
   Alist=accept_to_list(A)
   #Blist=accept_to_list(B)

   # THIS ASSUMES THE SECOND LIST IS A REJECT LIST!!!!!!!!
   nA=len(Alist[:,0,0,0])
   for i in range(0,nA):
      print 'on run: ', A.scans.items()[i][0]
      for module in range(0,19):
         for diode in range(0,4):
	    if Alist[i,module,diode,4] == 0:
		if B.isIncluded(Alist[i,module,diode,0],Alist[i,module,diode,1],Alist[i,module,diode,2],Alist[i,module,diode,3]):
		   accepts_common.write('%d %d %d %d \n' % (Alist[i,module,diode,0],Alist[i,module,diode,1],Alist[i,module,diode,2],Alist[i,module,diode,3]))
	  	else:
		   accepts_A.write('%d %d %d %d \n' % (Alist[i,module,diode,0],Alist[i,module,diode,1],Alist[i,module,diode,2],Alist[i,module,diode,3]))		
	    else:
		if B.isIncluded(Alist[i,module,diode,0],Alist[i,module,diode,1],Alist[i,module,diode,2],Alist[i,module,diode,3]):
                   accepts_B.write('%d %d %d %d \n' % (Alist[i,module,diode,0],Alist[i,module,diode,1],Alist[i,module,diode,2],Alist[i,module,diode,3]))
		else:
                   rejects_common.write('%d %d %d %d \n' % (Alist[i,module,diode,0],Alist[i,module,diode,1],Alist[i,module,diode,2],Alist[i,module,diode,3]))


# --------------------------
def slow_iterate(A,B,accepts_common,rejects_common,accepts_A,accepts_B):
 
  nA = A.acount()
  nB = B.acount()

  for i in range(0,nA):
     print 'on run:', A.scans.items()[i][0]
     for module in range(0,19):
        for diode in range(0,4):
           counterA = 0
           for j in range(0,nB):
              if A.scans.items()[i][0] == B.scans.items()[j][0] and A.scans.items()[i][1][module][diode] == B.scans.items()[j][1][module][diode]:
                 counterA=counterA+1
                 if A.scans.items()[i][1][module][diode]: # then it was common between lists, print out in common accept file
                   accepts_common.write('%d %d %d %d \n' % (A.scans.items()[i][0][0], A.scans.items()[i][0][1],module, diode))
                 else:
                   rejects_common.write('%d %d %d %d \n' % (A.scans.items()[i][0][0], A.scans.items()[i][0][1],module, diode))
           if counterA == 0:
              accepts_A.write('%d %d %d %d \n' % (A.scans.items()[i][0][0], A.scans.items()[i][0][1],module, diode))

  for j in range(0,nB):
     for module in range(0,19):
        for diode in range(0,4):
           counterB = 0
           for i in range(0,nA):
              if A.scans.items()[i][0] == B.scans.items()[j][0] and A.scans.items()[i][1][module][diode] == B.scans.items()[j][1][module][diode]:
                 counterB=counterB+1
           if counterB == 0:
              accepts_B.write('%d %d %d %d \n' % (B.scans.items()[j][0][0], B.scans.items()[j][0][1],module, diode))



###########################################################################
######### MAIN PART IS HERE ###############################################
###########################################################################

if __name__ =="__main__":

  # main part of the code starts

  accepts_common = open(string.join([sys.argv[1],'_vs_',sys.argv[2],'_accepts_common.txt'],''),"w+")
  rejects_common = open(string.join([sys.argv[1],'_vs_',sys.argv[2],'_rejects_common.txt'],''),"w+")
  accepts_A = open(string.join([sys.argv[1],'_vs_',sys.argv[2],'_accepts_onlylist1.txt'],''),"w+")
  accepts_B = open(string.join([sys.argv[1],'_vs_',sys.argv[2],'_accepts_onlylist2.txt'],''),"w+")

  A=quiet.AcceptedScanList(sys.argv[1])
  B=quiet.AcceptedScanList(sys.argv[2])
  #B=quiet.acceptedListFromRunList(sys.argv[2])
  
  #slow_iterate(A,B,accepts_common,rejects_common,accepts_A,accepts_B)
  fast_iterate_accept_accept(A,B,accepts_common,rejects_common,accepts_A,accepts_B)
 
  accepts_common.close()
  rejects_common.close()
  accepts_A.close()
  accepts_B.close()

############################################################################
############ END ###########################################################
############################################################################
