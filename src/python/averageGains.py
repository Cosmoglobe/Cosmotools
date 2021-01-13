#!/usr/bin/python

"""
Script to run through Columbia gains and, for all the CESs we have:
1. Identify the start MJD in the gain file
2. Start averaging from the start MJD
3. Stop averaging at the stop MJD
4. Write out to a new file

Run as
convertGains.py qband_pol_gain_091124.txt columbia-gains-converted.txt
Input: qband_pol_gain_091124.txt from /data4/quiet/calib_data/rob/

"""

import sys
from ces_translate import ces_translate
import secateur

def reverse(a_to_b):
        b_to_a = {}
        for a in a_to_b:
            b = a_to_b[a]
            b_to_a[b] = a
            return b_to_a

def print_dict(dict):
    for key in dict:
        print key, dict[key]

def online_variance(data):
    """
    from http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    """
	
    n = 0
    mean = 0
    M2 = 0

    for x in data:
        n = n + 1
	delta = x - mean
	mean = mean + delta/n
	M2 = M2 + delta*(x - mean) # This expression uses the new value of mean

    variance_n = M2/n
    variance = M2/(n - 1)
    return variance


if __name__ == "__main__":

    colf=sys.argv[1]
    newf=sys.argv[2]

    # Fetch all the CESs we have
    rundict,accdict=secateur.GenerateCMBRunlists()
    acclist=[]
    for patch in accdict.keys():
        acclist+=accdict[patch].scans.keys()
    acclist.sort()

    # Now build dictionary of patches keyed by each scan ID
    # Also build dictionary of mjd ranges keyed by scan ID
    patchdict={}
    mjddict={}
    for patch in accdict.keys():
        for (scan,seg) in accdict[patch].scans.keys():
            patchdict[(scan,seg)]=patch
            x,y,z,mjdstart,mjdstop \
                =secateur.FetchPointing(scan,seg,rundict[patch],patch)
	    mjddict[(scan,seg)]=[mjdstart,mjdstop]

    #print_dict(patchdict)

    f = open(colf,'r')
    f2= open(newf,'w')

    print 'Reading (slow!) %s' % colf

    i=1
    sample=50000
    # gain 54750.082488700   54750.082488700    -1  -1   0   0   2.01455949e+00
    # In the file, each diode is one after the other
    nscan=0
    nmodule=0
    ndiode=0
    lastnmodule=0
    lastndiode=0
    lastscan=(0,0)
    mean=0
    n=0
    startcounting=False
    stopcounting=False
    for line in f:
        word,mjd,x,y,z,nmodule,ndiode,gain=line.split()
	if word != 'gain': continue # skip non-gains
	# Reset to the beginning of the scan list for each change of mod,diode
	if (lastnmodule,lastndiode) != (int(nmodule),int(ndiode)):
		nscan=0
	# Assumes acclist has been sorted by scan:
	currentscan=acclist[nscan]
	# Finish this scan calculation and start again for the next CES
	if currentscan != lastscan:
	    line2='%s %s %2i %2i %f' \
		   % (str(currentscan),str(mjddict[currentscan]),\
		      int(nmodule),int(ndiode),mean)
	    f2.write(line2)
	    raw_input('press key')
	    n=0
	    mean=0
	# If MJD is in range for this scan
	if float(mjd) > mjddict[currentscan][0]: startcounting = True
	if float(mjd) < mjddict[currentscan][1]: stopcounting = False
	if startcounting and not stopcounting:
	    n+=1
	    delta = float(gain) - mean
	    mean = mean + delta/float(n)
	    print mjd,nmodule,ndiode,currentscan,mean,mjddict[currentscan]
	if stopcounting:
	    nscan+=1

	lastnmodule=nmodule
	lastndiode=ndiode
	lastscan=currentscan
	

#            line2= '%4i %1i %s %14.8f %14.8f %2i %2i %7.5f\n' \
#                   % (newscan,newseg,patch,mjdstart,mjdstop,\
#                      int(module),int(diode),float(gain))
#            #print line2
#            f2.write(line2)
#            f2.flush()
#        if (i % sample == 0): print 'Reached %i' % i
#        i+=1

    f.close()
    f2.close()

    sys.exit(0)
