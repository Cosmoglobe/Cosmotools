#!/usr/bin/python

"""
Script to run through Chicago gains from QUFS SQL and:
1. Convert the CES indexing into Oslo format
2. Remove any CESs we don't use in Oslo
3. Attach MJDs from runlist
4. Write out to a new file

Run as
convertGains.py Chicagogains_input.sorted.txt chicago-gains-converted.txt
Input: Chicagogains_input.sorted.txt comes from, in SQL:

select q19_gain_model_each_runs.run_id,
q19_gain_model_each_runs.run_subid, q19_gain_model_each_runs.module,
q19_gain_model_each_runs.diode, q19_gain_model_each_runs.G from
q19_gain_model_each_runs where q19_gain_model_each_runs.module < 16;

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

if __name__ == "__main__":

    chicf=sys.argv[1]
    newf=sys.argv[2]

    # Fetch all the CESs we have
    rundict,accdict=secateur.GenerateCMBRunlists()
    acclist=[]
    for patch in accdict.keys():
        acclist+=accdict[patch].scans.keys()
    acclist.sort()

    # Now build dictionary of patches keyed by each scan ID
    patchdict={}
    for patch in accdict.keys():
        for (scan,seg) in accdict[patch].scans.keys():
            patchdict[(scan,seg)]=patch

    #print_dict(patchdict)

    f = open(chicf,'r')
    f2= open(newf,'w')

    print 'Reading (slow!) %s' % chicf

    i=1
    sample=50000
    for line in f:
        #scan,seg,module,diode,gain=line.split() # TT style
	if line.startswith('#'): continue
	scan,seg,module,strdiode,gain=line.split()
	diode=secateur.DiodesLookup[strdiode]
        [(newscan,newseg)] = ces_translate([(int(scan),int(seg))])
        if newscan != -1 and (newscan,newseg) in acclist:
            #chicagogains[(newscan,newseg,int(module),int(diode))]=float(gain)
            patch=patchdict[(int(newscan),int(newseg))]
            x,y,z,mjdstart,mjdstop \
                =secateur.FetchPointing(newscan,newseg,rundict[patch],patch)
            line2= '%4i %1i %s %14.8f %14.8f %2i %2i %7.5f\n' \
                   % (newscan,newseg,patch,mjdstart,mjdstop,\
                      int(module),int(diode),float(gain))
            #print line2
            f2.write(line2)
            f2.flush()
        if (i % sample == 0): print 'Reached %i' % i
        i+=1

    f.close()
    f2.close()

    sys.exit(0)
