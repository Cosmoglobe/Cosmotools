#!/usr/bin/python

import secateur
from ces_translate import ces_translate
import quiet,copy,numpy,sys

if __name__=='__main__':

#    A=numpy.loadtxt(chicf)
#    n=len(A[:,0])
#    print 'Read    %s' % chicf
    #na=len(A[0,:])
#    for i in range(0,n):
#       [(newscan,newseg)] = ces_translate.ces_translate([(A[i,0],A[i,1])])
#       if newscan != -1:
#           chicagogains[(newscan,newseg,A[i,2],A[i,3])]=A[i,4]
#       else:
#           #print 'not here,fail',newscan,newseg

    # Set up default runlist and acceptlist if not user-supplied
    runlist=None
    patch='patch_6a'
    acceptlist=None
    if runlist is None:
        rundict,acceptdict=secateur.GenerateCMBRunlists()
        runlist=rundict[patch]
        if acceptlist is None:
            acceptlist=acceptdict[patch]

    # Remove any scans before #125:
    earlyscans=range(125)
    for (scan,seg) in acceptlist.scans.keys():
        if scan in earlyscans:
            acceptlist.removeScan(scan,seg)

    # For the acceptlist, get chicago gains
    readchicago=True
    if readchicago:
        chicf='Chicagogains_input.txt'
        #chicf='10000.txt'
        print 'Reading (slow!) %s' % chicf
        chicagogains={}
        file=open(chicf)
        i=1
        sample=20000
        for line in file:
            scan,seg,module,diode,gain=line.split()
            [(newscan,newseg)] = ces_translate([(int(scan),int(seg))])
            if newscan != -1 and (newscan,newseg) in acceptlist.scans.keys():
                chicagogains[(newscan,newseg,int(module),int(diode))]=float(gain)
            if (i % sample == 0):
                print 'Reached %i' % i
            i+=1

        print 'Read    %s' % chicf
        file.close()
    else:
        chicagogains={}

    # For the acceptlist, get MJDs at middles of CESs
    cestime={}
    for scan,seg in sorted(acceptlist.scans.keys()):
       az,el,dk,mjdstart,mjdstop=\
             secateur.FetchPointing(scan,seg,runlist,patch)
       cestime[(scan,seg)] = (mjdstart+mjdstop)/2.0
       # print '%i.%i el %2.1f' % (scan,seg,el[(scan,seg)])

    # Fetch Rob's gains
    mod={}
    diodelist=[(17,0),(17,3),(18,0),(18,3)]
    # read in module, diode gains:
    for (module,diode) in sorted(diodelist):
        filename = 'gainTT-%i-%i.txt' % (module,diode)
        print 'Reading %s' % filename
        # this is rob's gain in <mjd mjd mod dio gain>
        mod[(module,diode)]=numpy.loadtxt(filename)
        #m=len(mod[(module,diode)])

    # Open the output file
    outfile='output.txt'
    out=open(outfile,'w')

#    for (scan,seg) in sorted(acceptlist.scans.keys()):
#        for module in [17,18]:
#            for diode in [0,3]:
#                if (scan,seg,module,diode) in columbiagains.keys():
#                    print 'col',scan,seg,module,diode
#                if (scan,seg,module,diode) in chicagogains.keys():
#                    print 'chic',scan,seg,module,diode

    # Determine Rob's gains for the CESs in the acceptlist
    print 'Matching gains (slow)'
    columbiagains={}
    verbose=False
    #verbose=True
#    print mod[(17,3)][:,0] # This is the list of dates
    for (scan,seg) in sorted(acceptlist.scans.keys()):
        #if verbose: print scan,seg,cestime[(scan,seg)]
        dates=mod[(module,diode)][:,0]
        nearest=sorted(dates, key=lambda x: abs(x-cestime[(scan,seg)]))[0]
        index=numpy.where(mod[(module,diode)][:,0]==nearest)
        #print index,nearest,mod[(module,diode)][index,-1]
        for module in [17,18]:
            for diode in [0,3]:
                columbiagains[(scan,seg,module,diode)]=\
                    mod[(module,diode)][index,-1]
    # Remove any scans from acceptlist that don't have Chicago gains:
                if (scan,seg,module,diode) in columbiagains.keys() and \
                   (scan,seg,module,diode) in chicagogains.keys():
                    line='%i %i %i %i %12.8f %12.8f\n' % (scan,seg,module,diode,\
                        chicagogains[(scan,seg,module,diode)],\
                        columbiagains[(scan,seg,module,diode)])
                    if verbose: print line
                    out.write(line)
                else:
                    pass
                    #print 'Problem with ',scan,seg
        out.flush()

    out.close()

    sys.exit(0)
