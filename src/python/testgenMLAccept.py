#!/usr/bin/python

"""
Generate test acceptlist for ML null tests
"""

if __name__ == "__main__":

    import os,sys,copy
    import quiet,secateur

    print '***This is for testing only'

    # Write the acceptlist to this file
    outfile='mltestaccept.test'

    # Specify the cuts directory
    cutsdir='/data4/quiet/nulltests/cuts/'

    # This is the patch_6a acceptlist after Laura's weather cut
    keepfs=['weather_acceptlist_asrun_300.txt',\
            'weather_acceptlist_asrun_30.txt',\
            'weather_acceptlist_asrun_60.txt']
    #weather_keepf_6a='weather_cutlist_100101.txt'

    # Want the intersection of these three weather acceptlists
    acceptd={}
    for weather_keepf_6a in keepfs:
        acceptd[weather_keepf_6a]=\
            quiet.AcceptedScanList(os.path.join(cutsdir,weather_keepf_6a))

    a1=acceptd[keepfs[0]]
    a2=acceptd[keepfs[1]]
    a3=acceptd[keepfs[2]]
    print 'Merging weather cuts'
    for (scan,seg) in a1.scans.keys():
        if (scan,seg) in a2.scans.keys() and (scan,seg) in a3.scans.keys():
            pass
        else:
            a1.removeScan(scan,seg)

    acceptlist=copy.deepcopy(a1)

    print ' acceptlist has %i CESs with diodes %i in %i out' \
          % (acceptlist.acount(),acceptlist.countIncludedScanDiodes(),\
             acceptlist.countExcludedScanDiodes())

    # Fetch the runlist (only because needed by quiet.pruneAcceptList)
    patch='patch_6a'
    rundict,acceptdict=secateur.GenerateCMBRunlists()
    runlist=rundict[patch]

    # Apply the static cuts to this acceptlist:
    deadfiles=['dead.cut','static.cut']
    for deadfile in deadfiles:
        print 'Applying cut from %s' % (deadfile)
        dead=quiet.CutList(os.path.join(cutsdir,deadfile))
        dead.pruneAcceptList(runlist,acceptlist)
        print ' acceptlist has %i CESs with diodes %i in %i out' \
              % (acceptlist.acount(),acceptlist.countIncludedScanDiodes(),\
                 acceptlist.countExcludedScanDiodes())

    # Now cut out CES module 0 for which 60 deg <= el <= 65 deg
    verbose=False
    modules=[0]
    elmin=60.0
    elmax=65.0
    print 'Cutting out modules %i, %2.1f <= el <= %2.1f deg' % (modules[0],elmin,elmax)
    acceptlist=secateur.CutElRange(acceptlist,modules,elmin,elmax,runlist,verbose)
    print ' acceptlist has %i CESs with diodes %i in %i out' \
          % (acceptlist.acount(),acceptlist.countIncludedScanDiodes(),\
             acceptlist.countExcludedScanDiodes())

    # Now cut out any CESs for which scan length < 10 minutes
    verbose=False
    tmin=10.0
    print 'Cutting out CESs of duration <= %2.1f min' % tmin
    acceptlist=secateur.CutScanLength(acceptlist,tmin,runlist,verbose)
    print ' acceptlist has %i CESs with diodes %i in %i out' \
          % (acceptlist.acount(),acceptlist.countIncludedScanDiodes(),\
             acceptlist.countExcludedScanDiodes())

    # Now cut out any CESs that have been mispointed
    verbose=False
    print 'Cutting out CESs that are mispointed'
    acceptlist=secateur.CutMispointedScans(acceptlist,runlist,verbose)
    print ' acceptlist has %i CESs with diodes %i in %i out' \
          % (acceptlist.acount(),acceptlist.countIncludedScanDiodes(),\
             acceptlist.countExcludedScanDiodes())

    # Show summary stats
    #acceptlist.removeEmptyScans()
    print 'To summarize: from %i,' % acceptlist.acount()
    acceptlist.countSummary()

    # You get one life
    f=os.path.join(cutsdir,outfile)
    if os.path.exists(f):
        os.rename(f,f+'.old')
        print '(renamed %s -> .old)' % f

    acceptlist.toFile(f)

    print 'acceptlist written to %s' % f

    sys.exit(0)
