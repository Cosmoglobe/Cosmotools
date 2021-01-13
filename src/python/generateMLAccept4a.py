#!/usr/bin/python

"""
Generate the acceptlist for ML null tests
THIS VERSION IS FOR PATCH 4a
"""

#dry=True
dry=False

if __name__ == "__main__":

    import os,sys,copy
    import quiet,secateur

    # Write the acceptlist to this file
    outfile='mlaccept4a.txt'

    if dry:
        outfile='dry.txt'
        print 'dry run to %s' % outfile
        verbose=True

    verbose=False

    # Specify the cuts directory
    cutsdir='/data4/quiet/nulltests/cuts/'

    # Fetch the runlist (needed for several cuts)
    patch='patch_4a'
    rundict,acceptdict=secateur.GenerateCMBRunlists()
    runlist=rundict[patch]
    # And set up the acceptlist for this patch
    acceptlist=acceptdict[patch]
    ndiodes0=acceptlist.printSummary()

#### Cuts:

#1. *34128 - are removed by the weather cut (need Laura if want to adjust)
#2. 2571 - are lost to not-very-negotiable engineering-type cuts
#3. *297 - lost to a TT-based scan cut
#4. *2246 - 56 < EL < 66 deg cut on modules 0,1,2
#5. *3113 - Sun < 50 deg or Moon < 50 deg
#6. 0 - scans < 10 mins [already cut out]
#7. 66 - Mispointed scans
#8. 128 - Cut on Type-B goodness-of-fit [Chi^2 > 10; Chicago is thinking about an even tighter cut]


##### Start with non-negotiable cuts

    # Remove the TT diodes
    # Eventually try to replace this by modifying acceptlist methods
    # so that the TT diodes are still there but just ignored
    removeTT=True
    #removeTT=False
    if removeTT:
        print 'Removing TT diodes'
        for (scan,seg) in acceptlist.scans.keys():
            for module in range(quiet.N_MODULE_MAX-2,quiet.N_MODULE_MAX):
                acceptlist.excludeModule(module)
        acceptlist.printSummary()


    # Remove scans before #125
    print 'Removing any CESs before we had balanced the phase switches'
    earlyscans=range(125)
    for (scan,seg) in acceptlist.scans.keys():
        if scan in earlyscans:
            acceptlist.removeScan(scan,seg)
    acceptlist.printSummary()

    # Now cut out any CESs that have been mispointed
    print 'Cutting out CESs that are mispointed'
    acceptlist=secateur.CutMispointedScans(acceptlist,runlist,verbose,\
                                           patch=patch)
    acceptlist.printSummary()

    # Apply the static cuts to this acceptlist:
    # dead.cut
    deadfiles=['dead.cut']#,'static.cut']
    for deadfile in deadfiles:
        print 'Applying cut from %s' % (deadfile)
        dead=quiet.CutList(os.path.join(cutsdir,deadfile))
        dead.pruneAcceptList(runlist,acceptlist)
        acceptlist.printSummary()

    # static.cut
    print 'Removing CESs under static engineering cuts'
    # Supply patch, but not really used when acceptlist passed too
    acceptlist=secateur.CutStaticCESs(acceptlist=acceptlist,\
                                      verbose=verbose,patch=patch)
    ndiodes1=acceptlist.printSummary()

    # Remove RQ4_Q1
    print 'Removing RQ4_Q1'
    mod=4
    dinum=secateur.DiodesLookup['Q1']
    acceptlist.excludeDiode(mod,dinum)
    ndiodes1=acceptlist.printSummary()

    # Remove some more Chicago "bad" CESs
    print 'Removing more Chicago "bad" CESs'
    badlist=[(529,2),(951,1),(953,2),(953,3),(953,4),(954,1),(955,1),(955,2)]
    acceptlist.removeScansIn(badlist)
    ndiodes1=acceptlist.printSummary()

#### Slightly negotiable

    # Now cut out any CESs for which scan length < 10 minutes
    #tmin=10.0
    tmin=1000.0/60.0 # Changed scan length to 1000 sec
    print 'Cutting out CESs of duration <= %2.1f min' % tmin
    acceptlist=secateur.CutScanLength(acceptlist,tmin,runlist,verbose,\
                                      patch=patch)
    acceptlist.printSummary()

    # Now cut out CESs for which type-B parameter fits are poor
    print 'Cutting out CESs for which type-B Chi^2 > 10.0'
    acceptlist=secateur.CutTypeBBadness(acceptlist=acceptlist,\
                                        verbose=verbose)
    npreweather=acceptlist.printSummary()


#### Negotiable cuts

    # Toggle the cuts here
    doweathercut=True
    #doweathercut=False

    dottmapcut=True
    #dottmapcut=False

    doelcut=True
    #doelcut=False

    letsdolunch=True
    #letsdolunch=False

    #ensurebothqandu=True
    ensurebothqandu=False # Causes a problem for printSummary() at moment


    # Apply the weather cut
    if doweathercut:
        # Some of these have dependencies below
        #cuttype=1  # 1 30,60,300 1-sigma
        #cuttype=2  # 2 whole CESs (weather_cut_100114.txt)
        #cuttype=3  # 3 30,60 1.5-sigma 
        #cuttype=4  # 4 whole CESs (TTrejectlist_100114.txt)
        #cuttype=5  # 5 30,60 1-sigma
        #cuttype=6  # 6 30,60 2-sigma (asymmetric)
	#cuttype=7  # 7 10, 30 2-sigma (asymmetric)
        #cuttype=8   # 8 10,30 5-sigma (asymm with Gaussian)
        #cuttype=9   # 9 10,30 5-sigma (symm Gaussian)
        # For 6a and 2a:
        cuttype=10   # 9 10,30 5-sigma (symm Gaussian) but half-heaviside
        print 'Merging and applying weather cuts (type %i)' % cuttype
        acceptlist=secateur.CutFailedWeatherDiodes(\
            acceptlist=acceptlist,cuttype=cuttype,patch=patch)
        npostweather=acceptlist.printSummary()

        weatherfrac=100.0*float(npreweather-npostweather)/float(npreweather)
        print 'Weather cut has removed %3.1f per cent of data' \
              % weatherfrac

    #acceptlist.toFile('test.txt')

    # Apply the TT-map cut
    if dottmapcut:
        print 'Cutting CESs that failed the TT-map cut'
        acceptlist=secateur.CutTTFailedMapScans(acceptlist=acceptlist,\
                                         verbose=verbose,patch=patch)
        acceptlist.printSummary()


    # Now cut out CES module 0,1,2 for which 56 deg <= el <= 66 deg
    if doelcut:
        modules=[0,1,2]
        # modules=[0,1,2,3,4,5,6] #mab1
        elmin=56.0
        elmax=66.0
        print 'Cutting out modules %s %2.1f <= el <= %2.1f deg' \
              % (str(modules),elmin,elmax)
        # [[Should supply patch]]
        acceptlist=\
          secateur.CutElRange(acceptlist,modules,elmin,elmax,runlist,verbose,\
                              patch=patch)
        ndiodesfinal=acceptlist.printSummary()

#    # Now cut out CESs for which Sun or Moon within 50 deg
#    # Cut all modules
#    verbose=True
#    sunthresh=50.0
#    moonthresh=50.0
#    print 'Cutting out whole CESs for which Sun <= %3.1f or Moon <= %3.1f deg' \
#          % (sunthresh,moonthresh)
#    acceptlist=secateur.CutSunMoonDistance(acceptlist=acceptlist,\
#                        runlist=runlist,sunthresh=50.0,\
#                        moonthresh=50.0,verbose=verbose)

    # Remove hand-selected scans left in after removal of Sun/Moon cut
    #dobyeyecut=False
    #if cuttype==1 or cuttype==3:
    dobyeyecut=True

    if dobyeyecut and patch=='patch_6a':
        print 'Removing bad CESs found by inspection of relaxed Sun/Moon cut'
        eyefi='byeye-100114-removethese.txt'
        eyef=open(os.path.join(cutsdir,eyefi),'r')
        for line in eyef:
            strscan,strseg = line.split()
            scan=int(strscan)
            seg=int(strseg)
            if (scan,seg) in sorted(acceptlist.scans.keys()):
                acceptlist.removeScan(scan,seg)
                if verbose: print 'Removed %i.%i' % (scan,seg)
        acceptlist.printSummary()
        eyef.close()

    if cuttype == 6:
        print 'Removing bad CESs found by inspection of asym'
        eyefi='byeye-asym-removethese.txt'
        eyef=open(os.path.join(cutsdir,eyefi),'r')
        for line in eyef:
            strscan,strseg = line.split()
            scan=int(strscan)
            seg=int(strseg)
            if (scan,seg) in sorted(acceptlist.scans.keys()):
                acceptlist.removeScan(scan,seg)
                if verbose: print 'Removed %i.%i' % (scan,seg)
        acceptlist.printSummary()
        eyef.close()

    if cuttype == 10:
        print 'Removing bad CESs found by inspecting sub-maps for cuttype %i, patch %s' \
              % (cuttype,patch)
        if patch=='patch_6a':
            eyefi='byeye-100121a-removethese.txt'
        elif patch=='patch_2a':
            eyefi='byeye-100122b-removethese.txt'
        elif patch=='patch_4a':
            eyefi='byeye-100124a-removethese.txt'
        eyef=open(os.path.join(cutsdir,eyefi),'r')
        #verbose=True
        for line in eyef:
            if line == '' or line == '\n' or line.startswith('#'): continue
            strscan,strseg = line.split()
            scan=int(strscan)
            seg=int(strseg)
            if (scan,seg) in sorted(acceptlist.scans.keys()):
                acceptlist.removeScan(scan,seg)
                if verbose: print 'Removed %i.%i' % (scan,seg)
        acceptlist.printSummary()
        eyef.close()

    # Remove CESs which have the Sun passing through the triple sidelobe
    if letsdolunch:
        verbose=False
        print 'Removing lunch-spike CESs'
        acceptlist=secateur.CutSunspikeScans(acceptlist,verbose=verbose)
        acceptlist.printSummary()

    # Remove any CES-modules that have one of Q or U missing
    if ensurebothqandu:
        verbose=True
        print 'Removing any modules without both some Q and some U'
        acceptlist=secateur.EnsureBothQUinModule(acceptlist,verbose=verbose,\
                                                 patch=patch)
        acceptlist.printSummary()

    # Show summary stats
    print 'To summarize: from %i,' % acceptlist.acount()
    acceptlist.countSummary()

    ndiodesfinal=acceptlist.printSummary()

    print '\nCES-diode summary'
    print '====================================='
    print 'Before cuts:      %i                 ' % ndiodes0
    print 'After fixed cuts: %i (100 per cent)  ' % ndiodes1
    print 'After all cuts:   %i (%3.1f per cent)' \
          % (ndiodesfinal,100.0*float(ndiodesfinal)/float(ndiodes1))
    if doweathercut: print 'weather removed %3.1f per cent' % weatherfrac

    # You get one life
    f=os.path.join(cutsdir,outfile)
    if os.path.exists(f):
        os.rename(f,f+'.old')
        print '(renamed %s -> .old)' % f

    acceptlist.toFile(f)

    print 'acceptlist written to %s' % f

    sys.exit(0)
