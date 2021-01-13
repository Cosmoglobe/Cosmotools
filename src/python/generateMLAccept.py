#!/usr/bin/python

"""
Generate the acceptlist for ML null tests
This version is for ALL patches (work in progress)
"""

#-------------------------------------------------------------------

def finalize(acceptlist,cutsdir,outfile,header=None):

    "Write the acceptlist to outfile in cutsdir, and exit"

    import sys,os

    acceptlist.printSummary()

    # Output acceptlist to file
    f=os.path.join(cutsdir,outfile)
    if os.path.exists(f):
        os.rename(f,'%s.old' %f)
        print '(renamed %s -> .old)' % f
    acceptlist.toFile(f,header=header)
    print 'acceptlist written to %s' % f

    # Show summary stats
    #print 'To summarize: from %i,' % acceptlist.acount()
    #acceptlist.countSummary()

#    ndiodesfinal=acceptlist.printSummary()

#    print '\nCES-diode summary'
#    print '====================================='
#    print 'Before cuts:      %i                 ' % ndiodes0
#    print 'After fixed cuts: %i (100 per cent)  ' % ndiodes1
#    print 'After all cuts:   %i (%3.1f per cent)' \
#          % (ndiodesfinal,100.0*float(ndiodesfinal)/float(ndiodes1))
#    if doweathercut: print 'weather removed %3.1f per cent' % weatherfrac

    sys.exit(0)

    return

#-------------------------------------------------------------------

if __name__ == '__main__':

    import sys

    #if len(sys.argv) != 2:
    #    print 'Usage: generateMLAccept.py patch (2a|4a|6a|7b)'
    #    sys.exit(1)

    import os,copy
    import quiet,secateur
    from optparse import OptionParser

    print
    # Set up options
    usage = "./generateMLAccept.py --patch (2a|4a|6a|7b)' \n"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="make lots of noise")

    parser.add_option("-p", "--patch", dest="patch",\
           help="specify QUIET CMB patch", type="string",\
           metavar='(2a|4a|6a|7b)')

    parser.add_option("-s", "--static", action="store_true", \
                      dest="dostaticcuts",help="apply static cuts")

    parser.add_option("-d", "--duration",action="store_true",\
                      dest="dodurationcut",help="apply scan-duration cuts")

    parser.add_option("-b", "--typeb", action="store_true", \
                      dest="dotypebcut",help="apply type-B cut")

    parser.add_option("-w", "--weather", action="store_true",\
                      dest="doweathercut",help="apply weather cut")

    parser.add_option("-a", "--phsw", action="store_true",\
                      dest="dophswcut",help="apply phase-switch cuts")

    parser.add_option("-l", "--lunch", action="store_true", dest="dosuncut",
                      help="apply sidelobe cuts")

    parser.add_option("-y", "--lunch2", action="store_true", dest="dosuncut2",
                      help="apply sidelobe cuts (v31 or v41)")

    parser.add_option("-f", "--fft", action="store_true", dest="dofftcuts",
                      help="apply tod2map cuts")

    parser.add_option("-z", "--header", action="store_true", dest="doheader",
                      help="write a header")

    parser.add_option("-q", "--quit", dest="quit",\
           help="output acceptlist and exit after this cut", type="string",\
           metavar='(static|duration|typeb|weather|sun|sun2|phsw|fft)')

    parser.add_option("-o", "--outfile", dest="outfile",\
           help="specify outfile", type="string",\
           metavar='e.g. mlaccept2a-100418.txt')

    # Assert defaults
    #parser.set_defaults(\
    #    patch=None,verbose=False,dostaticcuts=True,dodurationcut=True,\
    #    dotypebcut=True,doweathercut=True,dosuncut=True,dophswcut=False,\
    #    doheader=False,quit=None)
    parser.set_defaults(\
        patch=None,verbose=False,dostaticcuts=False,dodurationcut=False,\
        dotypebcut=False,doweathercut=False,dosuncut=False,dophswcut=False,\
        dofftcuts=False,doheader=False,quit=None,outfile=None,dosuncut2=False)
    # Parse command line
    (options, args) = parser.parse_args()

    patch=options.patch
    verbose=options.verbose
    dostaticcuts=options.dostaticcuts
    dodurationcut=options.dodurationcut
    dotypebcut=options.dotypebcut
    doweathercut=options.doweathercut
    dosuncut=options.dosuncut
    dosuncut2=options.dosuncut2
    dophswcut=options.dophswcut
    dofftcuts=options.dofftcuts
    doheader=options.doheader
    quit=options.quit
    outfile=options.outfile

    # Set up the patch
    assert(patch is not None), 'patch is required!'
    patch = 'patch_%s' % patch
    print 'patch is %s' % patch[-2:]
    assert(patch in secateur.CMBPatches), '%s not available!' % patch

    # Check the quit option if it's there
    if quit is not None:
        assert(quit in ['static','duration','typeb','weather','sun','phsw',\
                        'fft','sun2']), 'invalid quit option %s' % quit

    # Sun cuts are mutually exclusive
    assert(not(dosuncut and dosuncut2)), "You can't have your cake and eat it!"

    # Set up the header if it's required
    print 'Options are %s' % str(options)
    if doheader:
        import re, glob
        from datetime import date,datetime
        user = os.getenv('USER')
        host = os.getenv('HOSTNAME').split('.')[0]
        host = re.sub('login','titan',host) # Switch login -> titan
        date = datetime.utcnow().strftime("%Y-%m-%d %A %I:%M%p UT")
        header='# %s\n# run by %s on %s\n# at %s\n'\
                % (str(options),user,host,date)
    else:
        header=None

    # Specify the cuts directory
    cutsdir='/data4/quiet/nulltests/cuts/'

    # Write the acceptlist to this file
    if outfile is None: outfile='mlaccept%s.txt' % patch[-2:]

    # Fetch the runlist and acceptlist for this patch
    rundict,acceptdict=secateur.GenerateCMBRunlists()
    runlist=rundict[patch]
    acceptlist=acceptdict[patch]
    ndiodes0=acceptlist.printSummary()

    if doheader: header='%s# initially: %i diodes\n' % (header,acceptlist.printSummary(silent=True))

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

    if dostaticcuts:
        removeTT=True
        # Remove the TT diodes
        # (Eventually try to replace this by modifying acceptlist methods)
        # (so that the TT diodes are still there but just ignored)
        if removeTT:
            print 'Removing TT diodes'
            for module in range(quiet.N_MODULE_MAX-2,quiet.N_MODULE_MAX):
                acceptlist.excludeModule(module)
            dummy=acceptlist.printSummary()
        else: print '**TT diodes not removed'

        # Remove scans before #125
        print 'Removing any CESs before phase switches were balanced'
        earlyscans=range(125)
        for (scan,seg) in acceptlist.scans.keys():
            if scan in earlyscans:
                acceptlist.removeScan(scan,seg)
        dummy=acceptlist.printSummary()

        # Now cut out any mislabelled CESs
        print 'Cutting out mislabelled CESs'
        acceptlist=secateur.CutMislabelledScans(\
            acceptlist=acceptlist,verbose=verbose,patch=patch)
        dummy=acceptlist.printSummary()

        # Remove dead modules/diodes
        print 'Removing RQ04_Q1, RQ08_U2 and RQ16'
        dead=[(4,'Q1'),(8,'U2')]
        for (mod,diode) in dead:
            dinum=secateur.DiodesLookup[diode]
            acceptlist.excludeDiode(mod,dinum)
        mod=16
        acceptlist.excludeModule(mod)
        dummy=acceptlist.printSummary()

        # 'Removing RQ12_Q2 for CESs before 766 inclusive'
        print 'Removing noisy RQ12_Q2 for CESs before 766 inclusive'
        (m,d)=(12,secateur.DiodesLookup['Q2'])
        for (scan,seg) in acceptlist.scans.keys():
            if scan <= 766: acceptlist.excludeSegmentDiode(scan,seg,m,d)
            else: pass
        dummy=acceptlist.printSummary()

        # Engineering cuts
        print 'Removing CESs under static engineering cuts'
        # Supply patch, but not really used when acceptlist passed too
        acceptlist=secateur.CutStaticCESs(acceptlist=acceptlist,\
                                          verbose=verbose,patch=patch)
        ndiodes1=acceptlist.printSummary()

    if doheader: header='%s# post-static: %i diodes\n' % (header,acceptlist.printSummary(silent=True))

    if quit=='static':
        finalize(acceptlist,cutsdir,outfile,header=header)

#### Slightly flexible cuts

    # Now cut out any CESs for which scan length < 10 min or < 1000 sec
    if dodurationcut:
        #tmin=10.0
        tmin=1000.0/60.0 # Changed scan length to 1000 sec
        print 'Cutting out CESs of duration <= %2.1f min' % tmin
        acceptlist=secateur.CutScanLength(\
            acceptlist=acceptlist,tmin=tmin,runlist=runlist,\
            verbose=verbose,patch=patch)
        dummy=acceptlist.printSummary()
    else: print '**CES-duration cut not applied'

    if quit=='duration':
        finalize(acceptlist,cutsdir,outfile,header=header)

    # Now cut out CESs for which type-B parameter fits are poor
    if dotypebcut:
        print 'Cutting out CESs for which type-B Chi^2 > 10.0'
        acceptlist=secateur.CutTypeBBadness(\
            acceptlist=acceptlist,verbose=verbose,patch=patch)
    else: print '**Type-B cut not applied'
    npreweather=acceptlist.printSummary()

    if doheader: header='%s# post-typeb: %i diodes\n' % (header,acceptlist.printSummary(silent=True))

    if quit=='typeb':
        finalize(acceptlist,cutsdir,outfile,header=header)

#### Negotiable cuts

    # Weather cut
    if doweathercut:
        # Some of these have dependencies below
        #cuttype=1  # 1 30,60,300 1-sigma
        #cuttype=2  # 2 whole CESs (weather_cut_100114.txt)
        #cuttype=3  # 3 30,60 1.5-sigma 
        #cuttype=4  # 4 whole CESs (TTrejectlist_100114.txt)
        #cuttype=5  # 5 30,60 1-sigma
        #cuttype=6  # 6 30,60 2-sigma (asymmetric)
	#cuttype=7  # 7 10, 30 2-sigma (asymmetric)
        #cuttype=8  # 8 10,30 5-sigma (asymm with Gaussian)
        #cuttype=9  # 9 10,30 5-sigma (symm Gaussian)
        cuttype=10 # 9 10,30 5-sigma (symm Gaussian) but half-heaviside, by patch
        print 'Merging and applying weather cuts (type %i)' % cuttype

        acceptlist=secateur.CutFailedWeatherDiodes(\
            acceptlist=acceptlist,cuttype=cuttype,patch=patch,verbose=verbose)
        npostweather=acceptlist.printSummary()

        weatherfrac=100.0*float(npreweather-npostweather)/float(npreweather)
        print 'Weather cut has removed %3.1f per cent of data' \
              % weatherfrac
    else: print '**Weather cut not applied'

    if doheader: header='%s# post-weather: %i diodes\n' % (header,acceptlist.printSummary(silent=True))

    if quit=='weather':
        finalize(acceptlist,cutsdir,outfile,header=header)

    # Sun-in-sidelobes cut
    if dosuncut:
        modules=[0,1,2]
        print 'Applying sun-in-sidelobe cut (modules %s)' % str(modules)
        acceptlist=secateur.CutSunspikes(\
            acceptlist=acceptlist,runlist=runlist,verbose=verbose,\
            patch=patch,modules=modules)
        dummy=acceptlist.printSummary()
    else: print '**Sun-in-sidelobe cut not applied (algorithm version)'

    if quit=='sun':
        finalize(acceptlist,cutsdir,outfile,header=header)

    # Sun-in-sidelobes cut (vers 31 or 41)
    #sunvers=31
    sunvers=41
    if dosuncut2:
        verbose=False
        print 'Applying sun-in-sidelobe cut (**version %i)' % sunvers
        acceptlist=secateur.CutSunspikes(\
            acceptlist=acceptlist,runlist=runlist,verbose=verbose,\
            patch=patch,usealgorithm=False,sunvers=sunvers)
        dummy=acceptlist.printSummary()
    else: print '**Sun-in-sidelobe cut (v%i) not applied' % sunvers

    if quit=='sun2':
        finalize(acceptlist,cutsdir,outfile,header=header)

    # Phase-switch humidity cut
    if dophswcut:
        print 'Applying phase-switch cut'
        acceptlist=secateur.CutPhaseSwitchHumidity(\
            acceptlist=acceptlist,verbose=verbose,patch=patch)
        dummy=acceptlist.printSummary()
    else: print '**Phase-switch cut not applied'

    if quit=='phsw':
        finalize(acceptlist,cutsdir,outfile,header=header)

    # tod2map cuts
    if dofftcuts:
        print 'Applying tod2map (FFT and TOD) cuts'
        t2mdir='/data4/quiet/jonathaz4/2010_jan/maps/maps_nocuts'
        t2mf='%s_reports/accepted_autocut_%s.txt' % (patch[-2:],patch[-2:])
        # BEGIN PATCH
        print '*** TEMPORARY PATCH FOR FFT FILES (Chi^2 < 10)'
        t2mdir='/data4/quiet/jonathaz4'
        t2mf='maps_nocuts_100317_%s/accepted_autocut_set001.txt' % patch[-2:]
        # END PATCH
        fftf=os.path.join(t2mdir,t2mf)
        print '** (postprocessing using %s)' % fftf
        at2m=quiet.AcceptedScanList(fftf)
        acceptlist=secateur.PostProcessAccept(\
            ain=acceptlist,aout=at2m,verbose=verbose)
        dummy=acceptlist.printSummary()
    else: print '**FFT cuts not applied'

    if quit=='fft':
        finalize(acceptlist,cutsdir,outfile,header=header)

    # TEMPORARY RQ12_Q2
    #dinum=secateur.DiodesLookup['Q2']
    #modnum=12
    #n=0
    #m=0
    #tot=0
    #split=765
    #for (scan,seg) in sorted(acceptlist.scans.keys()):
    #    if scan <= split: continue
    #    if acceptlist.isIncluded(scan,seg,modnum,dinum):
    #        n+=1
    #    if acceptlist.isExcluded(scan,seg,modnum,dinum):
    #        m+=1
    #    tot+=1
    #print 'RQ12_Q2: %i in + %i out = %i' % (n,m,tot)

    if doheader: header='%s# output: %i diodes\n' % (header,acceptlist.printSummary(silent=True))

    # Write out acceptlist and exit
    finalize(acceptlist,cutsdir,outfile,header=header)

    sys.exit(0)

#-----------------------------------------------------------------
