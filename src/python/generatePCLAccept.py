#!/usr/bin/python

"""
Generate an acceptlist from a specified combination of Chicago cuts
This version is for ALL patches
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

def subtractAcceptlist(patch=None,cutsdir=None,\
                       acceptlist=None,cuttype=None):

    import os,copy,quiet

    print 'Applying Chicago %s cut' % cuttype

    f='chicago/accept/chicago%s_oslo_%s.acc' % (cuttype,patch[-2:])
    aok=quiet.AcceptedScanList(os.path.join(cutsdir,f))

    retained=copy.deepcopy(acceptlist)

    mods=range(quiet.N_MODULE_MAX)
    diodes=range(quiet.NDIODE_MODULE)

    sok=set(aok.scans.keys())
    for (scan,seg) in sorted(acceptlist.scans.keys()):
        if (scan,seg) not in sok:
            retained.removeScan(scan,seg)
            continue
        else:
            for module in mods:
                for diode in diodes:
                    if aok.isIncluded(scan,seg,module,diode) and \
                           acceptlist.isIncluded(scan,seg,module,diode):
                        retained.includeSegmentDiode(scan,seg,module,diode)
                    else:
                        retained.excludeSegmentDiode(scan,seg,module,diode)

    retained.removeEmptyScans()

    return retained

#-------------------------------------------------------------------

if __name__ == '__main__':

    import sys,os,copy
    import quiet,secateur
    from optparse import OptionParser

    print
    # Set up options
    usage = "./generatePCLAccept.py --patch (2a|4a|6a|7b)' \n"
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

    parser.add_option("-f", "--fft", action="store_true", dest="dofftcuts",
                      help="apply fft cuts")

    parser.add_option("-g", "--glitch", action="store_true", dest="doglitchcut",
                      help="apply time-glitch cut")

    parser.add_option("-z", "--header", action="store_true", dest="doheader",
                      help="write a header")

    parser.add_option("-q", "--quit", dest="quit",\
           help="output acceptlist and exit after this cut", type="string",\
           metavar='(static|scanduration|typeb|weather|sunavoid|phsw)')

    # Assert defaults
    #parser.set_defaults(\
    #    patch=None,verbose=False,dostaticcuts=True,dodurationcut=True,\
    #    dotypebcut=True,doweathercut=True,dosuncut=True,dophswcut=False,\
    #    doheader=False,dofftcuts=True,doglitchcut=True,quit=None)
    parser.set_defaults(\
        patch=None,verbose=False,dostaticcuts=True,dodurationcut=False,\
        dotypebcut=False,doweathercut=False,dosuncut=False,dophswcut=False,\
        doheader=False,dofftcuts=False,doglitchcut=False,quit=None)
    # Parse command line
    (options, args) = parser.parse_args()

    patch=options.patch
    verbose=options.verbose
    dostaticcuts=options.dostaticcuts
    dodurationcut=options.dodurationcut
    dotypebcut=options.dotypebcut
    doweathercut=options.doweathercut
    dosuncut=options.dosuncut
    dophswcut=options.dophswcut
    dofftcuts=options.dofftcuts
    doglitchcut=options.doglitchcut
    doheader=options.doheader
    quit=options.quit

    #Some cuts are unsupported by Chicago
    assert(dophswcut==False), 'Phase-switch cut is not supported by Chicago!'

    # Set up the patch
    assert(patch is not None), 'patch is required!'
    patch = 'patch_%s' % patch
    print 'patch is %s' % patch[-2:]
    assert(patch in secateur.CMBPatches), '%s not available!' % patch

    # Check the quit option if it's there
    if quit is not None:
        assert(quit in ['static','scanduration','typeb',\
                        'weather','sunavoid','phsw','fftcuts','glitch']),\
                    'invalid quit option %s' % quit

    # Set up the header if it's required
    print 'Options are %s' % str(options)
    if doheader:
        header='# %s\n' % str(options)
    else:
        header=None

    # Specify the cuts directory
    cutsdir='/data4/quiet/nulltests/cuts/'

    # Write the acceptlist to this file
    outfile='pclaccept%s.txt' % patch[-2:]

    # Fetch the runlist and acceptlist for this patch
    rundict,acceptdict=secateur.GenerateCMBRunlists()
    acceptlist=acceptdict[patch]
    ndiodes0=acceptlist.printSummary()

##### Start with non-negotiable cuts

    cuttype='static'
    if dostaticcuts:
        acceptlist=subtractAcceptlist(\
            patch=patch,cutsdir=cutsdir,acceptlist=acceptlist,cuttype=cuttype)
        dummy=acceptlist.printSummary()
    else: print '**%s cuts not applied' % cuttype
    if quit==cuttype:
        finalize(acceptlist,cutsdir,outfile,header=header)

    cuttype='scanduration'
    if dodurationcut:
        acceptlist=subtractAcceptlist(\
            patch=patch,cutsdir=cutsdir,acceptlist=acceptlist,cuttype=cuttype)
        dummy=acceptlist.printSummary()
    else: print '**%s cuts not applied' % cuttype
    if quit==cuttype:
        finalize(acceptlist,cutsdir,outfile,header=header)

    cuttype='typeb'
    if dotypebcut:
        acceptlist=subtractAcceptlist(\
            patch=patch,cutsdir=cutsdir,acceptlist=acceptlist,cuttype=cuttype)
        dummy=acceptlist.printSummary()
    else: print '**%s cuts not applied' % cuttype
    if quit==cuttype:
        finalize(acceptlist,cutsdir,outfile,header=header)

    cuttype='weather'
    if doweathercut:
        acceptlist=subtractAcceptlist(\
            patch=patch,cutsdir=cutsdir,acceptlist=acceptlist,cuttype=cuttype)
        dummy=acceptlist.printSummary()
    else: print '**%s cuts not applied' % cuttype
    if quit==cuttype:
        finalize(acceptlist,cutsdir,outfile,header=header)

    cuttype='sunavoid'
    if dosuncut:
        acceptlist=subtractAcceptlist(\
            patch=patch,cutsdir=cutsdir,acceptlist=acceptlist,cuttype=cuttype)
        dummy=acceptlist.printSummary()
    else: print '**%s cuts not applied' % cuttype
    if quit==cuttype:
        finalize(acceptlist,cutsdir,outfile,header=header)

    cuttype='phsw'
    if dophswcut:
        print '**%s cut not implemented' % cuttype
    else:
        print '**%s cut not applied (or implemented)' % cuttype

    cuttype='fftcuts'
    if dofftcuts:
        acceptlist=subtractAcceptlist(\
            patch=patch,cutsdir=cutsdir,acceptlist=acceptlist,cuttype=cuttype)
        dummy=acceptlist.printSummary()
    else: print '**%s cuts not applied' % cuttype
    if quit==cuttype:
        finalize(acceptlist,cutsdir,outfile,header=header)

    cuttype='glitch'
    if doglitchcut:
        acceptlist=subtractAcceptlist(\
            patch=patch,cutsdir=cutsdir,acceptlist=acceptlist,cuttype=cuttype)
        dummy=acceptlist.printSummary()
    else: print '**%s cuts not applied' % cuttype
    if quit==cuttype:
        finalize(acceptlist,cutsdir,outfile,header=header)


    finalize(acceptlist,cutsdir,outfile,header=header)
    sys.exit(0)
