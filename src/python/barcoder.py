#!/usr/bin/python

"""
Wrapper to generate a difference barcode from two acceptlists,
for a given patch
"""
#-------------------------------------------------------------------

if __name__ == '__main__':

    import sys,os
    import quiet,secateur
    from optparse import OptionParser

    print
    # Set up options
    usage = "%s -h\n" % sys.argv[0]
    parser = OptionParser(usage=usage)
    #parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
    #                  help="make lots of noise")

    parser.add_option("-p", "--patch", dest="patch",\
           help="specify QUIET CMB patch", type="string",\
           metavar='(2a|4a|6a|7b)')

    parser.add_option("-a", "--infile1", dest="a1",\
           help="specify first acceptlist", type="string",\
           metavar='in1.txt')

    parser.add_option("-b", "--infile2", dest="a2",\
           help="specify second acceptlist", type="string",\
           metavar='in2.txt')

    parser.add_option("-o", "--out", dest="outfile",\
           help="specify outfile", type="string",\
           metavar='out.png')

    # Assert defaults
    parser.set_defaults(\
        patch=None,a1=None,a2=None,outfile='out.png',verbose=False)
    # Parse command line
    (options, args) = parser.parse_args()

    patch=options.patch
    #verbose=options.verbose
    a1=options.a1
    a2=options.a2
    outfile=options.outfile

    # Set up the patch
    assert(patch is not None), 'patch is required!'
    patch = 'patch_%s' % patch
    print 'patch is %s' % patch[-2:]
    assert(patch in secateur.CMBPatches), '%s not available!' % patch

    assert(a1 is not None), '-a required!'
    assert(a2 is not None), '-b required!'

    # Fetch the acceptlists
    acc1=quiet.AcceptedScanList(a1)
    acc2=quiet.AcceptedScanList(a2)

    # Fetch the runlist and acceptlist for this patch
    rundict,acceptdict=secateur.GenerateCMBRunlists()
    acceptlist=acceptdict[patch]

    x=secateur.DifferenceAcceptlistsToMatrix(\
        a1=acc1,a2=acc2,outf=outfile,master=acceptdict[patch])

    print 'infiles were'
    print 'a1 = %s - %s' % (a1,acc1.title)
    print 'a2 = %s - %s' % (a2,acc2.title)
    print
    print 'barcode written to %s' % outfile
    print
    
    sys.exit(0)
