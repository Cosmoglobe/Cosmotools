#!/usr/bin/python

"""
Take an acceptlist file and apply the specified jackknife
You need to specify the patch, but only for historical reasons
Usually run as ./jackknifeAccept.py -p2a -iacceptlist.txt -jall
"""
#-------------------------------------------------------------------

if __name__ == '__main__':

    import sys,os,copy
    import quiet,secateur
    from optparse import OptionParser

    print
    # Set up options
    usage = "./jackknifeAccept.py -h' \n"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="make lots of noise")

    parser.add_option("-p", "--patch", dest="patch",\
           help="specify QUIET CMB patch", type="string",\
           metavar='(2a|4a|6a|7b)')

    parser.add_option("-j", "--jackknife", dest="split",\
           help="specify jackknife", type="string",\
           metavar='(all)')

    parser.add_option("-i", "--infile", dest="acc",\
           help="specify input acceptlist", type="string",\
           metavar='accept1.txt')

    parser.add_option("-a", "--outfile1", dest="a",\
           help="specify second acceptlist", type="string",\
           metavar='output.a')

    parser.add_option("-b", "--outfile2", dest="b",\
           help="specify outfile", type="string",\
           metavar='output.b')

    # Assert defaults
    parser.set_defaults(split=None,\
        patch=None,acc=None,a='output.a',b='output.b',verbose=False)
    # Parse command line
    (options, args) = parser.parse_args()

    patch=options.patch
    verbose=options.verbose
    acc=options.acc
    a=options.a
    b=options.b
    split=options.split

    assert(split is not None), 'split is required!'
    assert(acc is not None), 'acceptlist required!'

    # Set up the patch
    assert(patch is not None), 'patch is required!'
    patch = 'patch_%s' % patch
    print 'patch is %s' % patch[-2:]
    assert(patch in secateur.CMBPatches), '%s not available!' % patch

    if split is not 'all':
        assert(a is not None), '-a required!'
        assert(b is not None), '-b required!'

    # Fetch the acceptlist
    acceptlist=quiet.AcceptedScanList(acc)

    if split=='all':
        flavdict={'module': ['radius', 'mab1','mab2','mab3','alternate', \
                             'horizontal','vertical1','vertical2', \
                             'diodes','frequency','typea','leakage',\
                             'modleak'], \
                  'time': ['interleaved', 'season', 'tandem', 'battenberg'], \
                  'point': ['elevation','elsidelobe','deck90','deck45','azimuth'], \
                  'weather': ['pwv','humidity'], \
                  'gain': ['reg'] }
#                  'frequency': ['fknee','chisqsspos','chisqssneg', \
#                                'chisqhighpos','chisqhighneg'] }
    else:
        print 'not implemented'
        sys.exit(0)

    pairs={}
    backup=copy.deepcopy(acceptlist)
    pairs=secateur.SplitAcceptlist(acceptlist=acceptlist,\
                   patch=patch,flavourdict=flavdict,verbose=verbose)
    assert(backup==acceptlist), '*acceptlist has been modified by jackknifing!'

    n=1
    report=''
    rf='accepted_scans.txt'
    reportf=open(rf,'w')
    line='NUM_SETS = %i\n' % (2*len(pairs.keys()))
    reportf.write(line)
    print 'infile was %s' % acc
    for jk in pairs.keys():
        a,b=pairs[jk]
        fa=a.title
        fb=b.title
        if os.path.exists(fa):
            os.rename(fa,'%s.old'%fa)
            print '  (renamed %s -> %s.old)' % (fa,fa)
        if os.path.exists(fb):
            os.rename(fb,'%s.old'%fb)
            print '  (renamed %s -> %s.old)' % (fb,fb)
        print '%s written to %s, %s' % (jk,fa,fb)
        line= "ACCEPTED_SCANS%i       = '%s'\n" %(n,fa)
        reportf.write(line)
        report+=line
        n+=1
        line= "ACCEPTED_SCANS%i       = '%s'\n" %(n,fb)
        reportf.write(line)
        report+=line
        n+=1
        a.toFile(fa)
        b.toFile(fb)

    reportf.close()

    print '%s can be included in your par file' % rf
    print

    #print 'jackknives written to %s' % outfile
    #print '...and you might want this:'
    #print report
    #print
    
    sys.exit(0)
