#!/usr/bin/python

"""
For a set of Type-B dirs, translate them into Oslo numbering

1. Fetch list of directories
2. Translate them using ces_translate
3. Report any tricky cases
4. Build a new list of directories via symlinks
"""
#-------------------------------------------------------------------

if __name__ == '__main__':

    import sys,os,scanf
    from optparse import OptionParser
    from ces_translate import ces_translate

    print
    # Set up options
    usage = "./translate_typeb_dirs.py -i inpath -o outpath' \n"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      help="make lots of noise")

    parser.add_option("-i", "--indir", dest="indir",\
           help="specify input dir", type="string",\
           metavar='indir')

#    parser.add_option("-o", "--outfile", dest="outdir",\
#           help="specify output dir", type="string",\
#           metavar='outdir')

    # Assert defaults
    parser.set_defaults(\
        indir='indir',verbose=False)
    # Parse command line
    (options, args) = parser.parse_args()

    verbose=options.verbose
    indir=options.indir
    #outdir=options.outdir

    assert(indir is not None), 'input directory required!'

    # From www.pixelbeat.org/talks/python/ls.py.html
    infiles=os.listdir(indir)
    # Ignore files starting with '.' using list comprehension
    infiles=[filename for filename in infiles if filename[0] != '.']
    inscans=[tuple([int(s) for s in infile.split('.')]) for infile in infiles]

    inscans.sort()

    #print inscans

    outscans=ces_translate(inscans,reverse=False)

    #print outscans

#    for n,(iscan,iseg) in enumerate(inscans):
#        (oscan,oseg)=outscans[n]
#    #    print '%5s -> %5s' % (str(inscans[n]),str(outscans[n]))
#        if oscan < 0:
#    #os.path.walk(indir,os.getcwd(),'.')

    outdirs=indir.split(os.sep)
    outdirs[-1]='%s_translated' % outdirs[-1]
    outdir=os.sep.join(outdirs)
    print outdir

    reportf=os.path.join(outdir,'nomapping.txt')
    if os.path.exists(reportf): os.rename(reportf,'%s.old'%reportf)
    report=open(os.path.join(outdir,reportf),'w')
    report.write('# These CESs could not be translated (oscan -1)\n')

    for n,(iscan,iseg) in enumerate(inscans):
        (oscan,oseg)=outscans[n]
        if oscan > -1:
            infile=os.path.join(indir,'%s.%s/anomaly.txt'%(iscan,iseg))
            #mkdir cmd
            outd=os.path.join(outdir,'%s.%s'%(oscan,oseg))
            if not os.path.exists(outd): os.makedirs(outd)
            #cmd='ln -s %s %s' % (infile,'%s/anomaly.txt'%outd)
            link='%s/anomaly.txt'%outd
            if os.path.exists(link):
                if os.path.islink(link) and (os.readlink(link) == infile):
                    pass
                else:
                    print 'File incorrectly linked! (%s.%s): %s --x--> %s, -----> %s' % \
                          (oscan,oseg,link,infile,os.readlink(link))
                    print 'Aborting...'
                    sys.exit(0)
            else:
                os.symlink(infile,link)
        else:
            print '%s has no mapping %s' % (str(inscans[n]),str(outscans[n]))
            report.write('%s %s\n' % (inscans[n][0],inscans[n][-1]))

    report.close()

    print 'Look in %s\n' % reportf

    sys.exit(0)
