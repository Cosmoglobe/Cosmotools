#!/usr/bin/python

"""
Translate an ASCII file containing a Chicago SQL query
(see for example /data4/quiet/nulltests/cuts/chicago/)
Renumber CES IDs to Oslo format
Remove any '-1' scans
Number diodes
Output to file
"""

if __name__ == '__main__':

    import os,sys
    import scanf
    import secateur
    from ces_translate import ces_translate

    if len(sys.argv) != 3:
        print 'Usage: ./processSQL.py infile outfile'
        sys.exit(1)

    inf=sys.argv[1]
    outf=sys.argv[2]

    print '%s ---> %s' % (inf,outf)

    infile=open(inf,'r')

    if os.path.exists(outf):
        os.rename(outf,outf+'.old')
    outfile=open(outf,'w')

    sqllinef='| %i | %i | %i | %s |'
    flushfreq=10000
    i=1
    for line in infile:
        # Grab the header
        if i==1 and line.startswith('mysql'):
            outfile.write('# '+line)
        if line.startswith('| run_id'): continue # Trap the header
        try:
            scan,seg,module,diname=scanf.sscanf(line,sqllinef)
        except scanf.IncompleteCaptureError:
            continue
        [(oscan,oseg)]=ces_translate([(scan,seg)],reverse=False)
        if oscan != -1:
            diode=secateur.DiodesLookup[diname]
            outline='%i %i %i %i\n' % (oscan,oseg,module,diode)
            outfile.write(outline)
        if (i % flushfreq == 0):
            outfile.flush()
            print 'Reached line %i' % i
        i+=1

    infile.close()
    outfile.close()

    print 'File written to %s' % outf

    sys.exit(0)
