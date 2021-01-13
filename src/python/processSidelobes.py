#!/usr/bin/python

"""
Take IDL-preprocessed beam-mapping data
Convert to theta,phi coordinates
Usage:
./plotSidelobes.py > ~/sidelobes/sidelobes_proc.dat
"""

if __name__ == '__main__':

    import os,sys
    import scanf
    import secateur

    #if len(sys.argv) != 3:
    #    print 'Usage: ./processSQL.py infile outfile'
    #    sys.exit(1)

    dk=float(sys.argv[1])
    inf=sys.argv[2]

    #print '%s ---> %s' % (inf,outf)

    #inf='/usit/titan/u1/jonathaz/sidelobes/sidelobes.dat'
    infile=open(inf,'r')

    #if os.path.exists(outf):
    #    os.rename(outf,outf+'.old')
    #outfile=open(outf,'w')

    linef=' %i %f %f %f %f %i %i %i %i '
    #dk=90.0
    azsource=20.0
    elsource=-5.0
    print '# index az el dk theta phi ddl ddl_db azsource elsource'
    for line in infile:
        index,az,el,ddl,ddl_db,dd0,dd1,dd2,dd3=scanf.sscanf(line,linef)
        theta,phi=secateur.Horiz2BoresightSystem(az,el,dk,azsource,elsource,\
                                       degrees=True)
        print index,az,el,dk,theta,phi,ddl,ddl_db,azsource,elsource

        #outline='%i %i %i %i\n' % (oscan,oseg,module,diode)
        #outfile.write(outline)

    infile.close()
    #outfile.close()

    #print 'File written to %s' % outf

    sys.exit(0)
