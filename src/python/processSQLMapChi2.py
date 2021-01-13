#!/usr/bin/python

"""
Do the following, but for a file containing Yuji's map chisq values, i.e.
quiet1:'/usit/titan/u1/jonathaz/sidelobes/mapchi2.txt'
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
    from math import cos

    if len(sys.argv) != 3:
        print 'Usage: ./processSQL.py infile outfile'
        sys.exit(1)

    inf=sys.argv[1]
    outf=sys.argv[2]

    print '%s ---> %s' % (inf,outf)

    infile=open(inf,'r')

    if os.path.exists(outf):
        os.rename(outf,outf+'.old')
        print '(renamed %s -> %s.old)' % (outf,outf)
    outfile=open(outf,'w')

    # Fetch the runlists and set up the patch dictionary
    rundict,acceptdict=secateur.GenerateCMBRunlists(\
        patch=secateur.AllPatches[:])
    patchdict=secateur.BuildPatchDict(acceptdict)

    #header='# run seg mod di patch mapchi2 mjd az el dk thetasun phisun\n'
    header='# run seg mod di patch mapchi2 el dk' +\
            ' mjd1 az1 azsun1 elsun1 thetasun1 phisun1' +\
            ' mjd az azsun elsun thetasun phisun' +\
            ' mjd2 az2 azsun2 elsun2 thetasun2 phisun2\n'
    outfile.write(header)

    sqllinef='| %i | %i | %i | %s | %f |'
    flushfreq=10000
    i=1
    daz=8.0 #AZ swing in deg
    for line in infile:
        # Grab the header (didn't quite work in this case!)
        if i==1 and line.startswith('mysql'):
            outfile.write('# '+line)
        if line.startswith('| run_id'): continue # Trap the header
        try:
            scan,seg,module,diname,mapchisq=scanf.sscanf(line,sqllinef)
        except scanf.IncompleteCaptureError:
            continue
        [(oscan,oseg)]=ces_translate([(scan,seg)],reverse=False)
        if oscan != -1:
            diode=secateur.DiodesLookup[diname]
            patch=patchdict[(oscan,oseg)]
            runlist=rundict[patch]
            (az,el,dk,mjd1,mjd2)=secateur.FetchPointing(\
                oscan,oseg,runlist,patch=patch)
            mjd=(mjd1+mjd2)/2.0
            dazel=daz/cos(el)
            azs={mjd1:az-dazel,mjd:az,mjd2:az+dazel}
            suntheta={}; sunphi={}; sunaz={}; sunel={}
            # EL, DK are constant throughout a CES
            for date in [mjd1,mjd,mjd2]:
                suntheta[date],sunphi[date],sunaz[date],sunel[date]=\
                                              secateur.CalculateSunCoords(\
                    azs[date],el,dk,date,verbose=False,returnsunpos=True)
            patch=str(0.0)
            outline=\
              '%i %i %i %i %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n' \
              % (oscan,oseg,module,diode,patch,mapchisq,el,dk,\
                 mjd1,azs[mjd1],\
                 sunaz[mjd1],sunel[mjd1],suntheta[mjd1],sunphi[mjd1],\
                 mjd,azs[mjd], \
                 sunaz[mjd],sunel[mjd],suntheta[mjd],sunphi[mjd],\
                 mjd2,azs[mjd2],\
                 sunaz[mjd2],sunel[mjd2],suntheta[mjd2],sunphi[mjd2])
            outfile.write(outline)
        if (i % flushfreq == 0):
            outfile.flush()
            print 'Reached line %i' % i
        i+=1

    infile.close()
    outfile.close()

    print 'File written to %s' % outf

    sys.exit(0)
