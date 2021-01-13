#!/usr/bin/python

import sys
import quiet, secateur

if len(sys.argv) != 4:
      print >> sys.stderr, '%s <cutlist> <flavour of jackknife> <list of dead diodes>' % sys.argv[0]
      sys.exit(2)

print sys.argv[1]
cutces = sys.argv[1]
jack = sys.argv[2]
deadfile = sys.argv[3]

runlist,fulllist=secateur.GenerateCMBRunlists()
remaining={}
cutlist={}
cutlist['patch_6a']=quiet.AcceptedScanList(cutces)

dead=quiet.CutList(deadfile)

quiet.CutList.pruneAcceptList(dead,runlist['patch_6a'],fulllist['patch_6a'])

print 'number of CESes in full list', secateur.acount(fulllist['patch_6a'])
print 'number of CESes in cut list',secateur.acount(cutlist['patch_6a'])

remaining['patch_6a']=secateur.multizap(fulllist['patch_6a'],cutlist['patch_6a'])
print 'number of CESes in runlists',secateur.acount(remaining)

secateur.CutTime(remaining,flavour=jack)

sys.exit(0)
