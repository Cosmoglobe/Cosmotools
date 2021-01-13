#!/usr/bin/python
import sys
import os

if len(sys.argv)!=4:
	sys.stderr.write("\nUsage: getSigma.py targetName L2Dir runlist\n")
	sys.stderr.write("e.g. getSigma.py moon /data4/quiet/level2 runlist_l2.txt\n")
	sys.stderr.write("\nFor all runs in the runlist the script finds the corresponding files and prints the standard deviation of all diodes\n\n")
	
	sys.exit(1)

import quiet
	
targetName,l2dir,runlistfile=sys.argv[1:4]

runlist=quiet.L2Runlist(runlistfile)

target=runlist.targets[targetName]
scans = target.scans.values()

scanIDs = target.scans.keys()
scanIDs.sort()

diode_names = ["Q1","U1","U2","Q2"]

for scanID in scanIDs:
	run = target.scans[scanID]
	for segment in run.segments:
		filename = l2dir + os.path.sep + segment.filename()
		l2data=quiet.l2_read(filename)
		for m,module in enumerate(l2data):
			for d,diode_timestream in enumerate(module.tod):
				sigma = diode_timestream.std()
				print "%d\t%d\t%d\t%s\t%e" % (run.id,segment.id,m,diode_names[d],sigma)
