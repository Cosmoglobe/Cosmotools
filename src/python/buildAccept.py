#!/usr/bin/python
import sys

try:
	infile=sys.argv[1]
	outfile=sys.argv[2]
except:
	sys.stderr.write("\nUsage: buildAccept.py runlist accepted_list\nRunlist should exist already, accepted_list is output\n\n")
	sys.exit(1)

import quiet

runlist = quiet.L2Runlist(infile)
acceptlist = quiet.acceptedListFromRunList(runlist)
acceptlist
print 'NB: Removing diodes 8/2 and 16/[0-3]'
acceptlist.excludeDiode(8,2)
acceptlist.excludeModule(16)
acceptlist.toFile(outfile)
