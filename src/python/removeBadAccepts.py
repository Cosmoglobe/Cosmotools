#!/usr/bin/python
import sys

def usage():
	sys.stderr.write("Syntax: removeBadAccepts.py runlist acceptlist new_acceptlist\n")
	sys.exit(1)

if not len(sys.argv)==4:
	usage()
import quiet

runlist = quiet.L2Runlist(sys.argv[1])
accept = quiet.AcceptedScanList(sys.argv[2])

okaySegments = runlist.segments()
accept.removeScansNotIn(okaySegments)
accept.toFile(sys.argv[3])
