#!/usr/bin/python
import quiet
import sys
import optparse

static_cut_file = "/data4/quiet/calib_data/laura/static_cuts.cut"
default_runlist = "/data4/quiet/runlist_l2.txt"
usage = "usage: %prog " + ( "[options] [cutfiles ...]\n\nWith no arguments this script uses loads the default runlist:\n%s\nand builds an accept list from it and then applies the static cuts in:\n%s\nAdditional arguments are further cutlist files that you want to apply." % (default_runlist,static_cut_file))
parser =  optparse.OptionParser(usage=usage)
parser.add_option("-r", dest="runlist", type="string",default=default_runlist, help="Use specified runlist (default: %s)"%default_runlist)
parser.add_option("-a", dest="acceptlist", type="string",default="", help="Start from the specified accept list.  Default is to build from the runlist")
parser.add_option("-x", action="store_true", dest="ignore_static", default=False, help="Omit the usual static cut in %s" % static_cut_file)
options,cut_files = parser.parse_args()

if not options.ignore_static:
	cut_files.append(static_cut_file)

runlist = quiet.L2Runlist(options.runlist)

if options.acceptlist:
	accept = quiet.AcceptedScanList(options.acceptlist)
else:
	accept = quiet.acceptedListFromRunList(runlist)

	
for filename in cut_files:
	cutlist = quiet.CutList(filename)
	cutlist.pruneAcceptList(runlist,accept)

accept.toFile(sys.stdout)
	
		
