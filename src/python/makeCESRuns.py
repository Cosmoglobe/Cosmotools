#!/usr/bin/env python
import optparse
import os
import shutil
import sys
EXE_DEFAULT="/usit/titan/u1/joez/quiet/oslo/src/f90/descart/descart_oslo"
TEMPLATEDIR_DEFAULT="/usit/titan/u1/joez/sample"
MODLIST_DEFAULT=TEMPLATEDIR_DEFAULT+os.path.sep+"module_list.txt"
NODES_DEFAULT=4

TASKS=16

PATCH_DEFAULT='patch_6a'
usage="%prog [options] accept_file run_name [additional parameters in the form NAME=val]\nRun with option -h for help"
parser=optparse.OptionParser(usage)

parser.add_option("-p","--patch",dest='patch',type='string',default=PATCH_DEFAULT,help='The name of the patch to map (%s)'%PATCH_DEFAULT)
parser.add_option("-n","--nodes",dest='nodes',type='int',default=NODES_DEFAULT,help='The number of nodes for the batch file (%d)' % NODES_DEFAULT)
parser.add_option("-x","--exe",dest='executable',type='string',default=EXE_DEFAULT,help='The mapmaker executable for the batch file (%s)' % EXE_DEFAULT)
parser.add_option("-m","--modules",dest='modules',type='string',default=MODLIST_DEFAULT,help='The list of modules to be used (%s)' %MODLIST_DEFAULT)
parser.add_option("-t","--template",dest='templatedir',type='string',default=TEMPLATEDIR_DEFAULT,help='The directory containing templates for the parameter files (%s)'%TEMPLATEDIR_DEFAULT)
#parser.add_option("-s","--output",dest='output',type='string',default="",help='The file to contain output from the batch job (run_name/run_name_output.txt)')
parser.add_option("-d","--dir",dest='output_dir',type='string',default="",help='Existing directory to contain run outputs pixels, map and covariance (run_name)')
parser.add_option("-L","--launch",dest='launch',default=False,action='store_true',help='Launch the job immediately after constructing it (False)')

options,args=parser.parse_args()

try:
	accept_src=os.path.abspath(args[0])
	name=args[1]
except:
	parser.error("Specify the accept file to use and a name for the run")
	
patch=options.patch
nodes=options.nodes
module_list_src=os.path.abspath(options.modules)
template_dir=os.path.abspath(options.templatedir)
exe=os.path.abspath(options.executable)

	
	

accept=name+'_accept.txt'
todproc=name+'_todproc.txt'
common=name+"_common.txt"
module_list=name+"_module_list.txt"


batch=open(template_dir+os.path.sep+'batch.CES.template').read()
start_dir=os.getcwd()

os.mkdir(name)
os.chdir(name)
os.mkdir("batch")
os.mkdir("params")
os.mkdir("accepts")
os.mkdir("results")
import quiet

accept=quiet.AcceptedScanList(accept_src)
runlist=quiet.L2Runlist("/data4/quiet/runlist_l2.txt")
runlist_segments=runlist.segments()
		

full_scans=accept.scans
accept.removeEmptyScans()
accepts={}
for (ID,scan) in full_scans.iteritems():
	accept.scans={ID:scan}
	filename=os.path.abspath("accepts/accept_%d_%d.txt" % ID)
	accepts[ID]=filename
	accept.toFile(filename)
accept.scans=full_scans
os.chdir(template_dir)

params=quiet.ParameterData("common.txt")
os.chdir(start_dir)
os.chdir(name)

params["PARFILE_TOD_PROC"]=todproc
params["MODULE_LIST"]=module_list
err=False
for arg in args[2:]:
	try:
		argname,val=arg.split('=')
		argname=argname.upper()
		if argname not in params.paramVals:
			sys.stderr.write("WARNING: Ignoring unkown parameter %s\n" % argname)
			err=True
		else:
			params[argname]=quiet.ParameterData.txt2param(val)
	except ValueError:
		sys.stderr.write("WARNING: Unable to parse argument parameter %s\n" % arg)
		err=True

batches=[]

for i in xrange(TASKS):
	batch_name = "batch/batch_%d.txt" % i
	output= "batch/output_%d.txt"%i
	open(batch_name,"w").write(batch % (output,os.getcwd()))
	batches.append(batch_name)
	
	outname=os.path.abspath("results") +os.path.sep
for i,ID in enumerate(full_scans):
	try:
		params["TARGET_NAME"]=runlist_segments[ID].targetName
		params["TARGET_TYPE"]=runlist_segments[ID].targetType
	except KeyError:
		print "Accepted segment (%d,%d) not found in runlist"%ID
		continue
	params["ACCEPTED_SCANS"]=accepts[ID]
	params["OUTPUT_MAP_FILENAME"]=outname + ("map_%d_%d.fits"%ID)
	params["COV_FILENAME"]=outname + ("cov_%d_%d.txt"%ID)
	params["PIXELS_FILENAME"]=outname + ("pix_%d_%d.txt"%ID)
	tod2map='params/tod2map_%d_%d.txt' % ID
	params["PARFILE_TOD2MAP"]=tod2map	
	common="params/common_%d_%d.txt" % ID
	params.toFiles(common)
	output=outname + ("output_%d_%d.txt"%ID)	
	open(batches[i%TASKS],"a").write("mpirun %s %s | tee %s\n" % (exe,common,output))
	


shutil.copy(module_list_src,module_list)



if options.launch:
	if not err:
		print "Launching job"
		for batch_name in batches:
			os.system("sbatch %s" % batch_name)
	else:
		sys.stderr.write("Cancelling immediate launch: errors noted")

os.chdir(start_dir)
