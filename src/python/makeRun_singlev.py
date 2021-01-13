#!/usr/bin/env python
import optparse
import os
import shutil
import sys
EXE_DEFAULT="/usit/titan/u1/joez/quiet/oslo/src/f90/descart/descart_oslo"
TEMPLATEDIR_DEFAULT="/usit/titan/u1/joez/sample"
MODLIST_DEFAULT=TEMPLATEDIR_DEFAULT+os.path.sep+"module_list.txt"
NODES_DEFAULT=4
PATCH_DEFAULT='patch_6a'
usage="%prog [options] accept_file run_name [additional parameters in the form NAME=val]\nRun with option -h for help"
parser=optparse.OptionParser(usage)

parser.add_option("-p","--patch",dest='patch',type='string',default=PATCH_DEFAULT,help='The name of the patch to map (%s)'%PATCH_DEFAULT)
parser.add_option("-n","--nodes",dest='nodes',type='int',default=NODES_DEFAULT,help='The number of nodes for the batch file (%d)' % NODES_DEFAULT)
parser.add_option("-x","--exe",dest='executable',type='string',default=EXE_DEFAULT,help='The mapmaker executable for the batch file (%s)' % EXE_DEFAULT)
parser.add_option("-m","--modules",dest='modules',type='string',default=MODLIST_DEFAULT,help='The list of modules to be used (%s)' %MODLIST_DEFAULT)
parser.add_option("-t","--template",dest='templatedir',type='string',default=TEMPLATEDIR_DEFAULT,help='The directory containing templates for the parameter files (%s)'%TEMPLATEDIR_DEFAULT)
parser.add_option("-s","--output",dest='output',type='string',default="",help='The file to contain output from the batch job (run_name/run_name_output.txt)')
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

if options.output_dir:
	outname=os.path.abspath(options.output_dir)+os.path.sep+name
else:
	outname=name
	
if options.output:
	output=os.path.abspath(options.output)
else:
	output=name+"_output.txt"
	
mapname=outname+'_map.fits'
pixels=outname+"_pixels.txt"
cov=outname+"_cov.dat"

accept=name+'_accept.txt'
tod2map=name+'_tod2map.txt'
todproc=name+'_todproc.txt'
common=name+"_common.txt"
module_list=name+"_module_list.txt"


batch=open(template_dir+os.path.sep+'batch.template').read()
start_dir=os.getcwd()

os.mkdir(name)

os.chdir(template_dir)

import quiet
params=quiet.ParameterData("common.txt")
os.chdir(start_dir)
os.chdir(name)

params["TARGET_NAME"]=patch
params["ACCEPTED_SCANS"]=accept
params["COV_FILENAME"]=cov
params["OUTPUT_MAP_FILENAME"]=mapname
params["PARFILE_TOD2MAP"]=tod2map
params["PARFILE_TOD_PROC"]=todproc
params["PIXELS_FILENAME"]=pixels
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

params.toFiles(common)

shutil.copy(module_list_src,module_list)
shutil.copy(accept_src,accept)
batch_name = name+"_batch.txt"
open(batch_name,"w").write(batch % (nodes,output,os.getcwd(),exe,common))

if options.launch:
	if not err:
		print "Launching job"
		os.system("sbatch %s" % batch_name)
	else:
		sys.stderr.write("Cancelling immediate launch: errors noted")

os.chdir(start_dir)
