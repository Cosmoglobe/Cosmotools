import pix2img
import sys

if "-x" in sys.argv:
	sys.argv.remove("-x")
	use_mpi=False
else:
	try:
		import mpi4py.MPI as MPI
		use_mpi=True
		if (MPI.COMM_WORLD.rank==0):
			print "MPI Enabled on %d procs" % MPI.COMM_WORLD.size
	except ImportError:
		use_mpi=False



if __name__=="__main__":
	import optparse
	parser=optparse.OptionParser(usage="Convert %prog infile outfile\nConvert a pixel file to an image.")
	parser.add_option("-p","--patch",dest="patch",type="string",default="patch_6a",help="The name of the patch [patch_6a]")
	parser.add_option("-m","--min",dest="min",type="float",default=None,help="The minimum plot range value [auto]")
	parser.add_option("-M","--max",dest="max",type="float",default=None,help="The maximum plot range value [auto]")

	opt,args=parser.parse_args()
	
	if use_mpi:
		rank=MPI.COMM_WORLD.rank
		size=MPI.COMM_WORLD.size
	else:
		rank=0
		size=1
	files=args[rank::size]
	
	for infile in files:
		outfile=infile+".png"
		print infile
		nside,pix,vals=pix2img.read_pixfile(infile)
		pix2img.plot_partial(nside,pix,vals,outfile,opt.patch,min=opt.min,max=opt.max)
	
