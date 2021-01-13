#!/bin/env python
import sys
import numpy
import pyfits
import healpy
import pylab
from quiet import patch_coords, read_pixfile



def plot_partial(nside,pix,vals,filename,patch,min=None,max=None):
	full_map=numpy.repeat(healpy.UNSEEN,healpy.nside2npix(nside))
	full_map[pix]=vals
	rot=patch_coords[patch]
	healpy.gnomview(full_map,rot=rot,reso=12.0,max=max,min=min)
	pylab.savefig(filename)


if __name__=="__main__":
	import optparse
	parser=optparse.OptionParser(usage="Convert %prog infile outfile\nConvert a pixel file to an image.")
	parser.add_option("-p","--patch",dest="patch",type="string",default="patch_6a",help="The name of the patch [patch_6a]")
	parser.add_option("-m","--min",dest="min",type="float",default=None,help="The minimum plot range value [auto]")
	parser.add_option("-M","--max",dest="max",type="float",default=None,help="The maximum plot range value [auto]")

	opt,args=parser.parse_args()
	if len(args)!=2:
		parser.print_help()
		parser.exit()
	infile,outfile=args
	

	nside,pix,vals=read_pixfile(infile)
#	save_fits(nside,pix,vals,outfile)
	plot_partial(nside,pix,vals,outfile,opt.patch,min=opt.min,max=opt.max)
	
