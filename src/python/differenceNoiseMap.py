#!/usr/bin/python
import sys
import numpy as np

UNSEEN=-1.6375e30
def nside2npix(nside):
	return 12*nside*nside


def read_map(filename):
	f=open(filename)
	line=f.readline()
	assert(line.startswith('nside'),"Map format not correct - first line should specify nside")
	nside = int(line.split('=')[1])
	line=f.readline()
	assert(line.startswith('npix'),"Map format not correct - first line should specify npix")
	npix = int(line.split('=')[1])
	f.seek(f.tell())  #Not sure if this is needed; it flushes the read-ahead buffer so we can mix the for line in f syntax and the readline syntax.
	pixels=[]
	values=[]
	for line in f:
		pixel,value = line.split()
		pixels.append(float(pixel))
		values.append(float(value))
	pixels=np.array(pixels,dtype=int)
	values=np.array(values,dtype=float)
	assert(len(pixels)==npix,"Map size is not the same as npix")
	return nside,pixels,values


##############################
#  Main program starts here  #
##############################

if __name__ == '__main__':

	if len(sys.argv)!=3:
           sys.stderr.write("Syntax: differenceNoiseMap map1 map2 > map3\n")
           sys.exit(1)

  	filename1 = sys.argv[1]
	filename2 = sys.argv[2]

        nside1,pix1,rms1=read_map(filename1)
	nside2,pix2,rms2=read_map(filename2)

	assert(nside1==nside2,"Maps have different nsides")
	npix=nside2npix(nside1)

	map1=np.repeat(UNSEEN,npix)
	map2=np.repeat(UNSEEN,npix)

	map1[pix1]=rms1
	map2[pix2]=rms2

	good_pixels = np.where(np.logical_and(map1!=UNSEEN,map2!=UNSEEN))[0]

	output_map = np.sqrt(map1**2+map2**2)
	print "nside=%d"%nside1
	print "npix=%d"%len(good_pixels)
	for p in good_pixels:
        	print "\t%d\t%le" % (p,output_map[p])


############################
#  Main program ends here  #
############################

