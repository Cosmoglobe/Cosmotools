#!/usr/bin/env python
import matplotlib
matplotlib.use('agg')
import struct
import healpy as hp
import quiet 
import numpy as np
import pylab
from numpy import zeros,arange,cos,sin,sqrt,where,ones,identity,array,any,pi,radians,diag,int32,float64
import os
from os import path

import sys
from quiet import NDIODE_MODULE
sys.path.append("/usit/titan/u1/joez/quiet/oslo/src/python/extensions")
from quietPyExt import accumulate_slice,get_alpha
import pprint
import optparse
QUIET_SAMPLING_FREQUENCY=50.

def parseLevel2Filename(filename):
	l2dir,filename=os.path.split(filename)
	filename,suffix=os.path.splitext(filename)
	parts=filename.split("_")
	targetType=parts[0]
	targetName="_".join(parts[1:-2])
	run=int(parts[-2])
	segString=parts[-1]
	seg=int(segString[3:])
	return l2dir,targetType,targetName,run,seg
	
	
class verbosity:
	silent=0
	quiet=1
	noisy=2
	debug=4

if "-x" in sys.argv:
	sys.argv.remove("-x")
	use_mpi=False
else:
	try:
		import mpi4py.MPI as MPI
		use_mpi=True
		print "%d: MPI Enabled" % MPI.COMM_WORLD.rank
	except ImportError:
		use_mpi=False

class MapMaker(object):
	def __init__(self,nside,nest=False,rank=0,verbosity=verbosity.debug):
		self.rank=rank
		self.verbosity=verbosity
		self.nside=nside
		npix=hp.nside2npix(self.nside)
		self.npix=npix
		self.C0=zeros(npix)
		self.C1=zeros(npix)
		self.C2=zeros(npix)
		self.q=zeros(npix)
		self.u=zeros(npix)
		self.hits=zeros(npix)
		self.NDIODE_MODULE=4
		self.nest=nest

	def report(self,threshold,message,*args):
		if self.verbosity>=threshold:
			message=message%tuple(args)
			print "%d: %s" % (self.rank,message)
		
	def add_level2(self,filename,data):
		"""
		Data should be an array of dictionaries of module properties.
		"""
		polarization_modules=range(17)
		modules=quiet.l2_read(filename,modules=polarization_modules)
		self.report(verbosity.quiet,"Adding file %s",filename)
		for m in modules:
			n=int(m.number)
			self.report(verbosity.noisy,"Processing module %d",n)
			self.process_module(m,data[n])
			self.add_module(m)
			
	def process_module(self,m,data):
		#Diode acceptance
		subtract_mean = data.get("subtract_mean",True)
		subtract_linear = data.get("subtract_linear",False)
		m.diode_flag=data.get("diode_flag",ones(NDIODE_MODULE))
		self.report(verbosity.debug,"	Diode flags: %s",m.diode_flag)
		
		#Apply gains and remove TOD means
		gains=data.get("gains",ones(NDIODE_MODULE))
		for i in range(NDIODE_MODULE):
			self.report(verbosity.debug,"    Processing diode %d",i)
			self.report(verbosity.debug,"        Applying gain %e",gains[i])
			tod=m.tod[i]
			tod/=gains[i]
			if subtract_mean: tod -= tod.mean()
			if subtract_linear: self.subtract_linear(tod)
		
		#Compute dpsi from fiducials
		m.dpsi=array([0,pi/4,-pi/4,-pi/2]) + data.get("dpsi",0)
		
		#Scale noise with gain and construct matrix.
		sigma=data.get("sigma",ones(NDIODE_MODULE))/gains
		self.report(verbosity.debug,"        Module noises = %s",sigma)
		m.invCw= diag(1/sigma**2)
		
	def subtract_linear(self,y):
		n=len(y)
		xsum=n*(n+1.0)/2.0
		x2sum=n*(n+1.0)*(2*n+1.0)/6.0
		x=arange(1,n+1)
		xysum=(x*y).sum()
		ysum=y.sum()
		m=(n*xysum - ysum*xsum) / (n*x2sum - xsum**2)
		c=ysum/n-m*xsum/n
		fit=m*x+c
		self.report(verbosity.noisy,"        Subtracting linear fit m=%.3e,c=%.3e",m,c)
		y-=fit
	def add_module(self,m):
		if not np.any(m.diode_flag): return
		pixels=hp.ang2pix(self.nside,m.pointing.theta,m.pointing.phi,self.nest)
		self.hits[pixels]=True
		NDIODE_MODULE=self.NDIODE_MODULE
		
		
		#Compute the rotation angles for each diode.
		alpha = get_alpha(m.pointing)
		c=[None]*NDIODE_MODULE
		s=[None]*NDIODE_MODULE
		for i in range(NDIODE_MODULE):
			if m.diode_flag[i]:
				c[i]=cos(2*(alpha+m.dpsi[i]))
				s[i]=sqrt(1-c[i]**2)
				
		#Accumulate the covariance matrix contribution from this module
		for i in range(NDIODE_MODULE):
			if m.diode_flag[i]:
				for j in range(NDIODE_MODULE):
					if m.diode_flag[j]:
						accumulate_slice(self.C0,pixels,c[i]*m.invCw[i,j]*c[j])
						accumulate_slice(self.C1,pixels,c[i]*m.invCw[i,j]*s[j])
						accumulate_slice(self.C2,pixels,s[i]*m.invCw[i,j]*s[j])

		#Weight by the diode noise covariance
		invCw_tod = self.inverseNoiseWeight(m.tod,m.diode_flag,m.invCw)
#		for i in range(NDIODE_MODULE):
#			tod=zeros(len(m.tod[i]))
#			for j in range(NDIODE_MODULE):
#				if m.diode_flag[j]:
#					tod+=m.invCw[i,j]*m.tod[j]
#			invCw_tod.append(tod)
			
		#Accumulate the maps from this module
		for i in range(NDIODE_MODULE):
			if m.diode_flag[i]:
				accumulate_slice(self.q,pixels,c[i]*invCw_tod[i])
				accumulate_slice(self.u,pixels,s[i]*invCw_tod[i])

	def inverseNoiseWeight(self,tods,flags,invCw):
		invCw_tod = []
		for i in range(NDIODE_MODULE):
			if flags[i]:
				tod=zeros(len(tods[i]))
				for j in range(NDIODE_MODULE):
					if flags[j]:
						tod+=invCw[i,j]*tods[j]
				invCw_tod.append(tod)
			else:
				invCw_tod.append(0)
		return invCw_tod
		

	def make(self):
		hit=where(self.hits)[0]

		C0=self.C0[hit]
		C1=self.C1[hit]
		C2=self.C2[hit]
		q=self.q[hit]
		u=self.u[hit]

		#Invert the covariance matrices
		DET=C0*C2-C1**2	
		C0,C1,C2 = (C2/DET, -C1/DET, C0/DET)
		#Generate the returned map.
		q,u  = ( C0*q+C1*u ), ( C1*q+C2*u)
		return hit,q,u

class OffsetsData(object):
	def __init__(self,filename):
		f=quiet.FortranFile(filename)
		self.lc=f.readArray(int32)
		self.na=f.readArray(int32)
		self.values=f.readArray(float64)
		f.close()


class DummyContainer(object):
	pass

class OffsetMapMaker(MapMaker):
	def process_module(self,m,data):
		MapMaker.process_module(self,m,data)
		self.report(verbosity.noisy,"Removing offsets")
		offsets_files=data["offsets_files"]

		offsets_data=[]
		offsets_values=[]
		
#		print m.invCw
		
		for i in xrange(self.NDIODE_MODULE):
			if m.diode_flag[i]:
				filename = offsets_files[i]
				off=OffsetsData(filename)
				offsets_data.append(off)
				offsets_values.append(off.values)
			else:
				offsets_data.append(0)
				offsets_values.append(0)
		weighted_values = self.inverseNoiseWeight(offsets_values,m.diode_flag,m.invCw)

		for i in xrange(self.NDIODE_MODULE):
			if not m.diode_flag[i]:
				continue
			offsets=offsets_data[i]
			na=offsets.na
			lc=offsets.lc
			self.report(verbosity.debug,"        Diode %d Offset std dev =  %f",i,weighted_values[i].std())
			tod=m.tod[i]
			for a in xrange(offsets.na):
				tod[a*lc:(a+1)*lc]-=weighted_values[i][a]
			#Cut down the TOD so we do not have unmodified part baselines at the end.
			tod=tod[0:(offsets.na+1)*lc]
			m.tod[i]=tod

class OffsetTester(OffsetMapMaker):  #For convenience we inherit from OffsetMapMaker, since we want to use its filters
	def __init__(self,output_dir,rank=0,verbosity=verbosity.debug):
		self.rank=rank
		self.verbosity=verbosity
		self.NDIODE_MODULE=4
		self.output_dir=output_dir
	def make(self):
		raise RuntimeError("The OffsetTester cannot make maps")
	def add_level2(self,filename,data):
		"""
		Data should be an array of dictionaries of module properties.
		"""
		polarization_modules=range(17)
		l2dir,targetType,targetName,run,seg = parseLevel2Filename(filename)
		modules=quiet.l2_read(filename,modules=polarization_modules)
		self.report(verbosity.quiet,"Adding file %s",filename)
		diode_names=["Q1","U1","U2","Q2"]
		for m in modules:
			n=int(m.number)
			output_file=self.output_dir+os.path.sep+"%d_%d_%d.png" % (run,seg,n)
			self.report(verbosity.noisy,"Processing module %d",n)
			raw_data=[]
			for tod in m.tod:
				length=len(tod)/50*50
				raw_data.append(abs(np.fft.fft(tod[:length]))**2)
			self.process_module(m,data[n])
			TODs=m.tod
			for i,(name,tod) in enumerate(zip(diode_names,TODs)):
				power=abs(np.fft.fft(tod))**2
				N=len(power)
				pylab.subplot(2,2,i+1)
				power_cut=power[1:N/2+1]
				f=(arange(len(power_cut))+1)*QUIET_SAMPLING_FREQUENCY/N
				pylab.loglog(f,power_cut)
				pylab.loglog(f,raw_data[i][1:N/2+1])
			pylab.savefig(output_file)
			pylab.close()
			

MAGIC_NUMBER="BINARYMAP"

def save_pix(filename,nside,P,X):
	npix=len(X)
	f=open(filename,"w")
	f.write("nside=%d\n"%nside)
	f.write("npix=%d\n"%npix)
	for (p,x) in zip(P,X):
		f.write("\t%d\t%e\n"%(p,x))
	f.close()

def save_combined_ascii(fname,nside,p,q,u,qq,qu,uu):
	f=open(fname,"w")
	f.write("%d  %d\n" % (nside,len(p)) ) 
	for data in zip(p,q,u,qq,qu,uu):
		f.write("%d  %e  %e  %e  %e  %e\n"%data)
	f.close()
	
def save_combined_binary(fname,nside,p,q,u,qq,qu,uu):
	f=open(fname,"w")
	f.write(MAGIC_NUMBER)
	npix=len(p)
	f.write(struct.pack('ii',nside,npix))
	np.save(f,p)
	np.save(f,q)
	np.save(f,u)
	np.save(f,qq)
	np.save(f,qu)
	np.save(f,uu)
	f.close()

def read_combined(fname,ascii=True):
	if ascii:
		return read_combined_ascii(fname)
	else:
		return read_combined_binary(fname)

def read_combined_ascii(fname):
	f=open(fname)
	nside,npix=f.next().split()
	nside=int(nside)
	npix=int(npix)
	data=np.loadtxt(f,unpack=True)
	p=data[0].astype(int)
	q,u,qq,qu,uu=data[1:]
	return nside,npix,p,q,u,qq,qu,uu

def read_combined_binary(fname):
	f=open(fname)
	magic=f.read(len(MAGIC_NUMBER))
	if magic!=MAGIC_NUMBER:
		raise IOError("File is not of correct combined format")
	nside,npix=struct.unpack('ii',f.read(8))
	p,q,u,qq,qu,uu=(np.load(f) for i in xrange(6) )
	return nside,npix,p,q,u,qq,qu,uu
	
	
def save_combined(fname,nside,p,q,u,qq,qu,uu,ascii=True):
	if ascii:
		save_combined_ascii(fname,nside,p,q,u,qq,qu,uu)
	else:
		save_combined_binary(fname,nside,p,q,u,qq,qu,uu)



		

def main(common_file,output_dir,naive=False,test=False,restart=False):
	if use_mpi:
		rank=MPI.COMM_WORLD.rank
		size=MPI.COMM_WORLD.size
	else:
		rank=0
		size=1
	
	print "%d: Running from common file: %s"%(rank,common_file)
	parameters=quiet.ParameterData(common_file)

	runlist_file=parameters["SCAN_DATABASE"]
	accept_file=parameters["ACCEPTED_SCANS"]
	calibration_file=parameters["PARAMETER_DATABASE"]
	nside=parameters["NSIDE_OUT"]
	scaling=True
	
	try:
		offset_dir=parameters["OFFSET_DIR"]
	except KeyError:
		offset_dir="saved_offsets"
	print "%d: Loading calibration from %s"%(rank,calibration_file)
	calibration=quiet.Calibration(calibration_file)
	print "%d: Loading runlist from %s"% (rank,runlist_file)
	runlist=quiet.L2Runlist(runlist_file)

	print "%d: Loading accepts from %s"%(rank,accept_file)
	accept=quiet.AcceptedScanList(accept_file)
	def filenamesForOutput(scanID,segID):
		qname = output_dir+("/map.q.%d.%d.pix"%(scanID,segID))
		uname = output_dir+("/map.u.%d.%d.pix"%(scanID,segID))
		return qname,uname
	segments=runlist.segments()
	scans=accept.scans.keys()
	for scan in scans[rank::size]:		
		segment=segments[scan]
		filename=segment.filename("/data4/quiet/level2/")
		calib_data=[]
		scanID=segment.runID
		segID=segment.id
		qname,uname=filenamesForOutput(scanID,segID)
		if restart and os.path.exists(qname) and os.path.exists(uname):
			print "%d Skipping %d,%d" % (rank,scanID,segID)
			continue
		print "%d: Making map from filename:%s"% (rank,filename)
		
		if scaling:
			print "%d: Scaling gains by 10^-6  to convert from mV/K to V/mK" % rank
		
		for i in xrange(17):
			module_data = build_module_info(i,accept,calibration,segment)
			calib_data.append(module_data)
		
		if test:
			mapper=OffsetTester(output_dir,rank=rank)
			mapper.add_level2(filename,calib_data)
		elif naive:
			mapper=MapMaker(nside,rank=rank)
			mapper.add_level2(filename,calib_data)
			p,q,u=mapper.make()
			save_pix(qname,nside,p,q)
			save_pix(uname,nside,p,u)
		else:
			mapper=OffsetMapMaker(nside,rank=rank)
			mapper.add_level2(filename,calib_data)
			p,q,u=mapper.make()
			save_pix(qname,nside,p,q)
			save_pix(uname,nside,p,u)


def build_module_info(module,accept,calibration,segment,scaling=True,offset_dir='saved_offsets',subtract_linear=False):
	
	date=(segment.startDate+segment.endDate)/2
	scanID=segment.runID
	segID=segment.id
	
	gains=[calibration.lookup_date("gain",module,diode,date,default=1.0) for diode in xrange(NDIODE_MODULE)]

	gains=array(gains)
	if scaling: #Convert from mV/K to V/mK
		gains*=1.0e-6

	dpsi=[calibration.lookup_date("diode_angle",module,diode,date,default=0) for diode in xrange(NDIODE_MODULE)]
	dpsi=array([radians(d) for d in dpsi])

	diode_flag=[accept.isIncluded(scanID,segID,module,diode) for diode in xrange(NDIODE_MODULE)]
	diode_flag=array(diode_flag)

	sigma=[calibration.lookup_scan("sigma0",module,diode,scanID,segID) for diode in xrange(NDIODE_MODULE)]
	sigma=array(sigma)
#			invCw=diag(1.0/array(sigma)**2)

	offsets_files = ["%s/savedOffset_%.5d_%.3d_%.3d_%.2d.off"% (offset_dir,scanID,segID,module,diode+1) for diode in xrange(NDIODE_MODULE)]

	module_data={"gains":gains,"dpsi":dpsi,"diode_flag":diode_flag,"sigma":sigma,"offsets_files":offsets_files,"subtract_linear":subtract_linear}  #add gain and other calibration here.
	return module_data
	

	
if __name__ == '__main__':


	parser=optparse.OptionParser(usage="Usage: %prog [options] common_file output_dir")
	parser.add_option("-n",default=False,action='store_true',dest='naive',help='Make naive CES maps [false]')
	parser.add_option("-d",default=True,action='store_true',dest='destripe',help='Make destriped CES maps [true]')
	parser.add_option("-t",default=False,action='store_true',dest='test',help='Make plots of destriped CES timestreams [false]')
	parser.add_option("-r",default=False,action='store_true',dest='restart',help='Restart from earlier run [do not use this for now]')
	parser.add_option("-i","--diodes",default=False,action='store_true',dest='diodes',help='Make diode maps')

	opt,args=parser.parse_args()
	
	if opt.diodes:
		runlist=quiet.L2Runlist("/data4/quiet/runlist_l2.txt")
		calibration=quiet.Calibration("/data4/quiet/joez/calibration/standard_calibration.dat")
		l2dir="/data4/quiet/level2/"
		root=args[0]
		make_diode_maps(runlist,calibration,root,l2dir,nside=512)
		sys.exit(0)
	
	if len(args)!=2:
		parser.error("Two arguments required: common_file and output_dir")
	common_file=args[0]
	output_dir=args[1]
   	main(common_file,output_dir,naive=opt.naive,test=opt.test,restart=opt.restart)
