import numpy as np
import quiet
import quietPyExt
import mapping
import healpy
import quiet_filters
from os import path
import os

CUT_FACTOR=1*2*3*5*7

def make_signal_tod(pix,psi,maps,module,dpsi):
	sky_q=maps[0][pix]
	sky_u=maps[1][pix]
	for dp,tod in zip(dpsi,module.tod):
		c = np.cos(2*(psi+dp))
		s = np.sqrt(1-c**2)
		tod[:]= c*sky_q + s*sky_u

def cut_module(module,cut_length):
	module.tod.Q1=module.tod.Q1[:cut_length]
	module.tod.U1=module.tod.U1[:cut_length]
	module.tod.Q2=module.tod.Q2[:cut_length]
	module.tod.U2=module.tod.U2[:cut_length]
	module.pointing.theta=module.pointing.theta[:cut_length]
	module.pointing.phi=module.pointing.phi[:cut_length]
	module.pointing.psi=module.pointing.psi[:cut_length]
	module.time=module.time[:cut_length]
	module.nsamp=cut_length
	


def make_diode_maps(runlist,calibration,root,l2dir,nside=256,use_mpi=False,nest=True,cut_factor=CUT_FACTOR):
	if use_mpi:
		rank=MPI.COMM_WORLD.rank
		size=MPI.COMM_WORLD.size
	else:
		rank=0
		size=1
	accept=quiet.acceptedListFromRunList(runlist)
	segments=runlist.segments().values()
	signal_maps=[]

	for i in xrange(2):
		filename="/xanadu/project/quiet/joez/signal_simulations/map_%d.fits"
		q=healpy.read_map(filename%i,1)
		u=healpy.read_map(filename%i,2)
		signal_maps.append((q,u))
		
	for segment in segments[rank::size]:
		infile = segment.filename(l2dir)
		scan_data = quiet.l2_read(infile)
		print "Loaded %s"%infile
		dirname = path.join(root,segment.targetType,segment.targetName)
		ntod = (scan_data[0].time.size/cut_factor)*cut_factor
		filter=quiet_filters.build_filter(ntod)
		print "Built filter length %d" % ntod
		if not path.exists(dirname):
			print "Making %s" % dirname
			os.makedirs(dirname)		
		for module_number in xrange(16):  #Skip broken module 16
			module_data = scan_data[module_number]
			cut_module(module_data,ntod)
			module_info = mapping.build_module_info(module_number,accept,calibration,segment)
			pix=healpy.ang2pix(nside,module_data.pointing.theta,module_data.pointing.phi,nest)
#			pix=quietPyExt.get_pixel(module_data.pointing)
			psi=quietPyExt.get_alpha(module_data.pointing)
			print "Got module data"
			for ns,signal_map in enumerate(signal_maps):
				make_signal_tod(pix,psi,signal_map,module_data,module_info['dpsi'])
				print "Faked signal"
				quiet_filters.apply_filter(filter,module_data)
				print "Filter complete"
				for diode in xrange(quiet.NDIODE_MODULE):
					mapmaker=mapping.MapMaker(nside)
					if diode==0: mapmaker.process_module(module_data,module_info)
					module_data.diode_flag[:]=False
					module_data.diode_flag[diode]=True
					mapmaker.add_module(module_data)
					p,q,u=mapmaker.make()
					C_qq=mapmaker.C0[p]
					C_qu=mapmaker.C1[p]
					C_uu=mapmaker.C2[p]
					filename = "%d_%.3d_%d_%d.map" % (segment.runID,segment.id,module_number,diode)
					outfile = path.join(dirname,filename)
					print outfile
					mapping.save_combined(outfile,nside,p,q,u,C_qq,C_qu,C_uu,ascii=False)

def main():
	runlist=quiet.L2Runlist("/data4/quiet/runlist_l2.txt")
	print "Loaded runlist" 
	calibration=quiet.Calibration("/data4/quiet/joez/calibration/standard_calibration.dat")
	print "Loaded calibration"
	l2dir="/data4/quiet/level2_old/"
	root='./test1'
	make_diode_maps(runlist,calibration,root,l2dir)
	
if __name__=="__main__":
	main()