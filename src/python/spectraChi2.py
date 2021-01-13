import quiet
import numpy
import fitPowerSpectrum
import sys



if "-x" in sys.argv:
	sys.argv.remove("-x")
	rank=0
	size=1
else:
	from mpi4py import MPI
	size=MPI.COMM_WORLD.size
	rank=MPI.COMM_WORLD.rank

runlist_file="/data4/quiet/runlist_l2.txt"
calib_file="/data4/quiet/joez/calibration/standard_calibration.dat"
l2_dir="/data4/quiet/level2/"
patch="patch_6a"
QUIET_SAMPLING=50.


runlist=quiet.L2Runlist(runlist_file)
calib=quiet.Calibration(calib_file)
segments=runlist.segments().items()
my_segments=segments[rank::size]
for (run,seg),segment in my_segments:
	if segment.targetName!=patch:
		continue
	filename=segment.filename(l2_dir)
	data=quiet.l2_read(filename)
	for module in data:
		number=module.number[0]
		for diode in xrange(quiet.NDIODE_MODULE):
			fknee=calib.lookup_scan("f_knee",number,diode,run,seg)
			alpha=-calib.lookup_scan("alpha",number,diode,run,seg) #NB Minus sign for alpha convention
			sigma=calib.lookup_scan("sigma0",number,diode,run,seg)
			if numpy.isnan(fknee) or numpy.isnan(alpha) or alpha==0 or fknee==0:
				continue
			chi2=fitPowerSpectrum.powerSpectrumChi2(module.tod[diode],QUIET_SAMPLING,fknee,sigma,alpha)
			print run,seg,number,diode,chi2,fknee,alpha,sigma
			
