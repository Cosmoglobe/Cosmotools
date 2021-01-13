import _quietPyExt
#import _getalpha
import numpy

QUIET_SAMPLING_FREQUENCY=50.0  #Hz

def get_alpha(*args):
	"""
	Pass in either a pointing structure (i.e. module.pointing) or three args:
	phi,theta,psi as arrays.
	"""
	if len(args)==1:
		pointing=args[0]
		phi=pointing.phi
		theta=pointing.theta
		psi=pointing.psi
	elif len(args)==3:
		phi,theta,psi=args
	
	assert(len(phi)==len(theta)==len(psi),"Arrays different lengths.")
		
	if not phi.dtype==numpy.float64:
		phi=phi.astype(numpy.float64)
	if not psi.dtype==numpy.float64:
		psi=psi.astype(numpy.float64)
	if not theta.dtype==numpy.float64:
		theta=theta.astype(numpy.float64)
	return _quietPyExt.get_alpha_arr(phi,theta,psi)


def accumulate_slice(target,index,source):
	if isinstance(source,numpy.ndarray):
		_quietPyExt.accumulate_slice1(target,index,source)
	else:
		_quietPyExt.accumulate_slice2(target,index,source)

def get_scan_freq(module,cut=50):
	"""
	Compute the scan freqency of the given module.
	"""
	az=_get_azimuth(module)
	cut=int(cut)
	n=len(az)//cut*cut  #Cut to a multiple of length 50
	az=az[:n]
	faz=numpy.fft.fft(az)
	saz=abs(faz)**2
	n=len(faz)
	nmaz=saz[2:n/2].argmax()+2
	df=QUIET_SAMPLING_FREQUENCY/n
	return nmaz*df

def _get_azimuth(*args):
	"""
	Pass in either a pointing structure (i.e. module.pointing) or three args:
	phi,theta,MJD as arrays.
	"""
	if len(args)==1:
		pointing=args[0].pointing
		phi=pointing.phi
		theta=pointing.theta
		time=args[0].time
	elif len(args)==3:
		phi,theta,time=args
	
	assert(len(phi)==len(theta)==len(time),"Arrays different lengths.")
	if not phi.dtype==numpy.float64:
		phi=phi.astype(numpy.float64)
	if not time.dtype==numpy.float64:
		time=time.astype(numpy.float64)
	if not theta.dtype==numpy.float64:
		theta=theta.astype(numpy.float64)
	az= _quietPyExt.get_azimuth_arr(phi,theta,time)
	return az
