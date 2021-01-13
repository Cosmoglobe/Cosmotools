"""
Fit a 1/f**a + white noise power spectrum to a data series 
using a close variant of the method described in 
Dunkley et al astro-ph/0405462
"""

from scipy.optimize import leastsq
from numpy import array , pi, arange, log, where,exp
from numpy.fft import fft,fftfreq
euler_mascheroni = 0.57721566490153286060


def _residuals(p,f,logP):
	model=p[1]+log(1.+(p[0]/f)**p[2])
	return logP-model-euler_mascheroni


def _cut_data(chain,n=50):
	m=len(chain)/n*n
	return chain[:m]

def powerSpectrumChi2(tod,fsamp,fknee,sigma,alpha,fit_length=0,cut_multiple=50):
	"""
	Obtain the chi^2 of the fit to the power spectrum of a TOD.
	The model is
	P(f) = sigma^2 (1+(fknee/f)^alpha )
	Args:
		tod is the timestream data
		fsamp is the sampling frequency
		fknee is the knee frequency (same units as the fsamp)
		sigma is the white noise standard deviation (i.e. close to tod.std())
		alpha is the low frequency noise index (should normally be positive ie 1/f noise has alpha=1)
	Keywords:
		fit_length - the maximum frequency index for the fit.  If not set or set False we use 10 *fknee/fsamp
		cut_multiple - specifies what multiple the chain is to be cut to to improve speed (default 50)
	"""
	chain=_cut_data(tod,cut_multiple)
	f= fftfreq(len(chain),1/fsamp)
	df=f[1]-f[0]
	#Compute the spacing and sizes off the Fourier space array
	n = len(chain) / 2 - 1
	if not fit_length:
		fit_length = int(10.*fknee/df)
	
	#Construct the x-axis of the power spectrum
	j = slice(1,n+1)
	k = f[j]

	#Compute the power spectrum.
	Pj = (abs(fft(chain))**2)[j]
	
	#We fit to the log of the power spectrum, and only fit to the initial part of it
	#so as not to weight the white noise tail to much. 
	log_Pj = log(Pj[:fit_length])
	f=k[:fit_length]

	parameters = [fknee,log(sigma),alpha]
	residuals = _residuals(parameters,f,log_Pj)
	return (residuals**2).sum() / fit_length

def fitPowerSpectrum(tod,fsamp,fit_length=5000,cut_multiple=50):
	"""
	Fit a 1/f + white noise power spectrum to the power spectrum data, of the form
	P(f) = sigma^2 (1+(fknee/f)^alpha )
	The arguments are the TOD itself and the sampling frequency.
	Optional keywords are fit_length, a first iteration of the range of f to fit to (default 5000)
	and cut_multiple, which specifies how the chain is to be cut to improve speed (default 50)
	
	If this runs too slowly you could pass tod[::n],fsamp/n as the arguments for some thinning fraction n.
	
	The fit uses the algorithm described in astro-ph/0405462
	In particular you should not just fit least squares to the power spectrum in real space - 
	you have to do it in log space and include a constant offset.  Just read the paper.
	"""
	#Cut the chain down so it does not a large prime length (speeds things up)
	#If this runs too slowly it might be worth considering passing a thinned chain.
	chain=_cut_data(tod,cut_multiple)

	#Compute the spacing and sizes off the Fourier space array
	f=fftfreq(len(chain),1./fsamp)
	df = f[1]-f[0]
	n = len(chain) / 2 - 1
	
	#Construct the x-axis of the power spectrum
	j = slice(1,n+1)
	k = f[j]

	#Compute the power spectrum.
	Pj = (abs(fft(chain))**2)[j]
	
	#We fit to the log of the power spectrum, and only fit to the initial part of it
	#so as not to weight the white noise tail to much. 
	log_Pj = log(Pj[:fit_length])

	#Our starting points for the white noise and alpha are simple enough to choose.
	#For the fknee we start with a low guess so the value is pushed up.
	parameters_guess = array([
		0.01,               #fknee
		log(chain.std()),   #white noise level
		1.0                #alpha
		])
	
	#The frequencies of the regime we fit to.
	f=k[:fit_length]
	
	#First least squares fit, to a large length	
	bestfit,flag = leastsq(_residuals, parameters_guess,args=(f,log_Pj) )
	if flag not in [1,2,3,4]:
		raise ValueError("Could not fit spectrum")
		
	#We iterate once, using the 10 * the fknee as the upper limit for the max f we fit to.
	#The slow part of this code is the fft, so this is no extra cost.
	fit_length=int(10*bestfit[0]/df)
	log_Pj = log(Pj[:fit_length])
	f=k[:fit_length]
	bestfit,flag = leastsq(_residuals, bestfit,args=(f,log_Pj) )
	if flag not in [1,2,3,4]:
		raise ValueError("Could not fit spectrum second time")
	fk,log_w,alpha=bestfit
	return fk,exp(log_w),alpha
		
