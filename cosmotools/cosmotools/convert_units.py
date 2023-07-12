import numpy as np
import astropy.units as u
from astropy.constants import c
from astropy.constants import h
from astropy.constants import k_B
#########################################################################################
# This file is written to import functions making conversions from one microwave/sub-mm
# data unit to another. The code is written to handle bandpasses additionally. 
# The three unit types taken into account here are:
#
# a - Antenna units     [uK_RJ]
# f - Flux units        [MJy/sr]
# t - Temperature units [uK_cmb]

# WARNING:
# The unit conversions are only supported in frequency space
# and the frequency passed in is assumed to be in GHz. Passing
# Hz (for example) will yield incorrect conversion factors.
#########################################################################################

T_cmb = 2.7255*u.K # Fixsen et al 1996


def convert_help():
    thegoods = """List of convert_units.py functions:
                  compute_bnu_prime
                  compute_bnu_prime_RJ
                  a2f (Converts from antenna to flux)
                  a2t (Converts from antenna to temperature(cmb))
                  f2t (Converts from flux to temperature)"""

    print(thegoods)
    
# The derivative of the Planck function wrt the temperature
def compute_bnu_prime(nu):
    y = (h*(nu*u.GHz)).to(u.J)/(k_B*T_cmb) 

    item = (2.0*h*(nu*u.GHz)**3)/(c**2*(np.exp(y)-1))*(np.exp(y)/(np.exp(y)-1))*h*(nu*u.GHz)/(k_B*T_cmb**2)

    return item.to(u.MJy/u.uK)

# Multiply this by some K_RJ value to get Flux Density
def compute_bnu_prime_RJ(nu):
    item =  (2.0*k_B*(nu*u.GHz)**2)/c**2

    return item.to(u.MJy/u.uK)/u.sr

# [uK_cmb/uK_RJ]
def a2t(nu,bp='none'):
    if bp == 'none':
        y = (h*(nu*u.GHz).to(1/u.s))/(k_B*T_cmb)
        convert = (np.exp(y)-1.0)**2/(y**2*np.exp(y)) 
    else:
        bandpass = np.loadtxt(bp)
        xs       = bandpass.T[0]
        weights  = bandpass.T[1]
        # Normalize the bandpass
        weights  = weights/np.sum(weights)
        nbp = len(xs)

        convert = 0
        for i in range(nbp):
            y = (h*(xs[i]*u.GHz).to(1/u.s))/(k_B*T_cmb)
            convert += weights[i]*(np.exp(y)-1.0)**2/(y**2*np.exp(y)) 
        
    return convert

# [MJy/sr / uK_RJ] Tested and confirmed
def a2f(nu,bp='none'):
    if bp == 'none':
        return compute_bnu_prime_RJ(nu)
    else:
        bandpass = np.loadtxt(bp)
        xs       = bandpass.T[0]
        weights  = bandpass.T[1]
        # Normalize the bandpass
        weights  = weights/np.sum(weights)
        nbp = len(xs)

        convert = 0
        for i in range(nbp):
            convert += weights[i]*compute_bnu_prime_RJ(xs[i])

    return convert

# [uK_cmb / MJy/sr]
def f2t(nu,bp='none'):
    if bp == 'none':
        return 1.0/(compute_bnu_prime(nu))
    else:
        bandpass = np.loadtxt(bp)
        xs       = bandpass.T[0]
        weights  = bandpass.T[1]
        # Normalize the bandpass
        weights  = weights/np.sum(weights)
        nbp = len(xs)

        convert = 0
        for i in range(nbp):
            if weights[i] < 1e-10:
                continue
            convert += weights[i]/(compute_bnu_prime(xs[i]))
    return convert
