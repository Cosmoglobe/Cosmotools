import numpy as np
import healpy as hp

# Simple importable function which simulates a noise realization
# given some passed in rms map. Input and output are strings pointing
# to the files to read in and write out

def sim_noise_map(map_in: str,map_out: str):
    rms_map     = hp.read_map(map_in)
    npix        = np.shape(rms_map)[0]
    realization = np.zeros(npix)
    for i in range(npix):
        realization[i] = np.random.normal(0.0,rms_map[i])

    hp.write_map(map_out,realization)
