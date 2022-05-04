import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

# Simple, importable function which returns a slope, offset, and the
# 2x2 correlation matrix between two maps, as well as a scatter plot
# with linear fit over-plotted.

def correlate_maps(filename1:str,filename2:str,fileout:str,maskfile='none'):

    map1 = hp.read_map(filename1)
    map2 = hp.read_map(filename2)

    npix = np.shape(map1)[0]

    if maskfile == 'none':
        m, b  = np.polyfit(map1,map2,1)
        model = m*map1+b
        rho   = np.corrcoef(map1,map2)

        print(f'Slope = {m}, offset = {b}')
        print(f'rho   = {rho}')

        x = np.linspace(np.min(map1),np.max(map1),1000)

        plt.figure(figsize=(6,4))
        plt.scatter(map1,map2,color='k')
        plt.plot(x,m*x+b,color='red',linestyle='--')
        plt.savefig(fileout,dpi=100,bbox_inches='tight')
        plt.show()
    else:
        mask = hp.read_map(maskfile)

        map1masked = []
        map2masked = []
    
        for i in range(npix):
            if mask[i] == hp.UNSEEN or mask[i] == 0.0:
                map1[i] = 0.0
                map2[i] = 0.0
            else:
                map1masked.append(map1[i])
                map2masked.append(map2[i])

        m, b  = np.polyfit(map1masked,map2masked,1)
        rho   = np.corrcoef(map1masked,map2masked)

        print(f'Slope = {m}, offset = {b}')
        print(f'rho   = {rho}')

        x = np.linspace(np.min(map1masked),np.max(map1masked),1000)
    
        plt.figure(figsize=(6,4))
        plt.scatter(map1masked,map2masked,color='k')
        plt.plot(x,m*x+b,color='red',linestyle='--')
        plt.savefig(fileout,dpi=100,bbox_inches='tight')
        plt.show()
