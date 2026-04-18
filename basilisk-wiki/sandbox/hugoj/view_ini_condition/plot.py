"""
# Investigation initial current profiles

The goal of this script is to understand why the profile of velocities under a
synthetic wave field is not 0

We generate a initial field using the spectrum.h lib and Basilisk code, then we
postprocess it using the present script 
"""

import numpy as np
from matplotlib import pyplot as plt
import xarray as xr

# add libpy
import os.path
import sys
dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, '../libpy/')
sys.path.append( filename )
from fftlib import get_spec_1D
from data_reader import read_bas_data



case = "mono"


"""
# Idealized case

## Monochromatic linear wave

"""
if case=="mono":
    # linear (sin) wave
    a = 1           # m, amplitude of the wave
    


"""
## Stokes wave
"""
if case=="Stokes":
    a=0




"""
# Synthetic wave field case

## General parameters

"""
if case=="synth":
    path = "./gen_ini/"

# a list of direction to test 
    list_theta = [0, np.pi, -np.pi, 2*np.pi]
    depth_max=-40. #m
    dds = []
    for k in range(len(list_theta)):
        dds.append(read_bas_data(path+"out_%d.nc" %k))

    clrs = ['b','g','orange','c']

    """
    ## Profiles with layer average
    """
    vmax = 0.05 # m/s
    fig, ax = plt.subplots(1,4,figsize = (10,5),constrained_layout=True,dpi=100)
    for k in range(len(list_theta)):
        ds,grid = dds[k]
        
        meanU=ds['u.x'].mean(['x','y','time'])
        meanV=ds['u.y'].mean(['x','y','time'])
        meanZ=ds['z'].mean(['x','y','time'])
        meanW=ds['u.z'].mean(['x','y','time'])

        # U
        ax[0].plot(meanU, meanZ, c=clrs[k], label=rf"$\theta_m$={k}$\pi$",
                marker='x')
        ax[0].set_ylim([depth_max, 0.])
        ax[0].set_xlim([-vmax,vmax])
        ax[0].set_xlabel('<U> (m/s)')
        # V
        ax[1].plot(meanV, meanZ, c=clrs[k],
                marker='x')
        ax[1].set_ylim([depth_max, 0.])
        ax[1].set_xlim([-vmax,vmax])
        ax[1].set_xlabel('<V> (m/s)')
        # W
        ax[2].plot(meanW, meanZ, c=clrs[k],
                marker='x')
        ax[2].set_ylim([depth_max, 0.])
        ax[2].set_xlim([-vmax,vmax])
        ax[2].set_xlabel('<W> (m/s)')

    ax[0].set_ylabel('depth (m)')
    ax[0].legend(loc='lower left')
    fig.savefig('velocity_profiles.svg')



    """
    ## Sea surface state
    """

    for k in range(len(list_theta)):
        fig, ax = plt.subplots(1,1,figsize = (6,5),constrained_layout=True,dpi=100)
        s = ax.pcolormesh(ds.x,ds.y,ds.eta.isel(time=0),cmap='Greys_r',vmin=-1,vmax=1)
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_title(fr'$\eta$ for $\theta_m$={k}$\pi$')
        plt.colorbar(s,ax=ax)
        fig.savefig(rf'eta_{k}pi.pdf')


#plt.show()
