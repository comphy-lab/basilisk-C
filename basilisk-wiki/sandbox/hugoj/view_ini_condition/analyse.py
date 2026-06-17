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
import scipy.integrate as integrate
import matplotlib as mpl

# add libpy
import os.path
import sys
dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, '../libpy/')
sys.path.append( filename )
from fftlib import get_spec_1D
from data_reader import read_bas_data
from tools_plots import plot_surface_3D
from solution_analytique import *

dopart = 1 # use 0 to do all

paths = {"mono":"./monochromatic/monoc.nc",
         "stokes":"./stokes/stokes.nc",
         "PM":"./synth_eta/synth_eta_0.nc"}

"""
## I. Mean current for different waves
"""
if dopart==1 or dopart==0:
    # linear (sin) wave
    H0 = 100        # m, depth of water
    dz = 0.5        # m
    a = 1           # m, amplitude of the wave
    L = 200.        # m, wavelength
    g = 9.81        # m.s-2
    k = 2*np.pi/L   # m-1, wavenumber
    
    omega = np.sqrt(g*k) # linear dispersion relation
    cphase = omega/k

    """
    We integrate numericaly the layer average from the linear solution
    """
    # Integration over 1 wavelength for a few layers
    # (see "Layer definitions")
    myX = np.arange(-L/2,L/2,L/1000)
    total_H = f_eta(myX, k, a) + H0
    H_frac = np.arange(0,1,0.05)
    z_layers = np.zeros((len(myX),len(H_frac)))
    # -> building an array with layer height(x)
    for l in range(len(H_frac)):
        z_layers[:,l] = total_H*H_frac[l]
    a_profile = np.zeros(len(H_frac))
    # -> we integrate
    borne1 = 0.
    borne2 = L
    for l in range(len(H_frac)):
        result, err = integrate.quad(lambda x: f_ul(x,d=(f_eta(x,k,a)+H0)*H_frac[l], k=k, a=a),
                                     borne1, borne2)
        result = 1/L * result
        if err>1e-6:
            print(result, err)
            raise Exception('Error in integration of function')
        a_profile[l] = result
    
    """
    Analytical solution
    """
    
    dzl = 0.05
    zl = np.arange(dzl,1.,dzl)
    a2_profile = np.zeros(zl.shape)
    for kz in range(len(zl)):
        a2_profile[kz] = analytical_u_profile_mode_m(n=10, 
                                                 sigma=zl[kz], 
                                                 t=0., 
                                                 omega=omega, 
                                                 m=1, 
                                                 H0=H0, 
                                                 L=L, 
                                                 a=a)

    """ 
    Verification that the serie converges to the numerical integration when
    we increase the number of terms for the mode m=1
    """
    # verification of recursion
    verif_first_terms_recursion()

    # nmax = 5
    # verif_analytical_solution(nmax)
    # plt.show()
    raise Exception

    """
    We compare with the initial field from the Basilisk simulation
    """
    kp_sth = 10*np.pi/L
    cphase_sth = np.sqrt(g*kp_sth)/kp_sth # using linear dispersion
    # monochromatic wave
    dsmono, grid = read_bas_data(paths["mono"])
    dsmono = dsmono.isel(time=0) # select the initial time
    bas_profile_u_m = dsmono['u.x'].mean(['x','y']).values
    bas_z_m = dsmono['z'].mean(['x','y']).values
    Nlayers = len(dsmono['zl'])
    # stokes wave
    ds_s, grid = read_bas_data(paths["stokes"])
    ds_s = ds_s.isel(time=0) # select the initial time
    bas_profile_u_s = ds_s['u.x'].mean(['x','y']).values
    bas_z_s = ds_s['z'].mean(['x','y']).values
    # synthetic wave field
    ds_PM, grid = read_bas_data(paths["PM"])
    ds_PM = ds_PM.isel(time=0) # select the initial time
    bas_profile_u_PM = ds_PM['u.x'].mean(['x','y']).values
    bas_z_PM = ds_PM['z'].mean(['x','y']).values

    """
    Plot the results for layer average
    """
    fig, ax = plt.subplots(1,1,figsize = (5,5),constrained_layout=True,dpi=100)
    ax.plot(a_profile/cphase,-np.mean(z_layers,axis=0)*k, label='analytical 1 linear mode', c='k')
    #TODO: analytical Stokes wave
    #TODO: analytical PM spectrum
    ax.scatter(bas_profile_u_m/cphase, bas_z_m*k, label=f'Basilisk 1 linear mode (nl={Nlayers})',c='g', marker='x')
    ax.scatter(bas_profile_u_s/cphase, bas_z_s*k, label=f'Basilisk Stokes wave (nl={Nlayers})',c='c', marker='x')
    ax.scatter(bas_profile_u_PM/cphase_sth, bas_z_PM*kp_sth, label=f'Basilisk PM Spectrum (nl={Nlayers})',c='orange', marker='s')
    ax.legend()
    ax.vlines(0,-H0,0,color='gray',ls='--')
    ax.set_xlabel('<u>/$c_p$ (m/s)')
    ax.set_ylabel(r'$z k_p$')
    ax.set_title(r'Layer average ($\sigma$ type)')
    ax.set_ylim([-H0*k, 0])
    fig.savefig('u_profiles_sigma.pdf')

    """
    3D surface field from basilisk
    """
    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100, subplot_kw={'projection': '3d'})
    plot_surface_3D(dsmono, 'u.x', cmin=-0.75, cmax=0.75,
                    zmin=-2,zmax=2, fig_tuple=(fig,ax), psave='3D_monochromatic_bas.pdf')

    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100, subplot_kw={'projection': '3d'})
    plot_surface_3D(ds_s, 'u.x',cmin=-0.75, cmax=0.75,
                    zmin=-2,zmax=2, fig_tuple=(fig,ax), psave='3D_stokes_bas.pdf')

    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100, subplot_kw={'projection': '3d'})
    plot_surface_3D(ds_PM, 'u.x',cmin=-0.75, cmax=0.75,
                    zmin=-25,zmax=25, fig_tuple=(fig,ax), psave='3D_PM0_bas.pdf')

    # """
    # Surface field analytical
    # """
    # fig, ax = plt.subplots(1,1,figsize = (5,5),constrained_layout=True,dpi=100)
    # ax.plot(myX, f_eta(myX, k, a), c='k', label='analytical')
    # for l in range(len(H_frac)):
    #     ax.plot(myX,z_layers[:,l]-H0, lw='0.5',c='k')
    # ax.scatter(dsmono.x.values, dsmono.eta[0,:],c='g',marker='x', label='Bas')
    # ax.set_xlim([-L/2,L/2])
    # ax.set_ylim([-50,2])
    # ax.legend()
    # ax.set_ylabel(r'$\eta$ (m)')
    # ax.set_xlabel("x (m)")
    # ax.set_title('Monochromatic wave')
    # ax.set_xlabel("x (m)")
    # ax.set_ylabel('depth (m)')

    """
    Illustration of the difference between layer definitions:
     i) fraction of total water height
     ii) constant distance from the surface
    """
    nl = 3
    H0_close = 10 # here I use 10m depth to exagerate the effect
    fig, ax = plt.subplots(1,1,figsize = (5,5),constrained_layout=True,dpi=100)
    ax.plot(myX, f_eta(myX, k, a), c='k', label='surface')
    ax.plot(myX, -np.ones(myX.shape)*H0_close, c='firebrick', label='flat bottom')
    total_H = f_eta(myX, k, a) + H0_close
    for l in range(1,nl):
        if l==1:
            ax.plot(myX, total_H*l/nl-H0_close, lw='0.5',c='k', label=r'fraction of $H_{total}$')
            ax.plot(myX, f_eta(myX, k, a)- H0_close*l/nl, lw='0.5',c='b', label=r'at $\eta$ - $d_i$')
        else:
            ax.plot(myX, total_H*l/nl-H0_close, lw='0.5',c='k')
            ax.plot(myX, f_eta(myX, k, a)- H0_close*l/nl, lw='0.5',c='b')
    ax.legend()
    ax.set_xlim([-L/2,L/2])
    ax.set_ylim([-H0_close,2])
    ax.set_title('Layer definitions')
    fig.savefig('layer_definitions.pdf')
    
    
    """
    Horizontal average at z = cst
    """
    # TODO: case with <u> at z=cst
    # z_havg = np.arange(0, -H0, -0.1)
    # if_dry1 = np.nan
    # if_dry2 = 0.
    # fig, ax = plt.subplots(1,1,figsize = (5,5),constrained_layout=True,dpi=100)
    



"""
# II. Direction of the synthetic wave field 


"""
if dopart==2 or dopart==0:
    path = "./synth_eta/"

    # a list of direction to test 
    list_theta = [0, np.pi, -np.pi, 2*np.pi]
    depth_max=-40. #m
    dds = []
    for k in range(len(list_theta)):
        dds.append(read_bas_data(path+"synth_eta_%d.nc" %k))

    clrs = ['b','g','orange','c']

    """
    ## Profiles with layer average
    """
    vmax = 0.07 # m/s
    fig, ax = plt.subplots(1,3,figsize = (10,5),constrained_layout=True,dpi=100)
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
        ds,grid = dds[k]

        fig, ax = plt.subplots(1,1,figsize = (6,5),constrained_layout=True,dpi=100)
        s = ax.pcolormesh(ds.x,ds.y,ds.eta.isel(time=0),cmap='Greys_r',vmin=-1,vmax=1)
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_title(fr'$\eta$ for $\theta_m$={k}$\pi$/4')
        plt.colorbar(s,ax=ax)
        fig.savefig(rf'eta_{k}pi.pdf')


plt.show()
