"""
# Description

    This script compares a Basilisk multi-layer simulation to [1,2].

[1] Wu, J., Popinet, S., & Deike, L. (2023). Breaking wave field statistics 
with a multi-layer model. Journal of Fluid Mechanics, 968, A12. 
https://doi.org/10.1017/jfm.2023.522

[2] Wu, J., Popinet, S., Chapron, B., Farrar, J. T., & Deike, L. (2025). 
Turbulence and Energy Dissipation from Wave Breaking. Journal of Physical 
Oceanography, 55(9), 1521–1534. https://doi.org/10.1175/JPO-D-25-0052.1

# Conda Environment

    conda activate ml_dec2025

# How to run

    python3 visualizer.py
"""
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import time
import xrft
from scipy.fft import fft2, fftfreq

# add libpy
import os.path
import sys
dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, '../../libpy/')
sys.path.append( filename )
from fftlib import *
from data_reader import read_data, build_grid
from tools import *
from diags import interpz, grad_velocities, vorticity, dissipation

"""
## General parameters
"""
# My data
#filename="/home/jacqhugo/basilisk/wiki/sandbox/hugoj/reproducing_jiarongs_plots/N10_P0.02_L15/out.nc"
#filename="/home/jacqhugo/basilisk/wiki/sandbox/hugoj/my_workflow_reproducing_jiarongs_plots/N10_L15_P0.02/out.nc"
filename="/home/jacqhugo/basilisk/wiki/sandbox/hugoj/my_workflow_reproducing_jiarongs_plots/using_pprf90/out.nc"
# getting back Jiarong's data
Jpath = "data_Jiarong/"
save_nc = './data.nc'

dpi=200

g = 9.81
L0 = 200.0 # m
kp = 10  * np.pi / L0
wp = np.sqrt(g*kp)
fp = wp/(2*np.pi)
Tp = 1/fp
N = 1024
L = 200.

print("--------------------")
print("Peak values:")
print("kp",kp)
print("wp",wp)
print("fp",fp)
print("Tp",Tp)
print("--------------------\n")


if False:
    print("* Surface elevation, spectra")

    """
    ## Surface elevation
    """

    # opening file
    ds, grid = read_data(filename, chunks=None)
    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=dpi)
    s=ax.pcolormesh(ds.x, ds.y, ds.eta.isel(time=-1)*kp, cmap='Greys_r', vmin=-0.15, vmax=0.3)
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    plt.colorbar(s,ax=ax,orientation='vertical',label=r'$\eta k_p$')
    fig.savefig('eta.png')
    

    Hs = 4*np.sqrt((ds.eta**2).mean(dim=['x','y']))
    print('Depth always under the surface: %f (m)' %ds.z.isel(zl=-1).min().values)
    print('kpHs ini = %f' %(kp*Hs[0]))
    print('Hs ini (m)= %f' %Hs[0])

    """
    ## Eta spectrum
    """

    # Jiarong's data
    Js_0 = np.loadtxt(Jpath+'2023_fig3c_0.txt',skiprows=1,delimiter=",")
    Js_20 = np.loadtxt(Jpath+'2023_fig3c_25.txt',skiprows=1,delimiter=",")
    Js_100 = np.loadtxt(Jpath+'2023_fig3c_124.txt',skiprows=1,delimiter=",")
    Js_120 = np.loadtxt(Jpath+'2023_fig3c_149.txt',skiprows=1,delimiter=",")
    Js = [Js_0, Js_20, Js_100, Js_120]
    at_t = [0, 20, 100, 120]

    # building a colormap
    colors = plt.get_cmap('plasma')(np.linspace(0, 1, len(at_t)))
    
    """
    Verification of Parceval equalities
    """
    if True:
        it = -1
        print("\n* Validation of the perceval theorem")
        print("timestep = %d" %it)
        
        print('- 1) using xrft')
        s_eta = xrft.isotropic_power_spectrum(ds.eta[it], dim=('x','y'), truncate=True)
        s_eta2D = xrft.xrft.cross_spectrum(ds.eta[it], ds.eta[it],dim=('x','y'))
        N = len(ds.eta.x)
        dk =(s_eta.freq_r[1] - s_eta.freq_r[0]).values
        dk = np.zeros(len(s_eta.freq_r))
        dk[0] = s_eta.freq_r[0].values
        dk[1:] = s_eta.freq_r[1:].values-s_eta.freq_r[:-1].values
        dkx = (s_eta2D.freq_x[1]-s_eta2D.freq_x[0]).values
        dky = (s_eta2D.freq_y[1]-s_eta2D.freq_y[0]).values
        var_sum = np.sum(ds.eta[it].values**2)/N**2
        var_integ1D = np.sum(s_eta.values*dk)
        var_2D = np.sum(s_eta2D.values)*dkx*dky
        print(f'Target: variance (sum of eta**2) = {var_sum}')
        print(f'variance (sum of F(kx,ky)dkxdky) = {var_2D}')
        print(f'variance (sum of F(k)dk) = {var_integ1D*2*np.pi}')
        
        print('- 2) using geostrokit')
        k,l, spec2D = get_spec_2D(ds.eta[it].values, ds.eta[it].values, Delta=200/1024)
        kr, spec_1D = get_spec_1D(ds.eta[it].values, ds.eta[it].values, Delta=200/1024)
        print(f'Target: variance (sum of eta**2) = {var_sum}')
        print(f'variance spectr1D = {np.sum(spec_1D*(kr[1]-kr[0]))}')
        print(f'variance spectr2D = {np.sum(spec2D*(k[0,1]-k[0,0])*(l[1,0]-l[0,0]))}')

        print("- 3) using Jiarong's code")
        print("     var(eta), sum(F), sum(F_polar)")
        
        F_150 = jspectrum_integration(ds.eta[it].values, L0, N, CHECK=True)
    
    if True:
        print('Evolution of wave energy <eta**2> ')
        for k in range(len(ds.eta.time)):
            print(f"t = {ds.time[k].values}, {np.sum(ds.eta[k].values**2)/N**2} m^2")

    """
    Spectrum evolution in time
    """
    # Computing spectrum at all timesteps
    # 1) xrft
    s_eta = xrft.isotropic_power_spectrum(ds.eta, dim=('x','y'), truncate=True)
    s_eta2D = xrft.xrft.cross_spectrum(ds.eta, ds.eta,dim=('x','y'))
    # 2) geostrokit
    f_eta = np.zeros((len(ds.eta.time),N//2-1))
    for k in range(len(ds.eta.time)):
        fr, f_eta[k] = get_spec_1D(ds.eta[k].values, ds.eta[k].values, Delta=L/N)
    # 3) Jiarong's code
    wavenumber = 2*np.pi*np.fft.fftfreq(n=N,d=L0/N)
    jiarong_k = wavenumber[0:int(N/2)]
    f_eta_jiarong = np.zeros((len(ds.eta.time),len(jiarong_k)))
    for k in range(len(ds.eta.time)):
        f_eta_jiarong[k] = jspectrum_integration(ds.eta[k].values, L0, N, CHECK=False)



    # Plotting
    fig, ax = plt.subplots(1,1,figsize = (3,3),constrained_layout=True,dpi=dpi)
    for k in range(len(at_t)):
        ax.loglog(fr*2*np.pi*L0, (f_eta[k]/(2*np.pi))*kp**3,
                    c=colors[k],
                label=r'$w_p t=$'+str(int(wp*at_t[k])))
        
    ax.set_xlabel(r'$kL$')
    ax.set_ylabel(r'$F_{\eta}(k).k_p^3$')
    ax.vlines(1,1e-8,1,ls='--', colors='gray') # TODO: modify this into 1/dx
    ax.vlines(kp*L0,1e-8,1,ls='--', colors='gray')
    ax.set_ylim([1e-7, 1e-2])
    ax.set_xlim([1e1,1e3])
    ax.legend()
    fig.savefig('eta_spectra_evolution.svg')
    fig.savefig('eta_spectra_evolution.png')

    # vs Jiarong
    for k,ttime in enumerate(at_t):
        fig, ax = plt.subplots(1,1,figsize = (3,3),constrained_layout=True,dpi=dpi)
        ax.loglog(jiarong_k*L0, (f_eta_jiarong[k])*kp**3,
                c='gray',
                label="My c code, Jiarong py")
        ax.loglog(fr*2*np.pi*L0, (f_eta[k]/(2*np.pi))*kp**3,
                c='b', ls='--',
                  label='my code code only') # using geostrokit
        ax.loglog(Js[k][:,0],Js[k][:,1], c='r', ls='-', label="Jiarong's paper")
        ax.set_title(r'$w_p t=$'+str(int(wp*at_t[k])))
        ax.set_ylim([1e-7,1e-2])
        ax.set_xlim([1e1,1e3])
        ax.set_xlabel(r'$kL$')
        ax.set_ylabel(r'$F_{\eta}(k).k_p^3$')
        plt.legend(fontsize=6)
        fig.savefig(r"eta_spectra_vs_J_%d.svg" % int(wp*at_t[k]))
        fig.savefig(r"eta_spectra_vs_J_%d.png" % int(wp*at_t[k]))
        

"""
## Profiles
"""

if True:
    print("* Profiles")
        
    time_ =  120 # s
    znew = np.arange(-8,0.1,0.1) # cartesian Z to interpolate on
    
    file_grads = ['dudz.nc','dudy.nc','dudx.nc','dvdz.nc','dvdy.nc','dvdx.nc','dwdz.nc','dwdy.nc', 'dwdx.nc']
    file_omeg = ['omegaxp.nc','omegayp.nc','omegazp.nc','enstrophy.nc']

    """### We first need to compute the required quantities: gradients, omega and dissipation"""
    # Note about performances (memory and compute):
    # We compute in a open-compute-close manner to avoid memory overload,
    # in combination with dask's chunks.

    chunks = {'time':5, 'x':64, 'y':64, 'zl':-1}
    # We need the grid for the following
    ds, grid = read_data(filename, chunks=chunks)    

    # main file
    if not os.path.isfile(save_nc):
        print('- I need to save the data in a proper format first !')
        ds.to_netcdf(save_nc)
    ds.close()

    # gradients
    if not os.path.isfile("dudz.nc"):
        print('- I need to compute grads ...')
        ds = xr.open_dataset(save_nc, chunks=chunks)
        ds, update = grad_velocities(ds, grid)
        if update:
            for name in file_grads:
                ds[name[:-3]].to_netcdf(name)
        ds.close()
    
    # enstrophy
    if not os.path.isfile("enstrophy.nc"):  
        print('- I need to compute enstrophy ...')
        ds = xr.open_mfdataset(file_grads+[save_nc], chunks=chunks)
        ds, update = vorticity(ds, grid)
        if update:
            for name in file_omeg:
                ds[name[:-3]].to_netcdf(name)
        ds.close() 
    
    # Dissipation
    if not os.path.isfile("epsilon.nc"):  
        print('- I need to compute dissipation ...')
        ds = xr.open_mfdataset(file_grads+[save_nc], chunks=chunks)
        ds, update = dissipation(ds, grid)
        if update:
            for name in ['epsilon']:
                ds[name].to_netcdf(name+'.nc')
        ds.close() 

    # We have all the pieces ! lets open everything
    print('- Required data: computed')
    ds = xr.open_mfdataset([save_nc]+file_grads+file_omeg+["epsilon.nc"], chunks=chunks)
    #ds = ds.chunk({'time':5, 'x':64, 'y':64, 'zl':-1})  # Dask performances
    
    # initial profiles
    ds0 = ds.sel(time=0)
    ux_lagr0 = ds0['u.x'].mean(['x','y']).compute()
    ens_lagr0 = ds0['enstrophy'].mean(['x','y']).compute()
    diss_lagr0 = ds0['epsilon'].mean(['x','y']).compute()
    z_lagr0 = ds0.z.mean(['x','y']).compute()
    

    ds1 = ds.sel(time=time_) # select the time at which to plot profiles
    

    # print('starting interp1')
    # t1 = time.time()
    # print(ds.z.isel(time=-1))
    # print(ds['u.x'].isel(time=-1))
    # u_interp = interpz(ds.z.isel(time=-1), ds['u.x'].isel(time=-1), znew).compute()
    # u_interp = u_interp.assign_coords(znew=znew) 
    # t2 = time.time()
    # print(f'interp1 done, {int(t2-t1)} s')

    """We can now compute the different avg methods"""

    z_lagr = ds1.z.mean(['x','y'])
    ux_interp1 = interpz(ds1.z, ds1['u.x'], znew, fill_value=np.nan).mean(dim=['x','y']).compute()
    ux_interp2 = interpz(ds1.z, ds1['u.x'], znew, fill_value=0.).mean(dim=['x','y']).compute()
    ux_lagr = ds1['u.x'].mean(['x','y'])

    ens_interp1 = interpz(ds1.z, ds1['enstrophy'], znew, fill_value=np.nan).mean(dim=['x','y']).compute()
    ens_interp2 = interpz(ds1.z, ds1['enstrophy'], znew, fill_value=0.).mean(dim=['x','y']).compute()
    ens_lagr = ds1['enstrophy'].mean(['x','y'])

    diss_interp1 = interpz(ds1.z, ds1['epsilon'], znew, fill_value=np.nan).mean(dim=['x','y']).compute()
    diss_interp2 = interpz(ds1.z, ds1['epsilon'], znew, fill_value=0.).mean(dim=['x','y']).compute()
    diss_lagr = ds1['epsilon'].mean(['x','y'])

    # z_lagr = ds.z.sel(time=time_).mean(['x','y'])
    # ux_interp1 = interpz(znew, ds1.z, ds1['u.x'], fill_value=np.nan).mean(dim=['x','y']).compute()
    # ux_interp2 = interpz(znew, ds1.z, ds1['u.x'], fill_value=0.).mean(dim=['x','y'])
    # ux_lagr = ds1['u.x'].mean(['x','y'])
    #
    # ens_interp1 = interpz(znew, ds1.z, ds1['enstrophy'], fill_value=np.nan).mean(dim=['x','y'])
    # ens_interp2 = interpz(znew, ds1.z, ds1['enstrophy'], fill_value=0.).mean(dim=['x','y'])
    # ens_lagr = ds1['enstrophy'].mean(['x','y'])
    #
    # diss_interp1 = interpz(znew, ds1.z, ds1['epsilon'], fill_value=np.nan).mean(dim=['x','y'])
    # diss_interp2 = interpz(znew, ds1.z, ds1['epsilon'], fill_value=0.).mean(dim=['x','y'])
    # diss_lagr = ds1['epsilon'].mean(['x','y'])
    #


    # Jiarong's data
    Jux_interp1 = np.loadtxt(Jpath+'2025_fig7a_abs1.txt',skiprows=1,delimiter=",")
    Jux_interp2 = np.loadtxt(Jpath+'2025_fig7a_abs2.txt',skiprows=1,delimiter=",")
    Jux_lagr = np.loadtxt(Jpath+'2025_fig7a_layer.txt',skiprows=1,delimiter=",")
    
    Jens_interp1 = np.loadtxt(Jpath+'2025_fig7b_abs1.txt',skiprows=1,delimiter=",")
    Jens_interp2 = np.loadtxt(Jpath+'2025_fig7b_abs2.txt',skiprows=1,delimiter=",")
    Jens_lagr = np.loadtxt(Jpath+'2025_fig7b_layer.txt',skiprows=1,delimiter=",")

    Jdiss_interp1 = np.loadtxt(Jpath+'2025_fig7c_abs1.txt',skiprows=1,delimiter=",")
    Jdiss_interp2 = np.loadtxt(Jpath+'2025_fig7c_abs2.txt',skiprows=1,delimiter=",")
    Jdiss_lagr = np.loadtxt(Jpath+'2025_fig7c_layer.txt',skiprows=1,delimiter=",")
    
    """
    ### Profiles from my data
    """
    if True:
        fig, ax = plt.subplots(1,3,figsize = (9,3),constrained_layout=True,dpi=dpi)
        # U
        ax[0].set_xlabel('<u> (m/s)')
        ax[0].plot(ux_lagr0, z_lagr0, c='gray', ls='-', label='layer t=0',  marker='s', markerfacecolor='None')
        ax[0].plot(ux_lagr, z_lagr, c='k', ls='-', label='layer',  marker='s', markerfacecolor='None')
        ax[0].plot(ux_interp1, znew, c='b', ls='-', label='abs 1')
        ax[0].plot(ux_interp2, znew, c='b', ls='--', label='abs 2')
        ax[0].set_xlim([0,0.2])

        # Vorticity
        ax[1].set_xlabel('<enstrophy>')
        ax[1].semilogx(ens_lagr0, z_lagr0, c='gray', ls='-', label='layer t=0',  marker='s', markerfacecolor='None')
        ax[1].semilogx(ens_lagr, z_lagr, c='k', ls='-', label='layer',  marker='s', markerfacecolor='None')
        ax[1].semilogx(ens_interp1, znew, c='b', ls='-', label='abs 1')
        ax[1].semilogx(ens_interp2, znew, c='b', ls='--', label='abs 2')
        ax[1].set_xlim([1e-4,1])

        # Dissipation
        ax[2].set_xlabel('<diss>')
        ax[2].semilogx(diss_lagr0, z_lagr0, c='gray', ls='-', label='layer t=0', marker='s', markerfacecolor='None')
        ax[2].semilogx(diss_lagr, z_lagr, c='k', ls='-', label='layer', marker='s', markerfacecolor='None')
        ax[2].semilogx(diss_interp1, znew, c='b', ls='-', label='abs 1')
        ax[2].semilogx(diss_interp2, znew, c='b', ls='--', label='abs 2')
        ax[2].set_xlim([1e-3,10])

        for axe in ax:
            axe.legend(frameon=False)
            axe.set_ylim([-10,0])
            axe.set_ylabel('z (m)')
        fig.savefig('avg_profiles.svg')
        fig.savefig('avg_profiles.png')

    """
    ### Profiles comparison vs Jiarong's papers
    """
    if True:
        fig, ax = plt.subplots(1,3,figsize = (9,3),constrained_layout=True,dpi=dpi)
        
        # U
        ax[0].set_xlabel('<u> (m/s)')
        ax[0].plot(Jux_lagr[:,0], Jux_lagr[:,1], c='gray', label='J layer', marker='x', markerfacecolor='None')
        ax[0].plot(ux_lagr, z_lagr, c='k', ls='-', label='layer', marker='s', markerfacecolor='None')
        ax[0].set_xlim([0,0.2])
        # ax[0].plot(ux_interp1, znew, c='b', ls='-', label='abs 1')
        # ax[0].plot(ux_interp2, znew, c='b', ls='--', label='abs 2')
        
        # Vorticity
        ax[1].set_xlabel('<enstrophy>')
        ax[1].semilogx(Jens_lagr[:,0], Jens_lagr[:,1], c='gray', label='J layer', marker='x', markerfacecolor='None')
        ax[1].semilogx(ens_lagr, z_lagr, c='k', ls='-', label='layer', marker='s', markerfacecolor='None')
        ax[1].set_xlim([1e-4,1])
        # ax[1].semilogx(ens_interp1, znew, c='b', ls='-', label='abs 1')
        # ax[1].semilogx(ens_interp2, znew, c='b', ls='--', label='abs 2')
        
        # Dissipation
        ax[2].set_xlabel('<diss>')
        ax[2].semilogx(Jdiss_lagr[:,0], Jdiss_lagr[:,1], c='gray', label='J layer', marker='x', markerfacecolor='None' )
        ax[2].semilogx(diss_lagr, z_lagr, c='k', ls='-', label='layer', marker='s', markerfacecolor='None')
        # ax[2].semilogx(diss_interp1, znew, c='b', ls='-', label='abs 1')
        # ax[2].semilogx(diss_interp2, znew, c='b', ls='--', label='abs 2')
        ax[2].set_xlim([1e-3,10])
        
        # nicer figures
        ax[1].set_title(r'$t=%d s$'%time_)
        for axe in ax:
            axe.legend(frameon=False)
            axe.set_ylim([-10,0])
            axe.set_ylabel('z (m)')
        fig.savefig('avg_profiles_vsJiarong.png')
        fig.savefig('avg_profiles_vsJiarong.svg')


plt.show()

