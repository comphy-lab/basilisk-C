"""
# My C code vs Jiarong's C code

I compare the results from my C code (my workflow, but same parameters) and the original code from Jiarong Wu
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# add libpy
import os.path
import sys
dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, '../../libpy/')
sys.path.append( filename )
from fftlib import *
from diags import interpz

"""
## Some parameters
"""
g = 9.81
L0 = 200.0 # m
kp = 10  * np.pi / L0
wp = np.sqrt(g*kp)
fp = wp/(2*np.pi)
Tp = 1/fp
N = 1024
L = 200.

dpi=200


# My data
pathH_C = './using_pprC/'
pathH_f = './using_pprf90/'
pathJ = "../../reproducing_jiarongs_plots/postprocess/"
filedata = "data.nc"
fileJ = pathJ + filedata

file_grads = ['dudz.nc','dudy.nc','dudx.nc','dvdz.nc','dvdy.nc','dvdx.nc','dwdz.nc','dwdy.nc', 'dwdx.nc']
file_omeg = ['omegaxp.nc','omegayp.nc','omegazp.nc','enstrophy.nc']
file_diss = ['epsilon.nc']


# Jiarong's paper data
Jpath = "data_Jiarong/"
Js_0 = np.loadtxt(Jpath+'2023_fig3c_0.txt',skiprows=1,delimiter=",")
Js_20 = np.loadtxt(Jpath+'2023_fig3c_25.txt',skiprows=1,delimiter=",")
Js_100 = np.loadtxt(Jpath+'2023_fig3c_124.txt',skiprows=1,delimiter=",")
Js_120 = np.loadtxt(Jpath+'2023_fig3c_149.txt',skiprows=1,delimiter=",")
Js = [Js_0, Js_20, Js_100, Js_120]
at_t = [0, 20, 100, 120]


# We build a dict, easier to process
d_data = {"Hugo_pprC":xr.open_mfdataset([pathH_C+file for file in ([filedata]+file_grads+file_omeg+file_diss)]),
          "Hugo_pprf90":xr.open_mfdataset([pathH_f+file for file in ([filedata]+file_grads+file_omeg+file_diss)]),
          "Jiarong":xr.open_mfdataset([fileJ] + [pathJ+file for file in (file_grads+file_omeg+file_diss)])}


"""
## Variance evolution with time
"""
if False:
    print('\nVariance check')
    clrs = ['b','g','orange']
    fig, ax = plt.subplots(1,1,figsize = (3,3),constrained_layout=True,dpi=dpi)
    for k,name in enumerate(d_data.keys()):
        print('-> '+name)
        var = np.zeros(len(at_t))
        for it,time_ in enumerate(at_t):
            ds=d_data[name].sel(time=time_)
            var[it] = (ds.eta**2).mean(['x','y']).values
            print('time=%f: %f m^2' %(ds.time,var[it]))
            
        ax.scatter(at_t/Tp, var/var[0], marker='x', color=clrs[k], label=name)
    ax.legend()
    ax.set_xlabel('t/Tp')
    ax.set_ylabel(r'E/E0')
    fig.savefig('comp_wave_energy.svg')


"""
## Spectrum evolution
"""
if False:
    print('\nSpectrum evolution')
    clrs = ['b','g','orange']
    for it,time_ in enumerate(at_t):
        fig, ax = plt.subplots(1,1,figsize = (3,3),constrained_layout=True,dpi=dpi)
        ax.loglog(Js[it][:,0],Js[it][:,1], c='gray', ls='-', label="paper")
        for k,name in enumerate(d_data.keys()):
            ds=d_data[name].sel(time=time_)
            kr, F_k = get_spec_1D(ds.eta.values, ds.eta.values, Delta=L0/N)
            ax.loglog(kr*2*np.pi*L0, (F_k/(2*np.pi)*kp**3), label=name, c=clrs[k])
        ax.vlines(kp*L0,1e-8,1,ls='--', colors='gray')
        ax.set_ylim([1e-7, 1e-2])
        ax.set_xlim([1e1,1e3])
        ax.legend()
        ax.set_title(f'time = {time_} s ($w_p t=$'+str(int(wp*time_))+')')
        fig.savefig(f'comp_spec_t{time_}.svg')

"""
## Average profiles
"""
if True:
    print('\nProfiles ')
    time_=120 # s
    clrs = ['b','g','orange']
    # Jiarong's data
    Jux_lagr = np.loadtxt(Jpath+'2025_fig7a_layer.txt',skiprows=1,delimiter=",")
    Jens_lagr = np.loadtxt(Jpath+'2025_fig7b_layer.txt',skiprows=1,delimiter=",")
    Jdiss_lagr = np.loadtxt(Jpath+'2025_fig7c_layer.txt',skiprows=1,delimiter=",")

    fig, ax = plt.subplots(1,3,figsize = (9,3),constrained_layout=True,dpi=dpi)
    ax[0].plot(Jux_lagr[:,0], Jux_lagr[:,1], c='gray', label='paper', marker='x', markerfacecolor='None')
    ax[1].semilogx(Jens_lagr[:,0], Jens_lagr[:,1], c='gray', label='paper', marker='x', markerfacecolor='None')
    ax[2].semilogx(Jdiss_lagr[:,0], Jdiss_lagr[:,1], c='gray', label='paper', marker='x', markerfacecolor='None' )

    for k,name in enumerate(d_data.keys()):
        ds1=d_data[name].sel(time=time_)
        color = clrs[k]

        z_lagr = ds1.z.mean(['x','y'])
        ux_lagr = ds1['u.x'].mean(['x','y'])
        ens_lagr = ds1['enstrophy'].mean(['x','y'])
        diss_lagr = ds1['epsilon'].mean(['x','y'])

        # U
        ax[0].plot(ux_lagr, z_lagr, c=color, ls='-', label=name,  marker='s', markerfacecolor='None')
        ax[0].set_xlim([0,0.2])

        # Vorticity
        ax[1].semilogx(ens_lagr, z_lagr, c=color, ls='-', label=name,  marker='s', markerfacecolor='None')
        ax[1].set_xlim([1e-4,1])

        # Dissipation
        ax[2].semilogx(diss_lagr, z_lagr, c=color, ls='-', label=name, marker='s', markerfacecolor='None')
        ax[2].set_xlim([1e-3,10])

    ax[0].set_xlabel('<u> (m/s)')
    ax[1].set_xlabel('<enstrophy>')
    ax[2].set_xlabel('<diss>')
    ax[1].set_title(fr'time = {time_} s ($\omega_p t=$'+str(int(wp*time_))+')')
    for axe in ax:
        axe.legend(frameon=False)
        axe.set_ylim([-10,0])
        axe.set_ylabel('z (m)')
    fig.savefig('comp_avg_profiles.svg')


"""
## Enstrophy differences between PPR lib and remapc.c
"""
if False:
    print('\nEnstrophy maps ')
    time_ = 120 # s
    atz = -9. #m
    znew = np.ones(1)*atz
    fig, ax = plt.subplots(1,2,figsize = (7,3),constrained_layout=True,dpi=dpi)
    for k,name in enumerate(["Hugo_pprC","Hugo_pprf90"]):
        ds1=d_data[name].sel(time=time_)
        e = interpz(ds1.z, ds1['enstrophy'], znew, fill_value=np.nan).compute()[:,:,0]
        s=ax[k].pcolormesh(ds1.x,ds1.y,e,cmap='plasma', norm=colors.LogNorm(vmin=1e-5, vmax=1e-2))
        ax[k].set_title(name)
        ax[k].set_aspect(1.)
        ax[k].set_xlabel('X (m)')
        ax[k].set_ylabel('Y (m)')
    plt.colorbar(s,ax=ax)
    fig.savefig(f'enstrophy_maps_ppr_vs_remapc_z{atz}_t{time_}.png')
plt.show()




