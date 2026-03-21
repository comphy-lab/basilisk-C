# script version of the jupyter notebook from Jiarong:
# https://github.com/jiarong-wu/multilayer_breaking/blob/main/specgen/Specgen_example.ipynb


import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from specgen import spectrum_PM, spectrum_gen_linear, eta_random
import xarray as xr
import xrft
"""
Generate kx-ky spectrum from a given unidirectional spectrum shape

We demonstrate with the PM spectrum shape
"""

''' Wrap the function around so that it only takes in kmod as an argument 
    This is a 200 meter box case '''

P = 0.02          # energy level (estimated so that kpHs is reasonable)
L = 200           # domain size
kp = 2*np.pi/40   # peak wavenumber
N_mode = 32       # number of nodes
N_power = 5       # directional spreading coeff

def shape (kmod):
    ''' Choose values here '''
    global P, kp
    F_kmod = spectrum_PM (P=P, kp=kp, kmod=kmod)
    return F_kmod

kmod, F_kmod, kx, ky, F_kxky_tile = spectrum_gen_linear(shape, N_mode=N_mode, L=L, N_power=N_power)




''' Generate a grid in x-y to visualize random eta '''
N_grid = 1024; L = 200
x = np.linspace(-L/2,L/2,N_grid); y = np.linspace(-L/2,L/2,N_grid)
x_tile, y_tile = np.meshgrid(x, y)
kx_tile, ky_tile = np.meshgrid(kx,ky)
t = 0
eta_tile, phase_tile = eta_random(t, kx_tile, ky_tile, F_kxky_tile, x_tile, y_tile)
print('kpHs = %g' %(kp*np.std(eta_tile)*4))



''' Visualization '''
# plt.imshow(eta_tile, cmap='RdBu_r', vmax=eta_tile.max(), vmin=-eta_tile.max())

fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
s=ax.pcolormesh(x,y,eta_tile*kp, cmap='RdBu_r', vmax=eta_tile.max(), vmin=-eta_tile.max())
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_title(r'$\eta k_p$')
plt.colorbar(s,ax=ax)
fig.savefig('eta_ini.png')

deta = xr.DataArray(eta_tile, coords=[x, y], dims=["x", "y"])

window = None #'hann'
nfactor = 4
truncate = True
detrend = None
window_correction = False #True
F_eta =  xrft.isotropic_power_spectrum(deta, dim=('x','y'), window=window, nfactor=nfactor, 
                                truncate=truncate, detrend=detrend, window_correction=window_correction)

fig, ax = plt.subplots(1,1,figsize = (5,5),constrained_layout=True,dpi=100)
ax.loglog(kmod*L, F_kmod*kp**3, c='k', label='PM')
ax.loglog(F_eta.freq_r*np.pi*2*L, F_eta*kp**3, c='g', label=r'computed from $\eta$')
ax.set_ylim([1e-7,1e-2])
ax.set_xlim([10,1000])
ax.set_xlabel('k*L')
ax.set_ylabel(r'$\phi(k)*kp^3$')
plt.legend()
fig.savefig('initial_spectrum_k.png')

kline = np.linspace(100,300)
ax.plot(kline,kline**(-3)*2000, c='grey')
ax.text(180, 6.3e-4, '$k^{-3}$', color='gray')

fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
s=ax.pcolormesh(kx_tile, ky_tile, F_kxky_tile*kp**3, cmap='viridis')
ax.set_xlabel('kx')
ax.set_ylabel('ky')
plt.colorbar(s,ax=ax)
fig.savefig('initial_spectra_kxky.png')

plt.show()



''' When happy, output the spectra '''
print('Spectrum array shape:', F_kxky_tile.shape)
path = "./" # path of choice
file = "F_kxky" #_P%g_200m" %(P)
print('Written to', path+file)
fF = open(path + file, 'bw')
F_output = F_kxky_tile.astype('float32'); F_output.tofile(fF)
# output kx
file = "kx" #_200m"
print('Written to', path+file)
fx = open(path + file, 'bw')
F_output = kx.astype('float32'); F_output.tofile(fx)
# output ky
file = "ky" #_200m"
print('Written to', path+file)
fy = open(path + file, 'bw')
F_output = ky.astype('float32'); F_output.tofile(fy)

