"""

first run the two pieces of code:
# python
python3 gen_ini.py

# C
make
./main


Then execute this code with
python3 compare.py

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import xarray as xr
import xrft

import os.path
import sys
# add libpy
dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, '../libpy/')
sys.path.append( filename )
from fftlib import get_spec_1D

L=200.
kp=10*np.pi/L
N_grid = 1024
x = np.linspace(-L/2,L/2,N_grid)
y = x

# for spectra compute with xrft
window = None #'hann'
nfactor = 4
truncate = True
detrend = None
window_correction = False #True

# files
f_C = "F_kxky_C"
f_py = "F_kxky_py"
fkx_C, fky_C = "kx_C", "ky_C"
fkx_py, fky_py = "kx_py", "ky_py"


# data from Jiarong's plot
Jpath = "./data_Jiarong/"
Js_0 = np.loadtxt(Jpath+'2023_fig3c_0.txt',skiprows=1,delimiter=",")
Js_20 = np.loadtxt(Jpath+'2023_fig3c_25.txt',skiprows=1,delimiter=",")
Js_100 = np.loadtxt(Jpath+'2023_fig3c_124.txt',skiprows=1,delimiter=",")
Js_120 = np.loadtxt(Jpath+'2023_fig3c_149.txt',skiprows=1,delimiter=",")
Js = [Js_0, Js_20, Js_100, Js_120]

# data from C code
kx_C = np.fromfile(fkx_C)
ky_C = np.fromfile(fky_C)
rawC = np.fromfile(f_C)
Fkmod_C = np.fromfile("F_k_C")
kmod_C = np.fromfile("kmod_C")
dataC = np.zeros((len(kx_C),len(ky_C)))
raw_eta_C = np.fromfile("eta_C")
eta_C = np.zeros((N_grid,N_grid))
for row in range(len(kx_C)):
    for col in range(len(ky_C)):
        index = row*len(kx_C) + col
        dataC[row, col] = rawC[index]
for i in range(N_grid):
    for j in range(N_grid):
        index = i*N_grid + j
        eta_C[i,j] = raw_eta_C[index]
# deta_C = xr.DataArray(eta_C, coords=[x, y], dims=["x", "y"])
# F_eta_C =  xrft.isotropic_power_spectrum(deta_C, dim=('x','y'), window=window, nfactor=nfactor, 
#                                 truncate=truncate, detrend=detrend, window_correction=window_correction)
kr_C, F_eta_C = get_spec_1D(eta_C, eta_C, L/N_grid)


# data from python code
kx_py = np.fromfile(fkx_py)
ky_py = np.fromfile(fky_py)
rawpy = np.fromfile(f_py)
Fkmod_py = np.fromfile("F_k_py")
kmod_py = np.fromfile("kmod_py")
datapy = np.zeros((len(kx_py),len(ky_py)))
raw_eta_py = np.fromfile("eta_py")
eta_py = np.zeros((N_grid,N_grid))
for row in range(len(kx_py)):
    for col in range(len(ky_py)):
        index = row*len(kx_py) + col
        #index = col*len(kx_py) + row
        datapy[row, col] = rawpy[index]
for i in range(N_grid):
    for j in range(N_grid):
        index = i*N_grid + j
        eta_py[i,j] = raw_eta_py[index]
# deta_py = xr.DataArray(eta_py, coords=[x, y], dims=["x", "y"])
# F_eta_py =  xrft.isotropic_power_spectrum(deta_py, dim=('x','y'), window=window, nfactor=nfactor, 
#                                 truncate=truncate, detrend=detrend, window_correction=window_correction)

kr_py, F_eta_py = get_spec_1D(eta_py, eta_py, L/N_grid)

# compare kx ky
if False:
    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
    ax.plot(kx_C, label='kx_C')
    ax.plot(kx_py, label='kx_py', ls='--')
    ax.plot(ky_C, label='ky_C')
    ax.plot(ky_py, label='ky_py', ls='--')
    plt.legend()
# -> they are the same ! 
# so now I use the one from the py code

# Compare F_k
if True:
    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
    ax.loglog(kmod_py*L, Fkmod_py*kp**3, label="PM (py)") 
    ax.loglog(kmod_C*L, Fkmod_C*kp**3, label="PM (C)", ls='--') 
    # ax.loglog(F_eta_C.freq_r*np.pi*2*L, F_eta_C*kp**3, c='r', label=r'computed from $\eta$ (C)')
    # ax.loglog(F_eta_py.freq_r*np.pi*2*L, F_eta_py*kp**3, c='b', label=r'computed from $\eta$ (py)', ls='--')
    ax.loglog(kr_C*np.pi*2*L, (F_eta_C/(2*np.pi))*kp**3, c='r', label=r'computed from $\eta$ (C)')
    ax.loglog(kr_py*np.pi*2*L, (F_eta_py/(2*np.pi))*kp**3, c='b', label=r'computed from $\eta$ (py)', ls='--')
    ax.loglog(Js[0][:,0],Js[0][:,1], c='k', ls='--', label="Jiarong's paper")
    ax.set_xlim([10,1000])
    ax.set_ylim([1e-7,1e-2])
    ax.set_ylabel(r'$F(k)_C.k_p^3$')
    ax.set_xlabel(r'$k_p L$')
    ax.vlines(kp*L,1e-8,1,ls='--', colors='gray') # kp
    plt.legend()
    fig.savefig('comparison_Fk.svg')

# compare F_kxky
if False:
    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
    s=ax.pcolormesh(kx_C,ky_C,dataC.T, shading='nearest', vmin=0,vmax=3)
    ax.set_xlabel('kx')
    ax.set_ylabel('ky')
    ax.set_title('C')
    plt.colorbar(s,ax=ax)
    fig.savefig('F_kxky_C.png')

    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
    s=ax.pcolormesh(kx_py,ky_py,datapy.T, shading='nearest', vmin=0,vmax=3)
    ax.set_xlabel('kx')
    ax.set_ylabel('ky')
    ax.set_title('C')
    plt.colorbar(s,ax=ax)
    fig.savefig('F_kxky_py.png')

    # relative error
    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
    s=ax.pcolormesh(kx_py,ky_py,(dataC.T-datapy.T)/np.amax(datapy.T)*100, shading='nearest', 
                    cmap='bwr', vmin=-1, vmax=1)
    ax.set_xlabel('kx')
    ax.set_ylabel('ky')
    ax.set_title('(C-py)/py.amax (%)')
    plt.colorbar(s,ax=ax)
    fig.savefig('diff_F_kxky.svg')

# compare eta
if True:

    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
    s=ax.pcolormesh(x,y,eta_C.T*kp,cmap='Greys_r',vmin = -0.15, vmax=0.3)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('eta_C * kp')
    plt.colorbar(s,ax=ax)
    fig.savefig('eta_C.png')

    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
    s=ax.pcolormesh(x,y,eta_py*kp,cmap='Greys_r',vmin = -0.15, vmax=0.3)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('eta_py * kp')
    plt.colorbar(s,ax=ax)
    fig.savefig('eta_py.png')

    
    


# s = xrft.isotropic_power_spectrum(eta, 
#                         dim=('x','y'), 
#                         window=window,
#                         nfactor=nfactor, 
#                         truncate=truncate, 
#                         detrend=detrend,
#                         window_correction=window_correction)


plt.show()
