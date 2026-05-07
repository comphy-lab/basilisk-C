import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os.path
import sys

# add libpy
dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, '../libpy/')
sys.path.append( filename )
from fftlib import get_spec_1D, get_spec_2D, azimuthal_integral

L=200.
kp=10*np.pi/L
N_grid = 512
x = np.linspace(-L/2,L/2,N_grid)
y = x

# files
path = "./parceval/"
f_C = "F_kxky_C"
fkx_C, fky_C = "kx_C", "ky_C"


# data from C code
kx_C = np.fromfile(path+fkx_C)
ky_C = np.fromfile(path+fky_C)
rawC = np.fromfile(path+f_C)
F_kxky_tile = np.reshape(rawC, (len(kx_C),len(ky_C)))
Fkmod_C = np.fromfile(path+"F_k_C")
kmod_C = np.fromfile(path+"kmod_C")
#F_kxky_tile = np.zeros((len(kx_C),len(ky_C)))
raw_eta_C = np.fromfile(path+"eta_C")
eta_C = np.reshape(raw_eta_C,(N_grid,N_grid))
# for row in range(len(kx_C)):
#     for col in range(len(ky_C)):
#         index = row*len(kx_C) + col
#         F_kxky_tile[row, col] = rawC[index]
# for i in range(N_grid):
#     for j in range(N_grid):
#         index = i*N_grid + j
#         eta_C[i,j] = raw_eta_C[index]


kr_C, F_eta_C = get_spec_1D(eta_C, eta_C, L/N_grid)
kxx,kyy,Fkxky_eta = get_spec_2D(eta_C,eta_C, L/N_grid)



dkx = kx_C[1]-kx_C[0]
dkmod = kmod_C[1]-kmod_C[0]
dkr = kr_C[1] - kr_C[0]
dkxx = kxx[0,1]-kxx[0,0]
dkyy = kyy[1,0]-kyy[0,0]
"""
Variance check
"""
print("var F_kmod (computed from PM(k)) = %f" %(np.nansum(Fkmod_C*dkmod)))
print("var F_kxky (interpolated from F_ktheta) = %f" %(np.nansum(F_kxky_tile)*dkx**2))
print("var F_k (azimuthal integration from eta) = %f" %(np.nansum(F_eta_C*dkr)))
print("var F_kxky (2D spectrum from eta) = %f" %(np.nansum(Fkxky_eta)*dkxx*dkyy))

"""
Spectrum plots
"""

kr_C = kr_C*2*np.pi
kxx = kxx*2*np.pi
kyy = kyy*2*np.pi

# Compare F_k
if True:
    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
    ax.loglog(kmod_C/kp, Fkmod_C*kp**3, label="PM (C)", ls='--') 
    ax.loglog(kr_C/kp, (F_eta_C/(2*np.pi))*kp**3, c='r', label=r'computed from $\eta$ (C)')
    #ax.set_xlim([10,1000])
    ax.set_ylim([1e-7,1e-2])
    ax.set_ylabel(r'$F(k).k_p^3$')
    ax.set_xlabel(r'$k$ /$k_p$')
    #ax.vlines(kp*L,1e-8,1,ls='--', colors='gray') # kp
    ax.vlines(1,1e-8,1,ls='--', colors='gray') # kp
    plt.legend()
    fig.savefig('comparison_Fk.pdf')

# plot F_kxky
if True:
    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
    s=ax.pcolormesh(kx_C/kp,ky_C/kp,F_kxky_tile.T*kp**3, shading='nearest', 
                    norm=matplotlib.colors.LogNorm(vmin=1e-7, vmax=1e-2))
    ax.set_xlabel(r'$k_x$/$k_p$')
    ax.set_ylabel(r'$k_y$/$k_p$')
    ax.set_xlim([-5,5])
    ax.set_ylim([-5,5])
    plt.colorbar(s,ax=ax,label=r'$F(k_x,k_y) kp^{3}$')
    fig.savefig('F_kxky_C.pdf')

    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
    s=ax.pcolormesh(kxx/kp,kyy/kp,(Fkxky_eta/(4*np.pi)).T*kp**3, shading='nearest', 
                    norm=matplotlib.colors.LogNorm(vmin=1e-7, vmax=1e-2))
    ax.set_xlabel(r'$k_x$/$k_p$')
    ax.set_ylabel(r'$k_y$/$k_p$')
    ax.set_xlim([-5,5])
    ax.set_ylim([-5,5])
    plt.colorbar(s,ax=ax,label=r'$F(k_x,k_y) kp^{3}$')
    fig.savefig('F_kxky_eta.pdf')
    # Why does this plot isnt scaled as the one before ??

"""
We plot also the surface elevation field
"""

# eta
if True:
    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
    s=ax.pcolormesh(x/L,x/L,eta_C.T*kp, shading='nearest', cmap='Greys_r')
    ax.set_xlabel(r'$x$/$L$')
    ax.set_ylabel(r'$y$/$L$')
    plt.colorbar(s,ax=ax,label=r'$\eta kp$')
    fig.savefig('eta.png')

#plt.show()
