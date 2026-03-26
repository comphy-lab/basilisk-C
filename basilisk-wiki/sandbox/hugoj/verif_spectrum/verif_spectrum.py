"""
# Azimuthal integration of 2D spectrum

The goals of this script are:

* compare several methods for computing azimuthal average of 2D spectra
* make sure that every methods can retrieve the variance of the signal
* use 3 tests cases:
      1) start from F(k) = 1
      2) start from F(k) = PM(k)
      3) start from a synthetic eta, generated from F(k,theta) = PM(k).Psi(theta)

You can go from 1) to 2) by changing the function 'spectrum_PM'

Reminder:

Parceval equalities:

$$
\begin{aligned}
	<\eta^{2}> = \int_{0}^{\infty} \int_{-\pi}^{\pi} E_1(k,\theta) k dk d\theta =
	\int_{0}^{\infty} \int_{0}^{\infty} E_2 (k_x,k_y) d k_x d k_y
\end{aligned}
$$

with $E_1(k,\theta) = E_2 (k_x,k_y) k $

The omnidirectionnal wavenumber spectrum is:

$$
\begin{aligned}
	\phi(k) = \int_{- \pi}^{\pi} E(k,\theta) k d \theta
\end{aligned}
$$

Notes:

- PM is Pierson-Moscowitz spectrum (sea state)
- geostrokit azimuthal integration is way faster than interpolation then integration
- interpolation introduces errors.
- intergration with the 'azimuthal_integration' is quite precise (error < interp error)
"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os.path
import sys
import time
import xrft
from scipy.fft import fft2, fftfreq
from scipy.interpolate import griddata

# add libpy
dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, '../libpy/')
sys.path.append( filename )
from fftlib import *

"""

## Some parameters

"""
dpi=200 # for the figures

g = 9.81
L0 = 200.0 # m
kp = 10  * np.pi / L0
wp = np.sqrt(g*kp)
fp = wp/(2*np.pi)
Tp = 1/fp
N = 1024
L = 200.
P=0.02

kmaxL = 200*2*np.pi
kminL = 2*np.pi

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)
def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)


def spectrum_case(case): 
    if case==1:
        return 

# analytical spectrum
def spectrum_PM (P, kp, kmod):
    k = np.where(kmod==0.,np.nan,kmod)
    F_kmod = P*k**(-2.5)*np.exp(-1.25*(kp/k)**2)
    return F_kmod

def spectrum_PM_flat(P, kp, kmod, flat=False):
    '''
    Wrapper around spectrum_PM with an optional flat output mode.
    
    Parameters:
        P, kp, kmod : same as spectrum_PM
        flat        : if True, returns an array of ones with the same shape
                      as the spectrum (equivalent to uncommenting *0 + 1)
    '''
    result = spectrum_PM(P, kp, kmod)
    return np.ones_like(result) if flat else result



"""

## PART 1

From an analytical spectrum, recover variances and plot omnidir spec

Note:
Special care should be taken for the evaluation of the function.
The definition of the k and theta arrays are crucial to verify parceval equalities

"""
####################################
print("\n=====================================")
print("PART 1")
print(">From an analytical spectrum, recover variances and plot omnidir spec")
'''

'''

for k,flat in enumerate([True, False]):
    case=k+1
    print(f"Working on case number {case}")
    dk_fine = (kmaxL/L - kminL/L)/200
    fine_k = np.arange(kminL/L, kmaxL/L, dk_fine)
    fine_F = spectrum_PM_flat(P,kp,fine_k,flat=flat)
    dk_fine = fine_k[1]-fine_k[0]
    print("var fine: %f" %(np.sum(fine_F)*dk_fine*2*np.pi))

    # 1) giving it no particular direction
    N_mode = 512

    # we build F(k,theta)
    dtheta = 2*np.pi/N_mode
    dr = (kmaxL/L - kminL/L)/N_mode
    radii = np.arange(kminL/L, kmaxL/L, dr)
    thetas = np.arange(-np.pi, np.pi+dtheta, dtheta)
    r_tile, theta_tile = np.meshgrid(radii, thetas)
    F_ktheta = spectrum_PM_flat(P, kp, r_tile, flat=flat)/r_tile
    F_k = np.sum((F_ktheta*r_tile)[:-1,:], axis=0)*dtheta # why [:-1,:] ? -> do not sum twice theta=pi
    print('var F_k = %f' %(np.nansum(F_k)*dr))

    # plot F_ktheta
    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=dpi,subplot_kw={'projection': 'polar'})
    s=ax.pcolormesh(thetas, radii, (F_ktheta).T*kp**3, vmin=0., vmax=0.003)
    ax.set_title('F_ktheta')
    plt.colorbar(s,ax=ax)

    # we build F(kx,ky) from the F(k,theta) with interpolation
    dkx = (kmaxL/L - kminL/L)/(N_mode)
    kx = np.arange(-kmaxL/L/np.sqrt(2), kmaxL/L/np.sqrt(2), dkx)
    ky = kx
    dkx = kx[1]-kx[0]
    kx_tile, ky_tile = np.meshgrid(kx,ky)
    rc_tile, tc_tile = cart2pol(kx_tile,ky_tile) # kx,ky points in polar
    J = r_tile # E(k,theta) = E(kx,ky)*J avec J = k
    F_kxky_tile = griddata((r_tile.ravel(), theta_tile.ravel()), (F_ktheta).ravel(),
                            (rc_tile, tc_tile), method='linear', fill_value=0)
    print("var F_kxky (interpolated from F_ktheta) = %f" %(np.nansum(F_kxky_tile)*dkx**2))

    '''
    # you could think that the following is right:

    k_norm = np.sqrt(kx_tile**2 + ky_tile**2)
    F_kxky = spectrum_PM(P, kp, k_norm)

    # But in fact you evaluate at ||k|| but place the point at (kx,ky), which is at a distance
    # from center of sqrt(kx**2+ky**2) > kx
    # and kx is used as the radial wavenumber
    '''


    # plot F_kxky
    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=dpi)
    s=ax.pcolormesh(kx_tile, ky_tile, F_kxky_tile*kp**3, vmin=0., vmax=0.003)
    ax.set_title('F_kxky')
    plt.colorbar(s,ax=ax)

    # We test the azimuth_integral function
    # goes from F_kxky to F_k
    # Note: this comes from https://github.com/bderembl/geostrokit/blob/master/geostrokit/fftlib.py
    #       but I modified it to accept a already built kx,ky (i.e. skip their function 'get_wavenumber')
    def azimuthal_integral2(spec_2D, kx, ky, all_kr=False):
        nd = spec_2D.ndim
        if nd == 2:
            spec_2D = spec_2D[None,...]
        nl,N,naux = spec_2D.shape
        
        k,l = np.meshgrid(kx,kx)
        K = np.sqrt(k**2 + l**2)
        if all_kr== False:
            kr = kx[0,int(N/2)+1:]
        elif all_kr == True:
            kmax = K.max()
            dk = np.abs(kx[2]-kx[1])
            kr = dk*np.arange(0,int(kmax/dk)+2)

        dk = kr[1] - kr[0]
        spec_1D = np.zeros((nl,len(kr)))
        for i in range(kr.size):
            kfilt =  (K>=kr[i] - 0.5*dk) & (K<kr[i] + 0.5*dk)
            #    the azimuthal integral is the average value*2*pi*k
            # but to get the same value of the integral for the 1d spetrum
            # and the 2d spectrum, it is better to just sum the cells*dk
            
            spec_1D[:,i] = (spec_2D[:,kfilt].sum(axis=-1))*dk #/(2*np.pi) #*kr[i]*2*np.pi/Nbin
            # Nbin = kfilt.sum()
            # spec_1D[:,i] = (spec_2D[:,kfilt].sum(axis=-1))*kr[i]/(Nbin)*2*np.pi
        return kr, spec_1D.squeeze()

    F1_kr, F1_k = azimuthal_integral2(F_kxky_tile, kx, ky, all_kr=True)
    dkr = F1_kr[1]-F1_kr[0]
    print("var from azimuthal_integral = %f" %(np.nansum(F1_k)*dkr))



    # Plot all azimuth integrated spectra
    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=dpi)
    ax.plot(fine_k*L, fine_F*kp**3, label='analytical', marker='x')
    ax.plot(F1_kr*L, (F1_k/(2*np.pi))*kp**3, label='azimuth integrated (geostrokit)', marker='s', markerfacecolor='None')
    ax.plot(radii*L, (F_k/(2*np.pi))*kp**3, label='azimuth integrated (gridata)', ls='--')
    ax.set_ylim([0.,0.003])
    ax.set_ylabel(r'$\phi(k).k_p^3$')
    ax.set_xlabel(r'kL')
    plt.legend()
    fig.savefig(f'azi_integ_specs_case{case}.svg')
    fig.savefig(f'azi_integ_specs_case{case}.png')


"""

## PART 2

Computing spectra from synthetic eta field

Now we generate a eta field from a given spectrum shape.
We use the specgen.py script from 
    https://github.com/jiarong-wu/multilayer_breaking/blob/main/specgen/specgen.py

The spectrum shape is the Pierson-Moscowitz spectrum, the directionnal function is
a cos^5 (normalized in way that variance is preserved).

"""

print("\n=====================================")
print("PART 2")
print(">Computing spectra from synthetic eta field")

def spectrum_gen_linear(shape, N_mode=32, N_power=5, L=200):
    
    ''' The function to generate a kx-ky spectrum based on uni-directional spectrum and a 
        cos^N (theta) directional spreading.
        Arguments: 
            shape: the spectrum shape (a function)
            N_mode: # of modes (Right now it's N_mode for kx and N_mode+1 for ky; 
                    has to much what's hard-coded in the spectrum.h headerfile.
            L: physical domain size
            '''
    
    '''# of modes for the uni-directional spectrum 
        (doesn't matter as much because of interpolation anyway) '''
    N_kmod = 64; N_theta = 64 # Uniform grid in kmod and ktheta, can be finer than N_mode 
    thetam = 0 # midline direction
    kmod = np.linspace(2*np.pi/L,1.41*100*2*np.pi/L,N_kmod) # Change 
    theta = np.linspace(-0.5*np.pi, 0.5*np.pi, N_theta) + thetam # Centered around thetam
    kmod_tile, theta_tile = np.meshgrid(kmod,theta)

    '''Pick the spectrum shape '''
    F_kmod = shape (kmod) # includes the spectral shape, peak, and energy level of choice
    D_theta = np.abs(np.cos(theta-thetam)**N_power) 
    dtheta = theta[1]-theta[0]
    D_theta = D_theta/np.trapezoid(D_theta, theta)  # Normalize so the sum equals one
    F_kmod_tile, D_theta_tile = np.meshgrid(F_kmod,D_theta) 
    F_kmodtheta_tile = F_kmod_tile*D_theta_tile/kmod_tile # Notice!! Normalize by k
    
    '''Uniform grid in kx,ky'''
    kx = np.arange(1,N_mode+1)*2*np.pi/L # based on the grid, interval can't go smaller then pi/L
    ky = np.arange(-N_mode/2,N_mode/2+1)*2*np.pi/L
    kx_tile, ky_tile = np.meshgrid(kx,ky)
    kxp_tile, kyp_tile = pol2cart(kmod_tile, theta_tile)
    
    ''' Project from uniform k to uniform kx,ky '''
    F_kxky_tile = griddata((kxp_tile.ravel(), kyp_tile.ravel()), F_kmodtheta_tile.ravel(), 
                           (kx_tile, ky_tile), method='linear', fill_value=0) 
    
    return kmod, F_kmod, kx, ky, F_kxky_tile

""" We define a random phase model to build a sea surface elevation"""
def eta_random(t, kx_tile, ky_tile, F_kxky_tile, x_tile, y_tile):

    ''' Function to generate a random field given kx-ky spectrum and the x-y-array. Only
        works with uniformly spaced kx-ky so far.
    '''    
    np.random.seed(0) 
    phase_tile = np.random.random_sample(kx_tile.shape)*2*np.pi 
    eta_tile = np.zeros(x_tile.shape)
    
    kmod_cart_tile, theta_cart_tile = cart2pol(kx_tile,ky_tile)
    omega_tile = (9.8*kmod_cart_tile)**0.5 # linear gravity wave dispertion relation
    dkx = kx_tile[0,1]-kx_tile[0,0]; dky = ky_tile[1,0]-ky_tile[0,0]
    N_grid = x_tile.shape[0]

    for i1 in range(0, N_grid):
        for i2 in range(0, N_grid):
            ampl = (2*F_kxky_tile*dkx*dky)**0.5
            a = (kx_tile*x_tile[i1,i2]+ky_tile*y_tile[i1,i2])-omega_tile*t+phase_tile
            mode = ampl*(np.cos(a)) # uniform spacing in kx and ky
            eta_tile[i1,i2] = np.sum(mode)    
    return eta_tile, phase_tile

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
N_grid = 1024; L = 200
x = np.linspace(-L/2,L/2,N_grid); y = np.linspace(-L/2,L/2,N_grid)
x_tile, y_tile = np.meshgrid(x, y)
kx_tile, ky_tile = np.meshgrid(kx,ky)
t = 0
eta_tile, phase_tile = eta_random(t, kx_tile, ky_tile, F_kxky_tile, x_tile, y_tile)
print('kpHs = %g' %(kp*np.std(eta_tile)*4))

print('var eta = %f' %(np.sum(eta_tile**2)/N_grid**2))
# different method to compute omnidir spec

# 1) geostrokit
k_geostrokit, F_geostrokit = get_spec_1D(eta_tile, eta_tile, L/N_grid, all_kr= True)
dk_g = k_geostrokit[1]-k_geostrokit[0]

# 2) xrft
deta = xr.DataArray(eta_tile, coords=[x, y], dims=["x", "y"])
window = None #'hann'
nfactor = 4
truncate = True
detrend = None
window_correction = False #True
F_xrft =  xrft.isotropic_power_spectrum(deta, dim=('x','y'), window=window, nfactor=nfactor, 
                                truncate=truncate, detrend=detrend, window_correction=window_correction)
dk_xrft = (F_xrft.freq_r[1]-F_xrft.freq_r[0]).values*2*np.pi
#       with k = F_xrft.freq_r*2*np.pi

# 3) Jiarong's code
def jspectrum_integration(eta, L, N, CHECK=False):
    ''' This function performs azimuthal integration of the 2D spectrum (Notice it's 2D instead of 3D with the frequency as well).
        When in doubt, enable CHECK so that the integration is printed out at each step to make sure that 
        units are consistent and we always recover the variance of the data. '''
    if CHECK: print (np.var(eta))
    spectrum = np.fft.fft2(eta) / (N*N)**0.5 # FFT normalization 
    F = np.absolute(spectrum)**2 / N**2 # Per area normalization
    if CHECK: print (np.sum(F))
    wavenumber = 2*np.pi*np.fft.fftfreq(n=N,d=L/N)
    kx = np.fft.fftshift(wavenumber); ky = kx
    kx_tile, ky_tile = np.meshgrid(kx,ky)
    theta = np.arange(-N/2,N/2)/(N)*2*np.pi
    k = wavenumber[0:int(N/2)]
    dkx = kx[1] - kx[0]; dky = ky[1] - ky[0]
    dk = k[1]-k[0]; dtheta = theta[1]-theta[0]
    F_center = np.fft.fftshift(F)/dkx/dky # Further normalization by independent variables
    k_tile, theta_tile = np.meshgrid(k,theta)
    kxp_tile, kyp_tile = pol2cart(k_tile, theta_tile)
    F_center_polar = griddata((kx_tile.ravel(),ky_tile.ravel()), F_center.ravel(), (kxp_tile, kyp_tile), method='nearest')
    F_center_polar_integrated = np.sum(F_center_polar*k_tile, axis=0)*dtheta # Azimuthal integration
    if CHECK: print (np.sum(F_center_polar_integrated)*dk)
    return F_center_polar_integrated
k_Jiarong = 2*np.pi*np.fft.fftfreq(n=N,d=L0/N)[0:int(N/2)]
dk_j = k_Jiarong[1]-k_Jiarong[0]
F_jiarong = jspectrum_integration(eta_tile, L, N_grid, CHECK=False)


print('variances check')
print("geostrokit: %f" %(np.sum(F_geostrokit)*dk_g))
print("xrft: %f" %(np.sum(F_xrft)*dk_xrft))
print('jiarong (=griddata): %f' %(np.sum(F_jiarong)*dk_j))


plt.show()


"""
# Expected output
```
=====================================
PART 1
>From an analytical spectrum, recover variances and plot omnidir spec
Working on case number 1
var fine: 39.281026
var F_k = 39.281026
var F_kxky (interpolated from F_ktheta) = 31.145410
var from azimuthal_integral = 31.145410
Working on case number 2
var fine: 1.040840
var F_k = 1.040843
var F_kxky (interpolated from F_ktheta) = 1.038498
var from azimuthal_integral = 1.038498

=====================================
PART 2
>Computing spectra from synthetic eta field
kpHs = 0.246694
var eta = 0.154155
variances check
geostrokit: 0.154155
xrft: 0.140847
jiarong (=griddata): 0.153773
```
"""

