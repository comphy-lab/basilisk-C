import xarray as xr
import xrft
import numpy as np
from scipy.interpolate import LinearNDInterpolator, griddata

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)
def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def Spectra2D(eta: xr.DataArray, compute=False):

    window = 'hann'
    nfactor = 4
    truncate = True
    detrend = None
    window_correction = True

    s = xrft.isotropic_power_spectrum(eta, 
                                dim=('x','y'), 
                                window=window,
                                nfactor=nfactor, 
                                truncate=truncate, 
                                detrend=detrend,
                                window_correction=window_correction)

    
    # to do: add checks.
    return s

def jspectrum_integration(eta, L, N, CHECK=False):
            """ This function performs azimuthal integration of the 2D spectrum (Notice it's 2D instead of 3D with the frequency as well).
            When in doubt, enable CHECK so that the integration is printed out at each step to make sure that 
            units are consistent and we always recover the variance of the data. """  
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
