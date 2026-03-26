#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
# basic spectral tools for data in square domains
# based on https://github.com/adekunleoajayi/powerspec
#
# from https://github.com/bderembl/geostrokit/blob/master/geostrokit/fftlib.py

def azimuthal_integral(spec_2D, Delta, all_kr=False):
    '''
    Compute the azimuthal integral of a 2D spectra

    Parameters
    ----------

    spec_2D: array [(nz,) ny,nx] 2d spectra
    Delta: float, grid step
    all_kr: bool. if False, only get radial wave number
    in the inner circle, if True, gets all wave numbers. 
    Default is False. You should only use True to check Parseval

    Returns
    -------

    kr: array [nkr]  radial wave number
    spec_1D: array [(nz,) nkr] 1d spectra
    '''
    
    #print('import is ok')

    nd = spec_2D.ndim
    if nd == 2:
        spec_2D = spec_2D[None,...]
    nl,N,naux = spec_2D.shape

    k,l,K,kr = get_wavenumber(N,Delta, all_kr)
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
        
        # the loop is missing the value at K=0:
        # add it with all_kr option
    return kr, spec_1D.squeeze()


def get_wavenumber(N, Delta, all_kr=False):
    '''
    Compute wavenumber and radial wavenumber 
    wave number units are in cycle per unit length

    Parameters
    ----------

    N : int
    Delta: float, grid step
    all_kr: bool. if False, only get radial wave number
        in the inner circle, if True, gets all wave numbers. 
        Default is False. You should only use True to check Parseval

    Returns
    -------

    k: array [ny,nx] zonal wave numbers
    l: array [ny,nx] meridional wave numbers
    K: array [ny,nx]  (k^2 + l^2)
    kr: array [nkr]  radial wave number
    '''

    kx = np.fft.fftshift(np.fft.fftfreq(N,Delta)) # two sided  
    k,l = np.meshgrid(kx,kx)
    K = np.sqrt(k**2 + l**2)
    if all_kr == False:
        kr = k[0,int(N/2)+1:]
    elif all_kr == True:
        kmax = K.max()
        dk = np.abs(kx[2]-kx[1])
        kr = dk*np.arange(0,int(kmax/dk)+2)
    return k,l,K,kr


def get_spec_2D(psi1, psi2, Delta, window=None):
  ''' 
  Compute the 2D cross power spectrum of psi1 and psi2
  If psi1=psi2, then this is simply the power spectrum

   normalization such that parseval is ok: 
   E = 1/V*(np.sum(psi**2)*Delta**2) = np.sum(psi**2)/N**2
     = np.sum(spec_2D)*dk**2
     ~ np.sum(spec_1D)*dk
  (the last equality is not exactly true because the 1d spectra only contains 
  the inner circle to the square k,l. Set all_kr=True in get_spec_1D)
  
  with V = (N*Delta)**2

  Parameters
  ----------

  psi1 : array [(nz,) ny,nx]
  psi2 : array [(nz,) ny,nx]
  Delta: grid step
  window: None or "hanning"

  Returns
  -------

  k: array [ny,nx] zonal wave numbers
  l: array [ny,nx] meridional wave numbers
  spec_2D: array [(nz,) ny,nx] 2d spectra
  '''
  
  N = psi1.shape[-1]
  k,l,K,kr = get_wavenumber(N,Delta)

  if window == 'hanning':
    win1d = np.hanning(N)
    window = win1d[None,:]*win1d[:,None]
  else:
    window = np.ones((N,N))

  psi1 = window*psi1
  psi2 = window*psi2

  psi1_hat = np.fft.fft2(psi1)
  psi2_hat = np.fft.fft2(psi2)
  spec_2D = (psi1_hat*psi2_hat.conj()).real*Delta**2/N**2
  spec_2D = np.fft.fftshift(spec_2D, axes=(-1,-2))
  return k, l, spec_2D


def get_spec_1D(psi1, psi2, Delta, window=None, all_kr= False, averaging='radial'):
  '''
  Compute the 1D power spectrum of the data by azimuthal integral
  of the 2D spectra

  Parameters
  ----------

  psi1 : array [(nz,) ny,nx]
  psi2 : array [(nz,) ny,nx]
  Delta: grid step
  window: None or "hanning"
  all_kr: bool. if False, only get radial wave number
    in the inner circle, if True, gets all wave numbers. 
    Default is False. You should only use True to check Parseval
  averaging: 'radial' (default), or 'xy' 

  Returns
  -------

  if averaging == 'radial'

  kr: array [N] radial wave number
  spec_1D: array [(nz,) N] 1d spectra

  elif averaging == 'xy'

  k: array [N-1] zonal wave numbers
  l: array [N-1] meridional wave numbers
  spec_x: array [(nz,) N-1] 1d spectra in zonal direction
  spec_y: array [(nz,) N-1] 1d spectra in meridional direction

  '''

  k, l, spec_2D = get_spec_2D(psi1, psi2, Delta, window)

  if averaging == 'radial':
    kr, spec_1D = azimuthal_integral(spec_2D, Delta, all_kr)
    return kr, spec_1D

  elif averaging == 'xy':
    si = psi1.shape
    N = si[-1]
    if all_kr == False:
      kx = k[0,int(N/2)+1:]
      ky = l[int(N/2)+1:,0]
    else:
      kx = -k[0,:int(N/2)+1][::-1]
      ky = -l[:int(N/2)+1,0][::-1]
      


    dk = kx[1] - kx[0]
    dl = ky[1] - ky[0]
    spec_x = np.sum(spec_2D,axis=-2)*dl
    spec_y = np.sum(spec_2D,axis=-1)*dk

    if all_kr == False:
      spec_x = spec_x[...,int(N/2)+1:] + spec_x[...,1:int(N/2)][::-1]
      spec_y = spec_y[...,int(N/2)+1:] + spec_y[...,1:int(N/2)][::-1]
    else:
      spec_x[...,1:int(N/2)] += spec_x[...,int(N/2)+1:][...,::-1]
      spec_y[...,1:int(N/2)] += spec_y[...,int(N/2)+1:][...,::-1]
      spec_x = spec_x[...,:int(N/2)+1][...,::-1]
      spec_y = spec_y[...,:int(N/2)+1][...,::-1]


    return kx, ky, spec_x, spec_y




def get_spec_flux(psi1, psi2, Delta, window=None):
  '''
  Compute spectral flux

  Parameters
  ----------

  psi1 : array [(nz,) ny,nx]
  psi2 : array [(nz,) ny,nx]
  Delta: grid step
  window: None or "hanning"

  Returns
  -------

  kr: array [N] radial wave number
  flux: array [(nz,) N] spectral flux

  '''
  k, l, spec_2D = get_spec_2D(psi1, psi2, Delta, window)

  # kr,spec_1D = azimuthal_integral(spec_2D,Delta)
  # flux = -np.cumsum(spec_1D)*dk # integrate from low wavenumbers
  # flux = np.cumsum(spec_1D[::-1])[::-1]*dk # integrate from high wavenumbers

  nd = spec_2D.ndim
  if nd == 2:
    spec_2D = spec_2D[None,...]
  nl,N,naux = spec_2D.shape

  k,l,K,kr = get_wavenumber(N,Delta)
  dk = kr[1] - kr[0]

  flux = np.zeros((nl,len(kr)))
  for i in range(kr.size):
    kfilt =  (kr[i] <= K ) 
    flux[:,i] = (spec_2D[:,kfilt]).sum(axis=-1)*dk*dk
    # next line is probably faster. check it is the same result
    #flux[:,i] = (spec_2D*kfilt).sum(axis=(1,2))*dk*dk

  return kr, flux.squeeze()


def convolve2D(psi1, psi2, norm=False, psi_s=1, kr=0, k_target=0, dk=1, K=0):
  '''
  Convolution of psi1 and psi2 (assume these two fields are in spectral space)

  Parameters
  ----------

  psi1 : array [ny,nx] Shiffted Fourier coefs of psi1
  psi2 : array [ny,nx] Shiffted Fourier coefs of psi1
  norm: Bool if True, normalize the convolution by 1/N**2 (optional)
  psi_s : array [ny,nx] Shiffted Fourier coefs of psi_s (for scalar product)(optional)
  kr: array [N] radial wave number (optional)
  k_target: float (optional)
  dk: float (optional)
  K: wave number magnitude (optional)

  Returns
  -------

  out: array [ny,nx], convolution of psi1*psi2
  kmag: array [ny,nx], magnitude of k (1st dimension is K of psi2 second dimension is K of psi1)

  '''

  # only convol odd arrays
  si = psi1.shape
  N = si[-1]
  Nc = int((N-1)/2)

  if norm:
    cnorm = 1/N**2  # N or N-1??
  else:
    cnorm = 1

  flag_reshape = 0
  if N % 2 == 0:
    flag_reshape = 1
    psi1 = psi1[1:,1:]
    psi2 = psi2[1:,1:]
    if isinstance(K,np.ndarray):
      K = K[1:,1:]
    if isinstance(psi_s,np.ndarray):
      psi_s = psi_s[1:,1:]

  si = psi1.shape
  N = si[-1]
  Nc = int((N-1)/2)

  if not isinstance(psi_s,np.ndarray):
    psi_s = np.ones((N,N))

  out = np.zeros_like(psi1)

  if k_target:
    nk = len(kr)
    kmag = np.zeros((nk,nk))
    kfilt = (K>=k_target - 0.5*dk) & (K<k_target + 0.5*dk)
  else:
    nk = 1
    kmag = np.zeros((nk,nk))
    kfilt = np.zeros((N,N))

  psif = np.flip(np.flip(psi2,-1),-2)

  for j,i in np.ndindex((N,N)):

    # indices of the flipped matrix
    i1 = None if i >= Nc else Nc-i
    i2 = None if i <= Nc else Nc-i
    j1 = None if j >= Nc else Nc-j
    j2 = None if j <= Nc else Nc-j

    # indices of the non-flipped matrix
    k1 = None if i <= Nc else i-Nc
    k2 = None if i >= Nc else i-Nc
    l1 = None if j <= Nc else j-Nc
    l2 = None if j >= Nc else j-Nc
    #print("(",i,",",j,")",i1,i2,j1,j2, "--", k1,k2,l1,l2)
    out[j,i] = psi_s[j,i].conj()*np.sum(psif[j1:j2,i1:i2]*psi1[l1:l2,k1:k2])*cnorm

    if kfilt[j,i]:
      k_ind1 = np.round((K[l1:l2,k1:k2]/dk).flatten()).astype(int)
      k_ind2 = np.round((K[j1:j2,i1:i2]/dk).flatten()).astype(int)
      kmag[k_ind2,k_ind1] += (psi_s[j,i].conj()*
                              (psif[j1:j2,i1:i2]*psi1[l1:l2,k1:k2]*cnorm).flatten()).real

  if flag_reshape:
    out = np.pad(out,((1,0),(1,0)))

  if k_target == 0:
    return out
  else:
    return out, kmag
