"""
# Diagnostics for multilayer simulation
"""

import numpy as np
import gc
import xarray as xr
from scipy.interpolate import interp1d

def grad_velocities(ds, grid):
    '''
    Computes the gradient of velocities from a velocity field named u.x u.y u.z
    
    INPUTS:
        ds: xr.Dataset
        grid: xgcm grid
    OUTPUTS:
        ds: xr.Dataset, the updated dataset
        update: bool, whether the dataset has been modified or not
    '''
    if 'dudx' not in ds.keys():
        update=True
    else:
        upate=False

    if update:
        delta = ds.x[1]-ds.x[0]
        dudx = grid.interp(grid.diff(ds['u.x'], 'X'), 'X')/delta
        dudy = grid.interp(grid.diff(ds['u.x'], 'Y'), 'Y')/delta
        dudzl = grid.interp(grid.diff(ds['u.x'], 'Z'), 'Z')
        dvdx = grid.interp(grid.diff(ds['u.y'], 'X'), 'X')/delta
        dvdy = grid.interp(grid.diff(ds['u.y'], 'Y'), 'Y')/delta
        dvdzl = grid.interp(grid.diff(ds['u.y'], 'Z'), 'Z')
        dwdx = grid.interp(grid.diff(ds['u.z'], 'X'), 'X')/delta
        dwdy = grid.interp(grid.diff(ds['u.z'], 'Y'), 'Y')/delta
        dwdzl = grid.interp(grid.diff(ds['u.z'], 'Z'), 'Z')
        
        dzdx = grid.interp(grid.diff(ds.z, 'X'), 'X')/delta
        dzdy = grid.interp(grid.diff(ds.z, 'Y'), 'Y')/delta
        dzdzl = grid.interp(grid.diff(ds.z, 'Z'), 'Z')
        #print('Fields interpolated!')

        ds['dudz'] = (dudzl/dzdzl).compute()
        ds['dudy'] = (dudy - ds['dudz']*dzdy).compute()
        ds['dudx'] = (dudx - ds['dudz']*dzdx).compute()
        ds['dvdz'] = (dvdzl/dzdzl).compute()
        ds['dvdy'] = (dvdy - ds['dvdz']*dzdy).compute()
        ds['dvdx'] = (dvdx - ds['dvdz']*dzdx).compute()
        ds['dwdz'] = (dwdzl/dzdzl).compute()
        ds['dwdy'] = (dwdy - ds['dwdz']*dzdy).compute()
        ds['dwdx'] = (dwdx - ds['dwdz']*dzdx).compute()
        
        # print('Fields gradient computed!')
        
        del(dudx, dudy, dudzl, dvdx, dvdy, dvdzl, dwdx, dwdy, dwdzl, dzdx, dzdy, dzdzl)
        gc.collect()

    return ds, update

def vorticity(ds, grid):   
    '''
    Computes vorticity (squared) from a velocity field named u.x u.y u.z and a xgcm grid

    INPUTS:
        ds: xr.Dataset
        grid: xgcm grid
    OUTPUTS:
        ds: xr.Dataset, the updated dataset
        update: bool, whether the dataset has been modified or not
    '''
    if 'dudx' not in ds.keys():
        raise Exception('Compute gradient first !')

    if 'omegaxp' not in ds.keys():
        update = True
    else:
        update = False

    if update:
        ds['omegaxp'] = ds.dwdy - ds.dvdz
        ds['omegayp'] = ds.dudz - ds.dwdx
        ds['omegazp'] = ds.dvdx - ds.dudy
        ds['enstrophy'] = 1/2 * (ds.omegaxp**2 + ds.omegayp**2 + ds.omegazp**2)
    return ds, update

def dissipation(ds, grid):
    '''
    Computes dissipation from a velocity field named u.x u.y u.z and a xgcm grid

    INPUTS:
        ds: xr.Dataset
        grid: xgcm grid
    OUTPUTS:
        ds: xr.Dataset, the updated dataset
        update: bool, whether the dataset has been modified or not
    '''
    if 'dudx' not in ds.keys():
        raise Exception('Compute gradient first !')
    if 'epsilon' not in ds.keys():
        update = True
    else:
        update = False

    if update:
        ds['epsilon'] = 2*(ds.dudx**2 + 
            2*((ds.dudy + ds.dvdx)/2.)**2 + 
            2*((ds.dudz+ds.dwdx)/2.)**2 + 
            ds.dvdy**2 + 
            2*((ds.dvdz+ds.dwdy)/2.)**2 + 
            ds.dwdz**2)
    return ds, update



# ufunc for z interpolation
# (from Jiarong)
# class interp1d_class():
#     def __init__(self, znew, fill_value):
#         self.znew = znew
#         self.fill_value = fill_value
#     def __call__(self, x, y):
#         f = interp1d(x, y, kind='linear', fill_value=self.fill_value, bounds_error=False)
#         return f(self.znew)
#
#
# ''' fill_value: array-like or (array-like, array_like) or 'extrapolate' '''
# interpz = lambda znew, z, ds, fill_value: xr.apply_ufunc (
#     interp1d_class(znew, fill_value),
#     #lambda x, y: interpolate_1d(x, y, znew),  # The function to apply if interpolate_1d is not defined as a class
#     z,  # Original vertical coordinates
#     ds,  # Data (Dataarray or Dataset) to interpolate
#     input_core_dims=[['zl'], ['zl']],  # Core dimensions for each input
#     output_core_dims=[['zl']],  # Core dimensions for the output
#     exclude_dims=set(('zl',)),
#     dask_gufunc_kwargs={'output_sizes':{'zl':len(znew)}},
#     output_dtypes=['float32'],
#     vectorize=True,  # Enable vectorization
#     dask="parallelized",  # Parallelize using Dask if the data is large
# )
#

def interp_profile(z_prof, data_prof, znew, fill_value):
    '''
        Interpolate a single profile. z_prof and data_prof are 1D numpy arrays.
    INPUTS:
        z_prof: 
        data_prof:
        znew:
        fill_value:
    OUTPUTS:
        numpy array
    '''
    return np.interp(znew, z_prof, data_prof, left=fill_value, right=fill_value)

def interpz(z, data, znew, fill_value):
    '''
        Wrapper of 'interp_profile' for a // use with xarray/dask.
        Interpolate data, located at points z on points znew.
        The fill_value is used when znew is out the range of z.

        z and data are typically 4D with dimensions such as [t,level,y,x].

    INPUTS:
        z: xr.DataArray, N-dimensions array of the position of data 
        data: xr.DataArray, N-dimensions array of the data to interpolate
        znew: numpy array, 1D array position to interpolate data on
    OUTPUTS:
        xr.DataArray: the data variable interpolated at znew.

    NOTE:
        This function doesnt trigger compilation.
    '''
    return xr.apply_ufunc(
        interp_profile,
        z, data,                                        # both (time, y, x, zl)
        kwargs={'znew': znew, 'fill_value':fill_value},
        input_core_dims=[['zl'], ['zl']],            # operate along zl
        output_core_dims=[['znew']],                 # output has new dim
        exclude_dims=set(('zl',)),                   # zl disappears
        dask='parallelized',
        vectorize=True,                              # loop over time,y,x
        output_dtypes=[float],
        dask_gufunc_kwargs={'output_sizes': {'znew': len(znew)}}
    )
