import xarray as xr
import numpy as np
from xgcm import Grid

# read netcdf, set up for xgcm
def read_bas_data(filename: string, chunks='auto', dtype='float32'):
    """
    This function reads the output of Basilisk (using bderembl/libs/netcdf_pas.h)
        and return a dataset with xgcm coordinates to allow for easy differentiation

    INPUT:
        filename: a string where the .nc is
    OUPUT:
        dataset: a xr.Dataset
    """
    ds = xr.open_mfdataset(filename, chunks=chunks)
    zb = ds.zb[0,0,0].values

    # building z and left side of z
    Nt, Nl, Ny, Nx = ds.w.shape
    z = np.zeros((Nt, Nl, Ny, Nx))
    z_l = np.zeros((Nt, Nl, Ny, Nx))
    z[:,0,:,:] += ds.h[:,0,:,:].values/2
    for k in range(1,Nl):
        z[:,k,:,:] += (z[:,k-1,:,:]
                    + ds.h[:,k-1,:,:].values/2
                    + ds.h[:,k,:,:].values/2 )
        z_l[:,k,:,:] = np.sum(ds.h[:,:k].values,axis=1)
    ds['z'] = xr.DataArray(z+zb,dims=['time','level','y','x']).astype(dtype)
    ds['z_l'] = xr.DataArray(z_l+zb,dims=['time','level','y','x']).astype(dtype)
    
    # adding left side coordinates
    #   and z and z_l
    dx = (ds.x[1] - ds.x[0]).values
    dy = (ds.y[1] - ds.y[0]).values
    coords = {
            "x_l":(["x_l"], ds.x.data - dx/2),
            "y_l":(["y_l"], ds.y.data - dy/2),
            #"zl":(["zl"], ds.level.data), # doing nothing but its cleaner
            "zl_l":(['zl_l'], ds.level.data-0.5),
            #"z":(['z'], ds.z.data),
            #"z_l":(['z_l'], ds.z_l.data),
            }
    ds = ds.assign_coords(coords)
    #ds = ds.rename({'u.z':'w'}) # this could be set in bderemble/netcdf_pas.h
    ds = ds.rename({'level':'zl', 'w':'u.z'})

    # xgcm grid
    grid = build_grid(ds)
    return ds, grid

def read_data(filename: string, chunks='auto', dtype='float32'):
    return read_bas_data(filename, chunks, dtype)

def build_grid(ds):
    # xgcm grid
    grid = Grid(ds,
                coords={'X':{'center':'x','left':'x_l'},
                        'Y':{'center':'y','left':'y_l'},
                        'Z':{'center':'zl','left':'zl_l'}},
                autoparse_metadata=False,
                # periodic={'X':'True','Y':'True','Z':'False'},
                boundary={'X':'periodic','Y':'periodic', 'Z':'fill'},
                fill_value={'Z':0})
    return grid


