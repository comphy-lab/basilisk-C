import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

path="lee/out.nc"
ds = xr.open_dataset(path)

# adding z to variables
ds["z"] = ds.h*0. 

h=ds.h[:,0,:,:]
ds['z'][:,0,:,:] = ds.zb.values + h/2
for l in np.arange(1,len(ds.level)):    
    ds['z'][:,l,:,:] = ds['z'][:,l-1,:,:] + h/2 
    h=ds.h[:,l,:,:]
    ds['z'][:,l,:,:] += h/2
        

# select only right part of domaine
ds=ds.sel(x=slice(1500,2000))
ds=ds.isel(y=-1)

# avg over layer ~ avg at z=cst
ds=ds.mean(dim='x')

timelist=[0,-1]
dpi=100
fig, ax = plt.subplots(1,1,figsize = (3,3),constrained_layout=True,dpi=dpi)
for itime in timelist:
    Z = ds.z.isel(time=itime)
    ax.plot(ds.isel(time=itime).T,Z, label=f"time=%.2fh" % (ds.time[itime].values/3600,))
ax.set_xlabel('T (°C)')
ax.set_ylabel('depth (m)')
ax.legend()
plt.show()
