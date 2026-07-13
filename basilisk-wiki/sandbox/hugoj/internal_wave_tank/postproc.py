import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

# add libpy
import os.path
import sys
dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, '../libpy/')
sys.path.append( filename )
from data_reader import read_data


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
        

ds=ds.isel(y=-1)
# select only right part of domaine
ds2=ds.sel(x=slice(1500,2000))

# avg over layer ~ avg at z=cst
ds2=ds2.mean(dim='x')

timelist=[0,-1]
dpi=100
fig, ax = plt.subplots(1,1,figsize = (3,3),constrained_layout=True,dpi=dpi)
for itime in timelist:
    Z = ds2.z.isel(time=itime)
    ax.plot(ds2.isel(time=itime).T,Z, label=f"time=%.2fh" % (ds2.time[itime].values/3600,))
ax.set_xlabel('T (°C)')
ax.set_ylabel('depth (m)')
ax.legend()
#plt.show()



# velocity field
jumpx=8
jumpz=2
ds = ds.isel(time=-1)
print(ds)
X,LEVEL = np.meshgrid(ds.x.values,ds.level.values)
Z = ds.z.values
C = ds.T.values
U = ds['u.x'].values
V = ds['w'].values
fig, ax = plt.subplots(1,1,figsize = (9,3),constrained_layout=True,dpi=dpi)
s=ax.pcolormesh(X,Z,C,cmap='jet')
plt.colorbar(s,ax=ax)
ax.quiver(X[::jumpz,::jumpx], 
          Z[::jumpz,::jumpx], 
          U[::jumpz,::jumpx],
          V[::jumpz,::jumpx],
          angles='xy',pivot='mid',headlength=4.5,scale=10)
ax.set_ylabel('depth (m)')
ax.set_xlabel('x (m)')
plt.show()



# plt.show()
