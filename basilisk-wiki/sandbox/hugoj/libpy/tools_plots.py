import numpy as np
from matplotlib import pyplot as plt
import xarray as xr
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable


def plot_surface_3D(ds: xr.DataArray, 
                    cvar: string, 
                    cmin: float,
                    cmax: float, 
                    zmin: float,
                    zmax: float,
                    fig_tuple, # type ?
                    psave: string):
    """
    TODO: 
    """
    fig, ax = fig_tuple

    x = ds.x.values
    y = ds.y.values
    X, Y = np.meshgrid(x, y)
    eta = ds.h.sum(dim='zl').values + ds.zb.values
    var   = ds[cvar].isel(zl=-1).values
    
    # Normalize cvar to [0,1] for the colormap
    #norm    = Normalize(vmin=var.min(), vmax=var.max())
    norm    = Normalize(vmin=cmin, vmax=cmax)
    cmap    = plt.cm.Greys_r
    fcolors = cmap(norm(var))           # RGBA array shaped like eta

    # fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100, 
    #                        subplot_kw={'projection': '3d'})
    surf = ax.plot_surface(X, Y, eta,
                        facecolors=fcolors,           # color from u, not from eta
                        rstride=1, cstride=1,
                        linewidth=0, antialiased=False)
    ax.set_zlim([zmin,zmax])
    ax.set_xlim([x[0],x[-1]])
    ax.set_ylim([y[0],y[-1]])
    # Colorbar (must be built manually when using facecolors)
    mappable = ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array(var)
    fig.colorbar(mappable, ax=ax, shrink=0.5, aspect=10, label='u', pad=0.15)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('η (eta)')
    ax.set_title(fr'Surface $\eta$ colored by {cvar}')
    if psave is not None:
        fig.savefig(psave)
