import numpy as np
import pyvista as pv


def compute_interfaces(zb, h):
    """
    Cell-center elevations + thickness -> interface elevations.

    Inputs
    ----------
    zb : (ny, nx)
        Bottom topography.
    h : (nz, ny, nx)
        Cell thickness.

    Returns
    -------
    zf : (nz+1, ny, nx)
        Interface elevations.
    """
    nz = h.shape[0]
    zf = np.empty((nz + 1,) + h.shape[1:], dtype=h.dtype)
    zf[0] = zb
    for k in range(nz):
        zf[k + 1] = zf[k] + h[k]
    return zf


def interfaces_to_vertices(zf):
    """
    Average interface elevations onto grid vertices.

    Inputs
    ----------
    zf : (nz+1, ny, nx)

    Returns
    -------
    zv : (nz+1, ny+1, nx+1)
    """

    nz1, ny, nx = zf.shape
    zv = np.empty((nz1, ny + 1, nx + 1), dtype=zf.dtype)
    # Interior
    zv[:, 1:-1, 1:-1] = 0.25 * (
        zf[:, :-1, :-1] + zf[:, :-1, 1:] + zf[:, 1:, :-1] + zf[:, 1:, 1:]
    )
    # South edge
    zv[:, 0, 1:-1] = 0.5 * (zf[:, 0, :-1] + zf[:, 0, 1:])
    # North edge
    zv[:, -1, 1:-1] = 0.5 * (zf[:, -1, :-1] + zf[:, -1, 1:])
    # West edge
    zv[:, 1:-1, 0] = 0.5 * (zf[:, :-1, 0] + zf[:, 1:, 0])
    # East edge
    zv[:, 1:-1, -1] = 0.5 * (zf[:, :-1, -1] + zf[:, 1:, -1])
    # Corners
    zv[:, 0, 0] = zf[:, 0, 0]
    zv[:, 0, -1] = zf[:, 0, -1]
    zv[:, -1, 0] = zf[:, -1, 0]
    zv[:, -1, -1] = zf[:, -1, -1]

    return zv


def vertex_coordinates(x):
    dx = x[1] - x[0]
    xv = np.empty(len(x) + 1)
    xv[1:-1] = 0.5 * (x[:-1] + x[1:])
    xv[0] = x[0] - dx / 2
    xv[-1] = x[-1] + dx / 2
    return xv


def build_3Dgrid(ds):
    # Build vertices
    xv = vertex_coordinates(ds.x.values)  # 1D
    yv = vertex_coordinates(ds.y.values)  # 1D
    zf = compute_interfaces(ds.zb.values, ds.h.values)  # 3D
    zv = interfaces_to_vertices(zf)  # 3D
    # Expand vertices to 3D
    X, Y = np.meshgrid(xv, yv, indexing="xy")  # 2D
    X = np.broadcast_to(X, zv.shape)  # 3D
    Y = np.broadcast_to(Y, zv.shape)  # 3D
    Z = zv
    # Pyvista grid object is X,Y,Z but my data is Z,Y,X
    # so I need to transpose
    X = np.broadcast_to(X, zv.shape).transpose(2, 1, 0)
    Y = np.broadcast_to(Y, zv.shape).transpose(2, 1, 0)
    Z = zv.transpose(2, 1, 0)
    # Build the grid object
    grid = pv.StructuredGrid(X, Y, Z)
    return grid


def surface_extractor(grid, size, which="top", tuple=None):
    nx, ny, nz = size  # Pyvista convention
    if which == "top":
        surface = grid.extract_subset((0, nx, 0, ny, nz, nz))
    elif which == "bottom":
        surface = grid.extract_subset((0, nx, 0, ny, 0, 0))
    elif which == "xmin":
        surface = grid.extract_subset((0, 0, 0, ny, 0, nz))
    elif which == "xmax":
        surface = grid.extract_subset((nx, nx, 0, ny, 0, nz))
    elif which == "ymin":
        surface = grid.extract_subset((0, nx, 0, 0, 0, nz))
    elif which == "ymax":
        surface = grid.extract_subset((0, nx, ny, ny, 0, nz))
    else:
        surface = grid.extract_subset(tuple)
    return surface


def surfaces_domain(grid, size):
    top = surface_extractor(grid, size, which="top")
    xmin = surface_extractor(grid, size, which="xmin")
    xmax = surface_extractor(grid, size, which="xmax")
    ymin = surface_extractor(grid, size, which="ymin")
    ymax = surface_extractor(grid, size, which="ymax")
    return top, xmin, xmax, ymin, ymax


def add_mesh_domain_side(grid, sizes, plotter, vartop, varside):
    nx, ny, nz = sizes
    top, xmin, xmax, ymin, ymax = surfaces_domain(grid, (nx, ny, nz))
    plotter.add_mesh(
        top,
        name="top_surface",
        scalars=vartop["name"],
        clim=vartop["clim"],
        cmap=vartop["cmap"],
        scalar_bar_args={
            "title": vartop["name"],
            "fmt": "%.3f",
            "position_x": 0.25,
            "width": 0.5,
        },
    )
    plotter.add_mesh(
        xmin,
        name="xmin_surface",
        scalars=varside["name"],
        clim=varside["clim"],
        cmap=varside["cmap"],
        scalar_bar_args={
            "title": varside["name"],
            "fmt": "%.3f",
            "position_x": 0.25,
            "width": 0.5,
        },
    )
    plotter.add_mesh(
        xmax,
        name="xmax_surface",
        scalars=varside["name"],
        clim=varside["clim"],
        cmap=varside["cmap"],
    )
    plotter.add_mesh(
        ymax,
        name="ymax_surface",
        scalars=varside["name"],
        clim=varside["clim"],
        cmap=varside["cmap"],
    )
    plotter.add_mesh(
        ymin,
        name="ymin_surface",
        scalars=varside["name"],
        clim=varside["clim"],
        cmap=varside["cmap"],
    )


def add_ticks(plotter, L0, H0, tick_len_min, tick_len_maj, label_offset):
    xmin = -L0 / 2
    xmax = L0 / 2
    ymin = xmin
    ymax = xmax
    zmin = -H0
    zmax = 0

    plotter.show_bounds(
        bounds=(xmin, xmax, ymin, ymax, zmin, zmax),
        xtitle="X (m)",
        ytitle="Y (m)",
        ztitle="Z (m)",
        location="outer",
        show_xlabels=False,
        show_ylabels=False,
        show_zlabels=False,
    )

    # plotter.show_axes()

    # x-axis
    for x in np.arange(xmin, xmax + 1, 10):
        p0 = [x, ymin, zmin]
        p1 = [x, ymin - tick_len_min, zmin]
        p2 = [x, ymin - tick_len_maj, zmin]
        p_label = [x, ymin - tick_len_maj - label_offset, zmin]
        if x % 20 == 0:
            plotter.add_lines(np.array([p0, p2]), color="black", width=1)
        else:
            plotter.add_lines(np.array([p0, p1]), color="black", width=1)
        if x % 20 == 0:
            plotter.add_point_labels(
                [p_label],
                [str(x)],
                show_points=False,
                shape=None,
                font_size=12,
            )
    # y-axis
    for y in np.arange(ymin, ymax + 1, 10):
        p0 = [xmin, y, zmin]
        p1 = [xmin - tick_len_min, y, zmin]
        p2 = [xmin - tick_len_maj, y, zmin]
        p_label = [xmin - tick_len_maj - label_offset, y, zmin]
        if y % 20 == 0:
            plotter.add_lines(np.array([p0, p2]), color="black", width=1)
        else:
            plotter.add_lines(np.array([p0, p1]), color="black", width=1)
        if y % 20 == 0:
            plotter.add_point_labels(
                [p_label],
                [str(y)],
                show_points=False,
                shape=None,
                font_size=12,
            )

    # z-axis
    for z in np.arange(zmin, zmax + 1, 10):
        p0 = [xmin, ymax, z]
        p1 = [xmin - tick_len_min, ymax, z]
        p2 = [xmin - tick_len_maj, ymax, z]
        p_label = [xmin - tick_len_maj - label_offset, ymax, z]
        if z % 20 == 0:
            plotter.add_lines(np.array([p0, p2]), color="black", width=1)
        else:
            plotter.add_lines(np.array([p0, p1]), color="black", width=1)

        if z % 20 == 0:
            plotter.add_point_labels(
                [p_label],
                [str(z)],
                show_points=False,
                shape=None,
                font_size=12,
                always_visible=True,
            )
