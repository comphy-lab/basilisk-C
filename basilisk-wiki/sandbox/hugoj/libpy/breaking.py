import numpy as np
from skimage.feature import hessian_matrix, hessian_matrix_eigvals
from scipy import ndimage as ndi
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from dask.diagnostics import ProgressBar
import xarray as xr
import dask
import os
from collections import Counter

"""
Base code from Jiarong Wu ("simple_mapping")
https://github.com/jiarong-wu/multilayer_breaking/blob/bab0b535be26eaf3021bf5a34c9bf20d1962e483/postprocessing/mlpython/breaking.py


"""


# The detection function
def detect_ridges(gray, sigma=1):
    """
    Detects ridge structures in a grayscale image using the Hessian matrix
    and its eigenvalues.

    INPUTS:
    - gray (numpy.ndarray): 2D grayscale image.
    - sigma (int): Standard deviation of the Gaussian filter used for
      smoothing.

    OUTPUTS:
    - maxima_ridges (numpy.ndarray): Eigenvalue corresponding to the maximum
      curvature, highlighting ridge-like structures.
    - minima_ridges (numpy.ndarray): Eigenvalue corresponding to the minimum
      curvature, highlighting valley-like structures.
    """
    H_elems = hessian_matrix(
        gray, sigma=sigma, order="rc", use_gaussian_derivatives=False
    )
    maxima_ridges, minima_ridges = hessian_matrix_eigvals(H_elems)
    return maxima_ridges, minima_ridges


# detect slope based on gradient threshold
def detect_slope(grad, threshold):
    """
    Detects regions where the gradient is below a given threshold.

    INPUTS:
    - grad (numpy.ndarray): Array containing gradient values.
    - threshold (float): Gradient threshold used to classify slope regions.

    OUTPUTS:
    - numpy.ndarray: Boolean array indicating the regions where the gradient
      is lower than or equal to the threshold.
    """
    return np.logical_not(grad > threshold)


def get_bins(kp):
    """
    Returns predefined bins corresponding to a given peak wavenumber.

    INPUTS:
    - kp (float): Peak wavenumber used to select the predefined bins.

    OUTPUTS:
    - bins (numpy.ndarray): Array containing the predefined bin values.
    """
    bins = None

    try:
        if kp == 2 * np.pi / 10:
            bins = np.concatenate(([0.01], np.arange(0.25, 5, 0.25)))
        elif kp == 2 * np.pi / 40:
            bins = np.concatenate(([0.01], np.arange(0.5, 10, 0.25)))
        elif kp == 2 * np.pi / 100:
            bins = np.concatenate(([0.01], np.arange(0.5, 15, 0.25)))
        else:
            raise ValueError(f"No predefined bins for kp={kp}.")

    except ValueError as e:
        print(e)

    return bins


# Detect breaking fronts based on geometry of surface elevation (eta) and compute the statistics
#    of the velocity of the detected front
def simple_mapping(
    eta,
    ux,
    uy,
    delta,
    bins=np.array([]),
    method=0,
    threshold=0,
    EXTRA_FILTER=True,
):
    """
    Detects wave crests from a 2D surface field and computes the velocity
    distribution along the detected crests.

    INPUTS:
    - eta (numpy.ndarray): 2D surface elevation field.
    - ux (numpy.ndarray): 2D x-component velocity field.
    - uy (numpy.ndarray): 2D y-component velocity field.
    - delta (float): Grid spacing used to normalize spatial derivatives.
    - bins (numpy.ndarray): Velocity bin edges used to compute the histogram.
    - method (int): Crest detection method:
        - 0: Hessian-based ridge detection.
        - 1: Gradient-based slope thresholding.
    - threshold (float): Threshold used to filter detected crests.
        For method=0, it is applied to the normalized Hessian eigenvalue.
        For method=1, it is applied to the x-gradient of the surface field.
    - EXTRA_FILTER (bool): If True and method=0, applies an additional
      elevation filter to remove low-amplitude features.

    OUTPUTS:
    - hist (numpy.ndarray): 1D histogram containing the velocity distribution
      along the detected crests. The values correspond to the number of
      detected crest pixels in each velocity bin.
    - a_ (numpy.ndarray): 2D boolean mask indicating the detected crest
      locations.
    """
    # A padding is added to get the periodic BC: an object exiting from one
    # side should be the same as the one entering in the other side
    PAD = 8
    eta, ux, uy = (np.pad(np.asarray(f), PAD, mode="wrap") for f in (eta, ux, uy))

    if method == 0:
        _, b = detect_ridges(eta, sigma=1)  # Maxima and minima ridges
        b_norm = b / delta**2  # Normalize the curvature by grid size
        b_ = np.less(b_norm, threshold)  # Is this value fixed???
        # Extra filtering by above 2.5 sigma (we can kind of use something between 2.5-3\sigma)
        if EXTRA_FILTER:
            height_filter = np.logical_not(eta < 2.5 * np.var(eta) ** 0.5)
            b_ = b_ * height_filter

    elif method == 1:
        eta_gradx = np.gradient(eta, axis=0) / delta
        b_ = detect_slope(eta_gradx, threshold=threshold)
    else:
        raise NameError(f"method={method} behavior is not defined")

    # A edge has two sides !
    #    Here we consider that the front is at the right edge

    a_ = np.zeros(b_.shape)

    # for i in range(0, N - 1):
    #     for j in range(1, N - 1):
    #         if (b_[i, j - 1] > 0) and (b_[i, j + 1] == 0):
    #             a_[i, j] = 1
    # vectorized version of the loop above.
    a_[:-1, 1:-1] = (np.asarray(b_[:-1, :-2]) > 0) & (np.asarray(b_[:-1, 2:]) == 0)

    # mask is when value = 1 (filters give detected = True) so
    # we invert the values
    a_ = np.logical_not(a_)

    ux_a = np.ma.masked_where(a_, ux, copy=True)
    uy_a = np.ma.masked_where(a_, uy, copy=True)
    # not in threshold = masked. Need to use .compressed() to avoid taking into account masked values
    # see https://stackoverflow.com/questions/3610040/how-to-create-the-histogram-of-an-array-with-masked-values-in-numpy

    a_ = np.asarray(a_)[PAD:-PAD, PAD:-PAD]  # True = hidden, as before
    ux_a = ux_a[PAD:-PAD, PAD:-PAD]
    uy_a = uy_a[PAD:-PAD, PAD:-PAD]

    hist, _ = np.histogram(
        ((ux_a**2 + uy_a**2) ** 0.5).compressed(), bins=bins, density=False
    )

    return hist, a_


def simple_mapping_parallel(
    eta,
    ux,
    uy,
    lvl_surf,
    delta,
    bins,
    threshold,
    output_file="crest_mask.nc",
):
    """
    Computes a breaking crest mask from surface elevation and velocity
    fields using Dask-parallelized processing.

    INPUTS:
    - eta (xarray.DataArray): surface elevation field (4D).
    - ux (xarray.DataArray): X direction velocity field (4D).
    - uy (xarray.DataArray): Y direction velocity field (4D).
    - lvl_surf (int): Index of the vertical level corresponding to the
      surface velocity fields.
    - delta (float): Grid spacing used by the crest detection algorithm.
    - bins (numpy.ndarray): Velocity bin edges used to compute the
      velocity histogram.
    - threshold (float): Threshold used for crest detection.
    - output_file (str, optional): Path of the NetCDF file where the
      resulting mask is saved. If None, the result is not saved.

    OUTPUTS:
    - out (xarray.Dataset): Dataset containing the crest mask as a
      uint8 DataArray. A value of 1 indicates a masked point and a
      value of 0 indicates a detected crest point.
    """
    if output_file is not None and os.path.exists(output_file):
        return xr.open_dataset(output_file)

    with ProgressBar():
        hist, l_mask = xr.apply_ufunc(
            simple_mapping,
            eta,
            ux.isel(level=lvl_surf),
            uy.isel(level=lvl_surf),
            input_core_dims=[["y", "x"]] * 3,
            output_core_dims=[["bin"], ["y", "x"]],
            kwargs=dict(
                delta=delta,
                bins=bins,
                threshold=threshold,
            ),
            vectorize=True,
            dask="parallelized",
            output_dtypes=[np.int64, bool],
            dask_gufunc_kwargs={"output_sizes": {"bin": len(bins) - 1}},
        )

        # Trigger the Dask computation.
        hist, l_mask = dask.compute(hist, l_mask)

        out = xr.Dataset(
            {
                "mask": l_mask.astype("uint8"),
            },
            attrs={
                "description": (
                    "Breaking crest mask. 1 = masked (not a crest), 0 = crest point."
                ),
                "threshold": float(threshold),
            },
        )

        out["mask"].attrs["long_name"] = "crest mask (1 = hidden, 0 = crest)"

        if output_file is not None:
            encoding = {
                "mask": {
                    "zlib": True,
                    "complevel": 4,
                    "dtype": "uint8",
                }
            }
            out.to_netcdf(output_file, encoding=encoding)

    return out


def label_periodic(points, structure):
    """
    Performs connected-component labelling on a 2D boolean array with
    periodic boundary conditions along both axes.

    INPUTS:
    - points (numpy.ndarray): 2D boolean array indicating the points to
      label. Axis 0 corresponds to y and axis 1 to x.
    - structure (numpy.ndarray): Odd-sized 2D boolean array defining the
      connectivity used for connected-component labelling.

    OUTPUTS:
    - labels (numpy.ndarray): 2D int32 array assigning a unique label to
      each connected component. Periodically connected regions receive
      the same label.
    - n_objects (int): Number of connected components, excluding the
      background.
    """
    ny, nx = points.shape
    py, px = structure.shape[0] // 2, structure.shape[1] // 2

    # Pad only by the structure's reach (1 px for a 3x3 structure)
    padded = np.pad(points, ((py, py), (px, px)), mode="wrap")
    lab = np.empty(padded.shape, dtype=np.int32)
    n = ndi.label(padded, structure=structure, output=lab)
    if n == 0:
        return lab[py : py + ny, px : px + nx].copy(), 0

    # Each padded cell is a copy of one original cell: find that cell's label
    iy = (np.arange(ny + 2 * py) - py) % ny + py
    ix = (np.arange(nx + 2 * px) - px) % nx + px
    src = lab[np.ix_(iy, ix)]

    # Copies carrying different labels belong to the same object
    m = lab > 0
    a, b = lab[m], src[m]
    graph = coo_matrix((np.ones(a.size, np.int32), (a, b)), shape=(n + 1, n + 1))
    ncomp, comp = connected_components(graph, directed=False)

    # comp[0] == 0 (background), other objects are numbered 1..ncomp-1
    out = comp.astype(np.int32)[lab[py : py + ny, px : px + nx]]
    return out, ncomp - 1


def match_frames(prev_pts, prev_lab, cur_pts, cur_lab, r_max, box):
    """
    Finds correspondences between objects in two consecutive frames using
    periodic nearest-neighbour distances and counts the matches between
    previous and current object labels.

    INPUTS:
    - prev_pts (numpy.ndarray): 2D array of shape (N, 2) containing the
      (y, x) coordinates of points in the previous frame, in metres.
    - prev_lab (numpy.ndarray): 1D integer array assigning an object label
      to each point in the previous frame.
    - cur_pts (numpy.ndarray): 2D array of shape (M, 2) containing the
      (y, x) coordinates of points in the current frame, in metres.
    - cur_lab (numpy.ndarray): 1D integer array assigning an object label
      to each point in the current frame.
    - r_max (float): Maximum distance in metres for matching points between
      the two frames.
    - box (tuple): Periodic domain size as (box_y, box_x), in metres.

    OUTPUTS:
    - pairs (numpy.ndarray): 2D integer array of shape (K, 2), where each
      row contains a (previous_label, current_label) pair with at least
      one matched point.
    - counts (numpy.ndarray): 1D integer array containing the number of
      matched points for each label pair in ``pairs``.
    """
    if len(prev_pts) == 0 or len(cur_pts) == 0:
        return np.empty((0, 2), int), np.empty(0, int)
    tree = cKDTree(prev_pts, boxsize=box)  # periodic distances
    d, i = tree.query(cur_pts, distance_upper_bound=r_max)
    ok = np.isfinite(d)  # no neighbour -> d = inf
    pairs = np.stack([prev_lab[i[ok]], cur_lab[ok]], axis=1)
    return np.unique(pairs, axis=0, return_counts=True)


def btrack(frames, dx, dy, structure, r_max, out_file=None, events=None, verbose=True):
    """
    Tracks connected objects across a sequence of 2D binary frames using
    periodic boundary conditions and assigns persistent track IDs to them.

    INPUTS:
    - frames (xarray.DataArray): 3D array (time + 2D space) of boolean arrays,
      where True values indicate points belonging to objects of interest.
    - dx (float): Grid spacing along the x-axis.
    - dy (float): Grid spacing along the y-axis.
    - structure (numpy.ndarray): Odd-sized 2D boolean array defining the
      connectivity used to identify connected objects.
    - r_max (float): Maximum physical distance allowed when matching
      objects between consecutive frames.
    - out_file (string): the name of the ouptput file
    - events (list, optional): List to which tracking events are appended.
      Events are represented as tuples:
        - ("birth", t, id): A new object appears with no parent.
        - ("split", t, child, parent): A new object appears from an
          existing tracked object.
        - ("merge", t, ended, into): A tracked object is absorbed into
          another object.
        - ("death", t, id): A tracked object disappears with no successor.
        - verbose (bool): wether to show progressBar and stats of events or not

    OUTPUTS:
    - tracks (xarray.DataArray): 3D int64 array, one timestep per input
      frame, containing the persistent track ID at each pixel. A value of
      0 represents the background.
    """
    next_id, prev, out = 1, None, []
    nframes = frames.shape[0]
    for t in range(nframes):
        mask = np.logical_not(frames[t])
        ny, nx = mask.shape
        box = (ny * dy, nx * dx)
        labels, n = label_periodic(mask, structure)
        ys, xs = np.nonzero(mask)
        pts = np.column_stack([ys * dy, xs * dx])
        lab = labels[ys, xs]

        id_of = np.zeros(n + 1, dtype=np.int64)  # label -> track id
        used, vote_parent, vote_child = set(), {}, {}
        if prev is not None and n > 0:
            pairs, counts = match_frames(prev["pts"], prev["lab"], pts, lab, r_max, box)
            for k in np.argsort(-counts):  # strongest claims first
                pl, cl = pairs[k]
                tid = prev["id_of"][pl]
                vote_parent.setdefault(cl, tid)  # best parent of each new object
                vote_child.setdefault(tid, cl)  # best successor of each old track
                if id_of[cl] == 0 and tid not in used:
                    id_of[cl] = tid
                    used.add(tid)

        born = []
        for cl in range(1, n + 1):
            if id_of[cl] == 0:
                id_of[cl] = next_id
                born.append(cl)
                next_id += 1

        if events is not None:
            for cl in born:
                if cl in vote_parent:
                    events.append(("split", t, int(id_of[cl]), int(vote_parent[cl])))
                else:
                    events.append(("birth", t, int(id_of[cl])))
            if prev is not None:
                for tid in prev["id_of"][1:]:
                    if tid not in used:
                        if tid in vote_child:
                            events.append(
                                ("merge", t, int(tid), int(id_of[vote_child[tid]]))
                            )
                        else:
                            events.append(("death", t, int(tid)))

        out.append(id_of[labels])
        prev = dict(pts=pts, lab=lab, id_of=id_of)

        if verbose:
            pct = 100 * (t + 1) / nframes
            print(f"  timestep {t + 1}/{nframes} ({pct:.0f}%)", end="\r", flush=True)

    if verbose:
        if events is not None:
            # births, splits, merges, deaths
            print("   -> Counter of events: " + str(Counter(e[0] for e in events)))

        # lifetime of each track, in frames
        life = Counter(i for a in out for i in np.unique(a[a > 0]))
        lengths = np.array(list(life.values()))
        print(
            f"   -> {len(life)} tracks, median lifetime {np.median(lengths):.0f} frames, "
            f"{(lengths == 1).mean():.0%} last a single frame"
        )

    # return a xarray.DataArray
    tracks = xr.DataArray(
        np.stack(out, axis=0),
        dims=frames.dims,
        coords=frames.coords,
        name="track_id",
        attrs={
            "long_name": "persistent object track ID",
            "description": (
                "Persistent integer ID assigned to each tracked object. "
                "0 represents the background."
            ),
        },
    )

    if out_file is not None:
        encoding = {
            "track_id": {
                "zlib": True,
                "complevel": 4,
                "dtype": "uint16",
            }
        }
        tracks.to_netcdf(out_file, encoding=encoding)

    return tracks
