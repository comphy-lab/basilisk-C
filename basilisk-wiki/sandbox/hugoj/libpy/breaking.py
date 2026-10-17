import numpy as np
from skimage.feature import hessian_matrix, hessian_matrix_eigvals
from scipy import ndimage as ndi
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

"""
Code from https://github.com/jiarong-wu/multilayer_breaking/blob/bab0b535be26eaf3021bf5a34c9bf20d1962e483/postprocessing/mlpython/breaking.py
"""


# The detection function
def detect_ridges(gray, sigma=1.0):
    H_elems = hessian_matrix(
        gray, sigma=sigma, order="rc", use_gaussian_derivatives=False
    )
    maxima_ridges, minima_ridges = hessian_matrix_eigvals(H_elems)
    return maxima_ridges, minima_ridges


# detect slope based on gradient threshold
def detect_slope(grad, threshold):
    return np.logical_not(grad > threshold)


# bins for predefined kps
def get_bins(kp):
    try:
        if kp == 2 * np.pi / 10:
            bins = np.array([0.01])
            bins = np.concatenate((bins, np.arange(0.25, 5, 0.25)))
        elif kp == 2 * np.pi / 40:
            bins = np.array([0.01])
            bins = np.concatenate((bins, np.arange(0.5, 10, 0.25)))
        elif kp == 2 * np.pi / 100:
            bins = np.array([0.01])
            bins = np.concatenate((bins, np.arange(0.5, 15, 0.25)))
    except NameError:
        print("No pre-defined bins for this peak wavenumber case!")
    return bins


# Detect breaking fronts based on geometry of surface elevation (eta) and compute the statistics
#    of the velocity of the detected front
def simple_mapping(
    eta,
    ux,
    uy,
    delta,
    bins=[],
    method=0,
    threshold=0,
    EXTRA_FILTER=True,
    return_mask=False,
):
    """Arguments:
        eta, ux, uy: the 2D fields
        delta: grid spacing
        bins: velocity bins where the crest is counted
        method: 0-ridge detection; 1-slope thresholding
        threshold: if method=0, an empirical (ad hoc) threshold to filter the ridges; if method=1, the slope threshold
    Outputs:
        hist: 1D array of velocity binned crest length (in unit of pixels), averaged over tseries.
    """
    N = eta.shape[0]  # square eta Nx=Ny

    if method == 0:
        _, b = detect_ridges(eta, sigma=1.0)  # Maxima and minima ridges
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
    """ A edge has two sides ! 
        Here consider that the front is at the right edge
    """
    a_ = np.zeros(b_.shape)

    # for i in range(0, N - 1):
    #     for j in range(1, N - 1):
    #         if (b_[i, j - 1] > 0) and (b_[i, j + 1] == 0):
    #             a_[i, j] = 1
    # vectorized version of the loop above. Using .values to work with numpy array
    a_[:-1, 1:-1] = (b_[:-1, :-2].values > 0) & (b_[:-1, 2:].values == 0)

    # mask is when value = 1 (filters give detected = True) so
    # we invert the values
    a_ = np.logical_not(a_)

    # Extra filtering by above 2.5 sigma (we can kind of use something between 2.5-3\sigma)
    # if EXTRA_FILTER == True:
    #     height_filter = np.logical_not(eta < 2.5*np.var(eta)**0.5)
    #     a_ = a_*height_filter
    ux_a = np.ma.masked_where(a_, ux, copy=True)
    uy_a = np.ma.masked_where(a_, uy, copy=True)
    # not in threshold = masked. Need to use .compressed() to avoid taking into account masked values
    # see https://stackoverflow.com/questions/3610040/how-to-create-the-histogram-of-an-array-with-masked-values-in-numpy
    hist, _ = np.histogram(
        ((ux_a**2 + uy_a**2) ** 0.5).compressed(), bins=bins, density=False
    )
    if return_mask:
        return hist, a_
    else:
        return hist


def label_periodic(points, structure):
    """Connected-component labelling, periodic in both axes.
    points: 2D bool array (axis 0 = y, axis 1 = x)
    structure: odd-sized bool array
    Returns (labels as int32, number of objects)."""
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
    """Vote count per (previous label, current label) pair.
    Points are (y, x) in metres, in [0, box)."""
    if len(prev_pts) == 0 or len(cur_pts) == 0:
        return np.empty((0, 2), int), np.empty(0, int)
    tree = cKDTree(prev_pts, boxsize=box)  # periodic distances
    d, i = tree.query(cur_pts, distance_upper_bound=r_max)
    ok = np.isfinite(d)  # no neighbour -> d = inf
    pairs = np.stack([prev_lab[i[ok]], cur_lab[ok]], axis=1)
    return np.unique(pairs, axis=0, return_counts=True)


def btrack(frames, dx, dy, structure, r_max, events=None):
    """frames: iterable of 2D bool arrays (True on points of interest).
    Returns a list of int arrays holding persistent track ids (0 = background).
    If `events` is a list, it is filled with tuples:
      ("birth", t, id)            new object with no parent
      ("split", t, child, parent) new id created from an existing track
      ("merge", t, ended, into)   track `ended` absorbed by track `into`
      ("death", t, id)            track disappears with no successor
    """
    next_id, prev, out = 1, None, []
    for t, mask in enumerate(frames):
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
    return out
