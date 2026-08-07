import os
import numpy as np
from skimage.feature import hessian_matrix, hessian_matrix_eigvals
import matplotlib.pyplot as plt

"""
Code from https://github.com/jiarong-wu/multilayer_breaking/blob/bab0b535be26eaf3021bf5a34c9bf20d1962e483/postprocessing/mlpython/breaking.py
"""

""" The detection function """


def detect_ridges(gray, sigma=1.0):
    H_elems = hessian_matrix(
        gray, sigma=sigma, order="rc", use_gaussian_derivatives=False
    )
    maxima_ridges, minima_ridges = hessian_matrix_eigvals(H_elems)
    return maxima_ridges, minima_ridges


def detect_slope(grad, threshold):
    return np.logical_not(grad > threshold)


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


""" Detect breaking fronts based on geometry of surface elevation (eta) and compute the statistics 
    of the velocity of the detected front. """


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
