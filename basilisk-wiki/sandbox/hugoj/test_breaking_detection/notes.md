# 1) Test plan

## a. Unit tests per function, in isolation

get_bins: exact input match for the three kp cases (use np.isclose on
2*np.pi/10 etc., since float equality is fragile — this itself might be a bug
worth fixing with a tolerance), plus a test that an unmatched kp raises rather
than silently returning None. detect_slope: simple synthetic gradient array
with known values above/below threshold, check boolean output shape and
correctness at the threshold boundary (> vs >=). detect_ridges: use a synthetic
analytic surface with a known ridge (e.g. a 1D Gaussian ridge extruded along
one axis, or -x² parabola), where you can derive the expected sign/location of
maxima and minima curvature eigenvalues analytically. This is the anchor test
that pins down which of the two returned arrays is "maxima" vs "minima" and
confirms axis/order conventions from skimage.

## b. simple_mapping integration tests

Build small synthetic eta, ux, uy (e.g. 20×20 or 32×32) with a known,
hand-placed ridge/breaking front and known velocity field, so the expected
histogram is computable by hand or via a slow reference implementation. Test
method=0 and method=1 paths separately (once method=1 is fixed). Test
EXTRA_FILTER=True/False produce different, predictable results. Test edge
cases: flat eta (no ridges → all-zero histogram), non-square array (to catch
bug #4), array with NaNs (does it propagate or crash?). Test bins edge cases
(empty, single bin, unsorted bins).

## c. Regression / golden-master tests

Before refactoring, run the current (bug-fixed) implementation on a couple of
real or realistic synthetic netCDF snapshots, save the output arrays as
fixtures (e.g. .npy or a small .nc), and assert the new vectorized/dask version
reproduces them within a numerical tolerance. This is the safety net that lets
you refactor aggressively without fear.

## d. Property-based tests (optional, hypothesis)

E.g. output histogram total count is monotonic/bounded by number of detected
crest pixels; ridge detection is invariant to a constant additive offset in eta
(since curvature shouldn't care about mean height, only the height filter
should).

## e. Performance regression test

A lightweight benchmark test (e.g. pytest-benchmark) on a moderately sized
field/time series, so future changes don't silently reintroduce an O(N²) Python
loop.

______________________________________________________________________________

# 2) Speed-up plan

## a. Kill the double Python loop first (biggest win, orthogonal to xarray) 

The "one-sided edge" detection (a_[i][j] = 1 if b_[i][j-1]>0 and b_[i][j+1]==0) is
a classic case for pure array slicing: compare b_ shifted left vs shifted right
along the relevant axis, using numpy slicing or np.roll/padding, no explicit
loop needed. This alone should give a large speedup with zero change in
dependencies.

## b. Decide the right axis for "time" in your netCDF 

For ridge detection over a
time series of eta(t, y, x), the natural parallelism is over t — each time
slice is independent. That's the dimension to chunk with dask, not y/x (since
hessian_matrix with sigma needs neighbouring pixels, so spatial chunking
requires halo/overlap handling — better avoided unless individual frames are
huge).

## c. Wrapping in xr.apply_ufunc

Load eta, ux, uy as xarray.DataArrays with dask-backed chunks along time (chunk
size 1 or a small batch per chunk, spatial dims unchunked). Wrap detect_ridges
(and the vectorized simple_mapping core) with xr.apply_ufunc(...,
input_core_dims=[["y","x"]], output_core_dims=[["y","x"],["y","x"]],
vectorize=True, dask="parallelized", output_dtypes=[...]) so each time slice
gets dispatched as a delayed task. Because skimage.feature.hessian_matrix isn't
dask-aware itself, apply_ufunc with dask="parallelized" + vectorize=True is the
pragmatic bridge — it runs the numpy function per chunk under the hood.

## d. Histogram binning step

np.histogram per time step in a loop is also parallelizable across time; could
either keep it as a fast numpy step inside the same apply_ufunc-wrapped
function (returning per-time histograms as a new velocity_bin dimension), or
use dask.array reductions if we want a running/aggregated histogram across the
whole series without materializing everything at once.

## e. Other cheap wins

Avoid recomputing np.var(eta) globally per call if it can be computed once and
reused/cached, or made explicitly per-frame via xr.DataArray.var(dim=["y","x"])
broadcast back — this also disambiguates bug #5. Use float32 instead of float64
for the fields if memory bandwidth is the bottleneck (common in these
pipelines), verified via the regression tests. Suggested order of work Lock
down intended behaviour with unit tests on the current logic (accepting some
bugs as "documented, to be fixed"). Fix the identified bugs, updating tests to
reflect intended (not buggy) behaviour, with your sign-off on each fix.
Vectorize the Python loop, verify against regression fixtures. Wrap in
xarray/dask for time-parallelism, verify again against the same fixtures. Add
the performance regression test to lock in the gains.
