#!/usr/bin/env python3
"""
VTKHDF inspector and checker.

`print_vtkhdf_info()` prints the coordinates and topology of a file, which is
what the `~~~pythonplot` blocks of the tests display. The `--check` mode below
goes further: the tests set analytical fields (`u.x[]=x`, `f[]=point.level`),
so the file can be *verified*, not merely inspected. The report is a compact,
deterministic summary suitable for use as a Basilisk reference (`.ref`) file:
every quantity is an integer, a rounded float, or a thresholded verdict, so it
is stable across builds and optimisation levels.
"""
import sys
import os

try:
  import h5py
except ImportError:  # pragma: no cover - depends on the local installation
  # The spack h5py is built against a different Python than the default
  # interpreter here, so fall back to a throw-away environment. The re-exec
  # keeps the callers (and the `~~~pythonplot` blocks) calling plain python3.
  import shutil
  if os.environ.get('VTKHDF_NO_UV') or not shutil.which('uv'):
    raise
  os.environ['VTKHDF_NO_UV'] = '1'
  os.execvp('uv', ['uv', 'run', '--no-project', '--with', 'h5py',
                   'python3', os.path.abspath(__file__)] + sys.argv[1:])

# VTK cell types we expect from the output routines
VTK_TYPE_NAMES = {3: 'LINE', 4: 'POLY_LINE', 7: 'POLYGON',
                  8: 'PIXEL', 9: 'QUAD', 11: 'VOXEL', 12: 'HEXAHEDRON'}


def read_vtkhdf(filename):
  """Read a VTKHDF file and return its points and topology.

  The file is written one partition per MPI rank: `NumberOfPoints`,
  `NumberOfCells` and `NumberOfConnectivityIds` hold one entry per partition,
  `Offsets` restarts from 0 in each of them, and the connectivity indices are
  local to their own partition. They are concatenated here into the single
  global mesh a reader would assemble, so that everything downstream can
  ignore the decomposition.
  """
  with h5py.File(filename, 'r') as h5f:
    # VTKHDF format stores data in the VTKHDF group
    if 'VTKHDF' not in h5f:
      raise ValueError(f"File {filename} is not a valid VTKHDF file "
                       "(missing VTKHDF group)")

    ghdf = h5f['VTKHDF']

    points_data = ghdf['Points'][:]
    points = [(p[0], p[1], p[2]) if len(p) == 3 else (p[0], p[1], 0.0)
              for p in points_data]

    raw_connectivity = [int(i) for i in ghdf['Connectivity'][:].flatten()]
    raw_offsets = [int(i) for i in ghdf['Offsets'][:].flatten()]
    types = [int(t) for t in ghdf['Types'][:].flatten()]

    n_points = [int(n) for n in ghdf['NumberOfPoints'][:]]
    n_cells = [int(n) for n in ghdf['NumberOfCells'][:]]
    n_ids = [int(n) for n in ghdf['NumberOfConnectivityIds'][:]]

  # Shift each partition into the global numbering: its connectivity indices
  # by the points already seen, its offsets by the ids already seen.
  connectivity, offsets = [], [0]
  point_base, id_base, off_base = 0, 0, 0
  for part in range(len(n_cells)):
    connectivity += [i + point_base
                     for i in raw_connectivity[id_base:id_base + n_ids[part]]]
    # each partition contributes NumberOfCells+1 offsets, the first one 0
    local = raw_offsets[off_base:off_base + n_cells[part] + 1]
    offsets += [o + id_base for o in local[1:]]
    point_base += n_points[part]
    id_base += n_ids[part]
    off_base += n_cells[part] + 1

  return {
    'num_points': len(points),
    'num_cells': len(offsets) - 1,
    'num_partitions': len(n_cells),
    'points': points,
    'connectivity': connectivity,
    'offsets': offsets,
    'types': types
  }


def print_vtkhdf_info(filename):
  """Print coordinates and topology for a VTKHDF file"""
  print(f"\n{'='*70}")
  print(f"File: {os.path.basename(filename)}")
  print(f"{'='*70}")

  try:
    data = read_vtkhdf(filename)
  except Exception as e:
    print(f"Error reading file: {e}")
    return

  # Print coordinates
  print(f"\nCOORDINATES ({data['num_points']} points):")
  print(f"{'Index':<6} | {'X':<15} | {'Y':<15} | {'Z':<15}")
  print("-" * 55)
  for i, (x, y, z) in enumerate(data['points']):
    print(f"{i:<6} | {x:<15.8f} | {y:<15.8f} | {z:<15.8f}")

  # Print topology
  print(f"\nTOPOLOGY ({data['num_cells']} cells):")
  print(f"{'Cell':<6} | {'Point Indices':<30} | {'Type':<6}")
  print("-" * 50)
  for c_idx in range(data['num_cells']):
    start_ptr = data['offsets'][c_idx]
    end_ptr = data['offsets'][c_idx+1]
    point_indices = data['connectivity'][start_ptr:end_ptr]

    # Check for invalid indices
    invalid = [pid for pid in point_indices
               if pid < 0 or pid >= data['num_points']]
    if invalid:
      indices_str = f"{list(point_indices)} ⚠ INVALID: {invalid}"
    else:
      indices_str = str(list(point_indices))

    print(f"{c_idx:<6} | {indices_str:<30} | {data['types'][c_idx]:<6}")

  print(f"{'='*70}\n")


def read_vtkhdf_arrays(filename):
  """Read a VTKHDF file together with its cell-centred attributes.

  `read_vtkhdf()` above deliberately reports only the mesh, because that is
  what the documentation blocks display. The checks below also need the field
  data, so it is collected here.
  """
  data = read_vtkhdf(filename)
  cell_data = {}
  with h5py.File(filename, 'r') as h5f:
    group = h5f['VTKHDF'].get('CellData')
    for name in (group or {}):
      values = group[name][:]
      # Scalars are stored as (N, 1) and vectors as (N, 3); flatten the
      # former so a scalar reads as a plain list of numbers.
      if values.ndim == 2 and values.shape[1] == 1:
        cell_data[name] = [float(v[0]) for v in values]
      else:
        cell_data[name] = [tuple(float(c) for c in row) for row in values]
  data['cell_data'] = cell_data
  return data


def cell_vertices(data, c_idx):
  """Return the list of vertex coordinates of cell `c_idx`."""
  start, end = data['offsets'][c_idx], data['offsets'][c_idx + 1]
  return [data['points'][i] for i in data['connectivity'][start:end]]


def check_topology(data):
  """Report point/cell counts, cell types and connectivity validity.

  Returns the number of failed checks.
  """
  failures = 0
  print(f"points                 {data['num_points']}")
  print(f"cells                  {data['num_cells']}")

  counts = {}
  for t in data['types']:
    counts[t] = counts.get(t, 0) + 1
  types_str = ', '.join(f"{VTK_TYPE_NAMES.get(t, t)}={n}"
                        for t, n in sorted(counts.items()))
  print(f"cell types             {types_str}")

  n_bad = sum(1 for i in data['connectivity']
              if i < 0 or i >= data['num_points'])
  print(f"connectivity in range  {'OK' if n_bad == 0 else f'FAIL ({n_bad} bad)'}")
  failures += n_bad > 0

  # The last offset must account for exactly the whole connectivity array
  consistent = (data['offsets'][-1] == len(data['connectivity'])
                and len(data['types']) == data['num_cells'])
  print(f"offsets consistent     {'OK' if consistent else 'FAIL'}")
  failures += not consistent

  # Every point must be referenced by at least one cell
  orphans = data['num_points'] - len(set(data['connectivity']))
  print(f"orphan points          {'OK' if orphans == 0 else f'FAIL ({orphans})'}")
  failures += orphans > 0

  failures += check_cell_ordering(data)
  return failures


def check_cell_ordering(data):
  """Verify the vertex ordering inside each QUAD/HEXAHEDRON cell.

  A centroid is permutation-invariant, so comparing centroids cannot detect a
  cell whose vertices were emitted in the wrong order -- yet such a cell is
  rendered as a bow-tie by ParaView. Grid cells are axis-aligned here, so the
  expected VTK ordering can be reconstructed from the cell bounding box and
  compared vertex by vertex.

  Returns the number of failed checks.
  """
  n_bad = 0
  n_checked = 0
  for c in range(data['num_cells']):
    t = data['types'][c]
    if t not in (9, 12):  # QUAD, HEXAHEDRON
      continue
    verts = cell_vertices(data, c)
    lo = [min(v[d] for v in verts) for d in range(3)]
    hi = [max(v[d] for v in verts) for d in range(3)]

    if t == 9:
      # A quad is flat in one direction: in a 2D run that is z, in a slice
      # file it is the slice normal. Express the expected ordering in the
      # two remaining axes, so every quad is checked rather than skipped.
      # The convention is counter-clockwise from the low corner.
      plane = [d for d in range(3) if hi[d] > lo[d]]
      if len(plane) != 2:
        continue  # degenerate cell, nothing meaningful to check
      d0, d1 = plane
      corners = [(lo[d0], lo[d1]), (hi[d0], lo[d1]),
                 (hi[d0], hi[d1]), (lo[d0], hi[d1])]
      expected = []
      for c0, c1 in corners:
        point = list(lo)
        point[d0], point[d1] = c0, c1
        expected.append(tuple(point))
    else:
      if any(hi[d] == lo[d] for d in range(3)):
        continue  # degenerate hexahedron
      # bottom face counter-clockwise, then the top face
      expected = [(lo[0], lo[1], lo[2]), (hi[0], lo[1], lo[2]),
                  (hi[0], hi[1], lo[2]), (lo[0], hi[1], lo[2]),
                  (lo[0], lo[1], hi[2]), (hi[0], lo[1], hi[2]),
                  (hi[0], hi[1], hi[2]), (lo[0], hi[1], hi[2])]

    n_checked += 1
    if len(verts) != len(expected) or any(
        abs(verts[i][d] - expected[i][d]) > 1e-12
        for i in range(len(expected)) for d in range(3)):
      n_bad += 1

  if n_checked:
    print(f"vertex ordering        "
          f"{'OK' if n_bad == 0 else f'FAIL ({n_bad}/{n_checked} cells)'}")
  return n_bad > 0


def check_analytical_fields(data, box=None):
  """Verify the analytical fields written by `test_output_vtkhdf[_box].c`.

  Those tests set `u.x[]=x`, `u.y[]=y`, `u.z[]=z` and `f[]=point.level`, so
  the cell-centred vector must reproduce the centroid of the cell computed
  from its own vertices. This checks the field data and the topology at once:
  if the connectivity were scrambled, the centroids would no longer match.

  Returns the number of failed checks.
  """
  import math
  failures = 0
  cd = data['cell_data']

  # u must equal the cell centroid, recomputed from the connectivity
  if 'u.x' in cd:
    max_err = 0.0
    for c in range(data['num_cells']):
      verts = cell_vertices(data, c)
      centroid = [sum(v[d] for v in verts) / len(verts) for d in range(3)]
      written = cd['u.x'][c]
      # All three components are compared, in 2D as well as in 3D: a 2D run
      # writes u.z as 0 and its cells are flat in z, so the centroid is 0
      # too and the comparison stays meaningful.
      for d in range(3):
        max_err = max(max_err, abs(written[d] - centroid[d]))
    ok = max_err < 1e-12
    print(f"u == cell centroid     {'OK' if ok else f'FAIL (max err {max_err:.3e})'}")
    failures += not ok

  # f[]=point.level, so f + log2(dx) must be the same constant for every cell,
  # whatever the refinement. This validates f across all levels present.
  if 'f' in cd:
    residuals, levels = set(), set()
    for c in range(data['num_cells']):
      verts = cell_vertices(data, c)
      # Use the largest extent rather than the x one: cells of a slice file
      # are flat in the slice-normal direction, whose extent is 0.
      dx = max(max(v[d] for v in verts) - min(v[d] for v in verts)
               for d in range(3))
      residuals.add(round(cd['f'][c] + math.log2(dx), 6))
      levels.add(int(cd['f'][c]))
    ok = len(residuals) == 1
    print(f"levels present         {sorted(levels)}")
    print(f"f == cell level        "
          f"{'OK' if ok else f'FAIL (residuals {sorted(residuals)})'}")
    failures += not ok

  # p[]=pid(); in serial this is 0 everywhere, under MPI one value per rank
  if 'p' in cd:
    print(f"p (pid) values         {sorted(set(int(v) for v in cd['p']))}")

  # the box variant must not emit anything outside its region
  if box is not None:
    failures += check_cells_inside_box(data, box)
  return failures


def check_cells_inside_box(data, box):
  """Verify that no cell's centroid lies outside `box = (lo, hi)`.

  Shared by `check_analytical_fields()` (regular field export) and
  `check_cardioid()` (VOF facets export) -- both `output_vtkhdf_box()` and
  `output_facets_vtkhdf_box()` restrict on the same cell-center-in-box test,
  so one check covers both writers.

  Returns the number of failed checks (0 or 1).
  """
  outside = 0
  for c in range(data['num_cells']):
    verts = cell_vertices(data, c)
    centroid = [sum(v[d] for v in verts) / len(verts) for d in range(3)]
    ndim = len(box[0])
    if any(not (box[0][d] <= centroid[d] <= box[1][d]) for d in range(ndim)):
      outside += 1
  print(f"cells inside box       "
        f"{'OK' if outside == 0 else f'FAIL ({outside} outside)'}")
  return outside > 0


def check_cardioid(data, a=0.15, box=None):
  """Verify the interface written by `test_output_vtkhdf_facets.c`.

  That test builds a VOF fraction from the cardioid
  ``(x^2+y^2-2*a*x)^2 = 4*a^2*(x^2+y^2)``, so every facet vertex must lie on
  that curve to within about one cell width. The distance is estimated as
  ``|g|/|grad g|``, the first-order distance to the zero level set.

  Returns the number of failed checks.
  """
  import math
  failures = 0

  # cell width of the finest facet, used as the tolerance scale
  widths = []
  for c in range(data['num_cells']):
    verts = cell_vertices(data, c)
    widths.append(max(max(v[d] for v in verts) - min(v[d] for v in verts)
                      for d in range(2)))
  dx_max = max(widths)

  worst = 0.0
  for (x, y, _z) in data['points']:
    common = x * x + y * y - 2 * a * x
    g = common * common - 4 * a * a * (x * x + y * y)
    # analytical gradient of g
    gx = 2 * common * (2 * x - 2 * a) - 8 * a * a * x
    gy = 2 * common * (2 * y) - 8 * a * a * y
    norm = math.hypot(gx, gy)
    if norm > 0:
      worst = max(worst, abs(g) / norm)

  # Report the verdict rather than the raw distance: the distance is a
  # floating-point quantity that may drift with the optimisation level, which
  # would make the reference file brittle.
  ok = worst < dx_max
  print(f"facet width (max)      {dx_max:.6f}")
  print(f"vertices on cardioid   "
        f"{'OK (within one cell)' if ok else f'FAIL (max distance {worst:.2e})'}")
  failures += not ok

  # The curvature is written alongside the facets; it must at least be finite
  # everywhere, which `curvature()` does not guarantee for interface-free cells.
  kappa = next((v for k, v in data['cell_data'].items()
                if k.lower().startswith('kappa')), None)
  if kappa is not None:
    n_finite = sum(1 for k in kappa if math.isfinite(k))
    print(f"kappa finite           "
          f"{'OK' if n_finite == len(kappa) else f'FAIL ({n_finite}/{len(kappa)})'}")
    failures += n_finite != len(kappa)
  else:
    print("kappa finite           FAIL (no curvature attribute)")
    failures += 1

  # the box variant must not emit anything outside its region
  if box is not None:
    failures += check_cells_inside_box(data, box)
  return failures


def check_plane(data, axis, value):
  """Verify that a slice file is genuinely planar.

  `output_vtkhdf_slice()` is asked for the plane `axis = value`, so every
  vertex it writes must lie exactly on it. This is the property that
  distinguishes a slice from a full domain dump, so it is checked separately.

  Returns the number of failed checks.
  """
  d = 'xyz'.index(axis)
  worst = max(abs(p[d] - value) for p in data['points'])
  ok = worst < 1e-12
  label = f"points on plane {axis}={value:g}"
  print(f"{label:<22} {'OK' if ok else f'FAIL (max deviation {worst:.3e})'}")
  return 0 if ok else 1


def check(filename, box=None, cardioid=False, plane=None):
  """Run every applicable check on `filename` and report a single verdict.

  Returns 0 when all checks pass, 1 otherwise, so the caller can propagate an
  exit status.
  """
  print(f"--- {os.path.basename(filename)}")
  if not os.path.exists(filename):
    print("MISSING FILE")
    return 1

  data = read_vtkhdf_arrays(filename)
  failures = check_topology(data)
  if plane is not None:
    failures += check_plane(data, plane[0], plane[1])
  if cardioid:
    failures += check_cardioid(data, box=box)
  else:
    failures += check_analytical_fields(data, box=box)

  print(f"verdict                {'PASS' if failures == 0 else f'FAIL ({failures})'}")
  return 0 if failures == 0 else 1


if __name__ == "__main__":
  args = [a for a in sys.argv[1:] if not a.startswith('--')]
  flags = {a for a in sys.argv[1:] if a.startswith('--')}

  if not args:
    print("Usage: python3 test_output_vtkhdf.py [--check] [--box=X0,Y0,X1,Y1] "
          "[--plane=AXIS=VALUE] [--cardioid] <file1.hdf> ...")
    sys.exit(1)

  if '--check' in flags or '--cardioid' in flags or any(
      f.startswith(('--box=', '--plane=')) for f in flags):
    box, plane = None, None
    for f in flags:
      if f.startswith('--box='):
        v = [float(t) for t in f[len('--box='):].split(',')]
        half = len(v) // 2
        box = (v[:half], v[half:])
      if f.startswith('--plane='):
        axis, value = f[len('--plane='):].split('=')
        plane = (axis, float(value))
    status = 0
    for filename in args:
      status |= check(filename, box=box, plane=plane,
                      cardioid='--cardioid' in flags)
    sys.exit(status)

  for filename in args:
    if os.path.exists(filename):
      print_vtkhdf_info(filename)
    else:
      print(f"Error: {filename} not found")
