#!/usr/bin/env python3
"""
Simple XDMF/HDF5 inspector: prints coordinates and topology

Run without flags it reports the mesh, as before. Run with `--check` it
verifies the file against the analytical fields set by `test_output_xmf*.c`
and exits non-zero on the first failure, so it can be used as the pass/fail
criterion of a `make <test>.tst` run.
"""
import os
import sys
import xml.etree.ElementTree as ET

# h5py is needed to read the heavy data. Where it is missing but `uv` is
# available, re-run ourselves under `uv run --with h5py` rather than failing:
# the test is then usable on a machine whose system Python has no h5py,
# without the .c file having to know anything about it.
try:
  import h5py
except ImportError:  # pragma: no cover - depends on the local installation
  import shutil
  if os.environ.get('XMF_NO_UV') or not shutil.which('uv'):
    raise
  os.environ['XMF_NO_UV'] = '1'
  os.execvp('uv', ['uv', 'run', '--no-project', '--with', 'h5py',
                   'python3', os.path.abspath(__file__)] + sys.argv[1:])

# XDMF topology names, with the number of vertices each cell carries
XDMF_CELL_POINTS = {
  'Polyvertex': 1, 'Polyline': 2, 'Triangle': 3, 'Quadrilateral': 4,
  'Tetrahedron': 4, 'Hexahedron': 8,
}

# A `Mixed` topology stores one flat stream of integers instead of a table of
# fixed width: each cell starts with its own type code, and the three variable
# sized types then give their vertex count before the vertices themselves.
XDMF_MIXED_CODES = {
  1: 'Polyvertex', 2: 'Polyline', 3: 'Polygon', 4: 'Triangle',
  5: 'Quadrilateral', 6: 'Tetrahedron', 7: 'Pyramid', 8: 'Wedge',
  9: 'Hexahedron',
}
XDMF_MIXED_COUNTED = ('Polyvertex', 'Polyline', 'Polygon')


def decode_mixed_topology(stream, num_cells):
  """Split a `Mixed` connectivity stream into one vertex list per cell.

  Returns `(connectivity, cell_types)`. Raises `ValueError` when the stream
  does not decode into exactly `num_cells` cells, which is itself a defect
  worth reporting: it means the topology is unreadable, not merely wrong.
  """
  connectivity, cell_types, i = [], [], 0
  while i < len(stream) and len(connectivity) < num_cells:
    kind = XDMF_MIXED_CODES.get(stream[i])
    if kind is None:
      raise ValueError(f"unknown mixed topology code {stream[i]} at offset {i}")
    i += 1
    if kind in XDMF_MIXED_COUNTED:
      n = stream[i]
      i += 1
    else:
      n = XDMF_CELL_POINTS[kind]
    connectivity.append(stream[i:i + n])
    cell_types.append(kind)
    i += n
  if len(connectivity) != num_cells or i != len(stream):
    raise ValueError(f"mixed topology decoded into {len(connectivity)} cells "
                     f"and {i}/{len(stream)} entries, expected {num_cells}")
  return connectivity, cell_types


def read_xdmf(xmf_filename):
  """Read XDMF file and return points and topology"""

  # Parse the XDMF XML file
  tree = ET.parse(xmf_filename)
  root = tree.getroot()

  # Find the Grid element
  domain = root.find('.//Domain')
  grid = domain.find('.//Grid')

  # Get the HDF5 filename from the XML
  h5_filename = None
  for data_item in grid.findall('.//DataItem'):
    text = data_item.text.strip()
    if '.h5:' in text:
      h5_filename = text.split(':')[0]
      break

  if h5_filename is None:
    raise ValueError("Could not find HDF5 filename in XDMF file")

  # Make path relative to XDMF file location
  xmf_dir = os.path.dirname(xmf_filename)
  h5_path = os.path.join(xmf_dir, h5_filename) if xmf_dir else h5_filename

  # Read data from HDF5 file
  with h5py.File(h5_path, 'r') as h5f:
    # Read points (geometry)
    points_data = h5f['/Geometry/Points'][:]
    num_points = len(points_data)
    points = [(p[0], p[1], p[2]) for p in points_data]

    # Read topology (connectivity)
    # This is stored as a 2D array: num_cells x points_per_cell
    topology_data = h5f['/Topology'][:]
    num_cells = len(topology_data)
    points_per_cell = topology_data.shape[1]

    # Convert to list of lists for easier handling (convert to plain Python ints)
    connectivity = [[int(x) for x in row] for row in topology_data]

  return {
    'num_points': num_points,
    'num_cells': num_cells,
    'points': points,
    'connectivity': connectivity,
    'points_per_cell': points_per_cell,
    'h5_file': h5_filename
  }


def print_xdmf_info(xmf_filename):
  """Print coordinates and topology for an XDMF file"""
  print(f"\n{'='*70}")
  print(f"File: {os.path.basename(xmf_filename)}")
  print(f"{'='*70}")

  data = read_xdmf(xmf_filename)

  print(f"HDF5 file: {data['h5_file']}")
  print(f"Points per cell: {data['points_per_cell']}")

  # Print coordinates
  print(f"\nCOORDINATES ({data['num_points']} points):")
  print(f"{'Index':<6} | {'X':<15} | {'Y':<15} | {'Z':<15}")
  print("-" * 55)
  for i, (x, y, z) in enumerate(data['points']):
    print(f"{i:<6} | {x:<15.8f} | {y:<15.8f} | {z:<15.8f}")

  # Print topology
  print(f"\nTOPOLOGY ({data['num_cells']} cells):")
  print(f"{'Cell':<6} | {'Point Indices':<40}")
  print("-" * 50)
  for c_idx in range(data['num_cells']):
    point_indices = data['connectivity'][c_idx]

    # Check for invalid indices
    invalid = [pid for pid in point_indices if pid < 0 or pid >= data['num_points']]
    if invalid:
      indices_str = f"{point_indices} ⚠ INVALID: {invalid}"
    else:
      indices_str = str(point_indices)

    print(f"{c_idx:<6} | {indices_str:<40}")

  print(f"{'='*70}\n")


def read_xdmf_arrays(xmf_filename):
  """Read an XDMF file together with its cell-centred attributes.

  `read_xdmf()` above deliberately reports only the mesh, because that is what
  the documentation blocks display. The checks below also need the field data
  and the declared topology type, so they are collected here. The dataset
  paths are taken from the XML rather than assumed, so that a writer emitting
  its attributes elsewhere in the file is followed rather than missed.
  """
  tree = ET.parse(xmf_filename)
  grid = tree.getroot().find('.//Domain').find('.//Grid')

  topology = grid.find('.//Topology')
  topology_type = topology.get('TopologyType')

  def h5_path_of(data_item):
    """Return the dataset path of a DataItem, e.g. `/Cells/f`."""
    text = data_item.text.strip()
    # `&HeavyData;` expands to `domain.h5:`, and ElementTree leaves the
    # entity unexpanded, so split on whichever separator is present.
    return text.split(';')[-1].split(':')[-1].strip()

  xmf_dir = os.path.dirname(xmf_filename)
  h5_name = read_xdmf(xmf_filename)['h5_file']
  h5_path = os.path.join(xmf_dir, h5_name) if xmf_dir else h5_name

  cell_data = {}
  with h5py.File(h5_path, 'r') as h5f:
    points = [tuple(float(c) for c in p) for p in h5f['/Geometry/Points'][:]]
    raw = h5f[h5_path_of(topology.find('.//DataItem'))][:]
    if topology_type == 'Mixed':
      connectivity, cell_types = decode_mixed_topology(
        [int(i) for i in raw.reshape(-1)],
        int(topology.get('NumberOfElements')))
    else:
      connectivity = [[int(i) for i in row] for row in raw]
      cell_types = [topology_type] * len(connectivity)
    for attribute in grid.findall('.//Attribute'):
      name = attribute.get('Name')
      item = attribute.find('.//DataItem')
      if item is None:
        continue
      values = h5f[h5_path_of(item)][:]
      # Scalars are stored as (N, 1) and vectors as (N, 3); flatten the
      # former so a scalar reads as a plain list of numbers.
      if values.ndim == 2 and values.shape[1] == 1:
        cell_data[name] = [float(v[0]) for v in values]
      else:
        cell_data[name] = [tuple(float(c) for c in row) for row in values]

  return {
    'num_points': len(points),
    'num_cells': len(connectivity),
    'points': points,
    'connectivity': connectivity,
    'topology_type': topology_type,
    'cell_types': cell_types,
    'cell_data': cell_data,
  }


def cell_vertices(data, c_idx):
  """Return the list of vertex coordinates of cell `c_idx`."""
  return [data['points'][i] for i in data['connectivity'][c_idx]]


def check_topology(data):
  """Report point/cell counts, topology type and connectivity validity.

  Returns the number of failed checks.
  """
  failures = 0
  print(f"points                 {data['num_points']}")
  print(f"cells                  {data['num_cells']}")
  print(f"topology type          {data['topology_type']}")

  # The declared topology type must agree with the width of the connectivity
  # array, otherwise a reader would walk the cells with the wrong stride.
  # In a Mixed topology each cell declares its own type, so the width has to
  # agree with that type cell by cell rather than globally. Report the types
  # that were decoded, since that is what the writer actually chose.
  if data['topology_type'] == 'Mixed':
    counts = {}
    for kind in data['cell_types']:
      counts[kind] = counts.get(kind, 0) + 1
    print("cell types             "
          + ' '.join(f'{k}={n}' for k, n in sorted(counts.items())))
    bad = sum(1 for kind, row in zip(data['cell_types'], data['connectivity'])
              if kind not in XDMF_MIXED_COUNTED
              and len(row) != XDMF_CELL_POINTS[kind])
    ok = bad == 0
    print(f"points per cell        {'OK' if ok else f'FAIL ({bad} cells)'}")
  else:
    expected = XDMF_CELL_POINTS.get(data['topology_type'])
    widths = set(len(row) for row in data['connectivity'])
    ok = expected is not None and widths == {expected}
    print(f"points per cell        "
          f"{'OK' if ok else f'FAIL (expected {expected}, got {sorted(widths)})'}")
  failures += not ok

  n_bad = sum(1 for row in data['connectivity'] for i in row
              if i < 0 or i >= data['num_points'])
  print(f"connectivity in range  {'OK' if n_bad == 0 else f'FAIL ({n_bad} bad)'}")
  failures += n_bad > 0

  # Every point must be referenced by at least one cell
  used = set(i for row in data['connectivity'] for i in row)
  orphans = data['num_points'] - len(used)
  print(f"orphan points          {'OK' if orphans == 0 else f'FAIL ({orphans})'}")
  failures += orphans > 0

  failures += check_cell_ordering(data)
  return failures


def check_cell_ordering(data):
  """Verify the vertex ordering inside each quadrilateral/hexahedral cell.

  A centroid is permutation-invariant, so comparing centroids cannot detect a
  cell whose vertices were emitted in the wrong order -- yet such a cell is
  rendered as a bow-tie by ParaView. Grid cells are axis-aligned here, so the
  expected ordering can be reconstructed from the cell bounding box and
  compared vertex by vertex.

  Returns the number of failed checks.
  """
  kind = data['topology_type']
  if kind not in ('Quadrilateral', 'Hexahedron'):
    return 0

  n_bad = 0
  n_checked = 0
  for c in range(data['num_cells']):
    verts = cell_vertices(data, c)
    lo = [min(v[d] for v in verts) for d in range(3)]
    hi = [max(v[d] for v in verts) for d in range(3)]

    if kind == 'Quadrilateral':
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
  """Verify the analytical fields written by `test_output_xmf[_box].c`.

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
  `check_cardioid()` (VOF facets export) -- both `output_xmf_box()` and
  `output_facets_xmf_box()` restrict on the same cell-center-in-box test, so
  one check covers both writers.

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
  """Verify the interface written by `test_output_xmf_facets.c`.

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

  # The attributes written alongside the facets must be finite everywhere
  if not data['cell_data']:
    print("attributes finite      FAIL (no cell attributes)")
    failures += 1
  else:
    for name, values in sorted(data['cell_data'].items()):
      if isinstance(values[0], (int, float)):
        n_finite = sum(1 for v in values if math.isfinite(v))
        label = f"{name} finite"
        ok_f = n_finite == len(values)
        print(f"{label:<22} {'OK' if ok_f else f'FAIL ({n_finite}/{len(values)})'}")
        failures += not ok_f

  # the box variant must not emit any facet whose cell lies outside its region
  if box is not None:
    failures += check_cells_inside_box(data, box)
  return failures


def check_plane(data, axis, value):
  """Verify that a slice file is genuinely planar.

  `output_xmf_slice()` is asked for the plane `axis = value`, so every vertex
  it writes must lie exactly on it. This is the property that distinguishes a
  slice from a full domain dump, so it is checked separately.

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

  data = read_xdmf_arrays(filename)
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
    print("Usage: python3 test_output_xmf.py [--check] [--box=X0,Y0,X1,Y1] "
          "[--plane=AXIS=VALUE] [--cardioid] <file1.xmf> ...")
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
      print_xdmf_info(filename)
    else:
      print(f"Error: {filename} not found")
