#!/usr/bin/env python3
"""
Simple VTU inspector: prints coordinates and topology
"""
import struct
import sys
import os

def read_vtu(filename):
  """Read VTU file and return points and topology"""
  with open(filename, 'rb') as f:
    content = f.read()

  # Parse XML metadata
  xml_part = content.split(b'<AppendedData')[0].decode('utf-8')
  pi_idx = xml_part.find('<Piece')
  
  np_idx = xml_part.find('NumberOfPoints="', pi_idx) + 16
  np_end = xml_part.find('"', np_idx)
  num_points = int(xml_part[np_idx:np_end])
  
  nc_idx = xml_part.find('NumberOfCells="', pi_idx) + 15
  nc_end = xml_part.find('"', nc_idx)
  num_cells = int(xml_part[nc_idx:nc_end])

  # Helper to find data offsets
  def get_offset(name, parent_tag=None):
    if parent_tag:
      start = xml_part.find(f'<{parent_tag}>')
      end = xml_part.find(f'</{parent_tag}>')
      section = xml_part[start:end]
    else:
      section = xml_part
    tag_start = section.find(f'Name="{name}"') if name else section.find('<DataArray')
    off_idx = section.find('offset="', tag_start) + 8
    off_end = section.find('"', off_idx)
    return int(section[off_idx:off_end])

  offset_points = get_offset(None, "Points")
  offset_offsets = get_offset("offsets", "Cells")
  offset_types = get_offset("types", "Cells")
  offset_connectivity = get_offset("connectivity", "Cells")
  
  # Find binary data start
  appended_start = content.find(b'<AppendedData encoding="raw">')
  data_start = content.find(b'_', appended_start) + 1
  
  # Read Points
  point_data_start = data_start + offset_points
  coords_start = point_data_start + 8
  points = []
  for i in range(num_points):
    p_idx = coords_start + i * 24
    x, y, z = struct.unpack('<ddd', content[p_idx:p_idx+24])
    points.append((x, y, z))
  
  # Read Offsets
  off_block_start = data_start + offset_offsets
  off_len = struct.unpack('<Q', content[off_block_start:off_block_start+8])[0]
  offsets = list(struct.unpack(f'<{off_len//8}q', content[off_block_start+8:off_block_start+8+off_len]))
  
  # Read Types
  type_block_start = data_start + offset_types
  type_len = struct.unpack('<Q', content[type_block_start:type_block_start+8])[0]
  cell_types = list(struct.unpack(f'<{type_len}b', content[type_block_start+8:type_block_start+8+type_len]))
  
  # Read Connectivity
  conn_block_start = data_start + offset_connectivity
  conn_len = struct.unpack('<Q', content[conn_block_start:conn_block_start+8])[0]
  connectivity = list(struct.unpack(f'<{conn_len//8}q', content[conn_block_start+8:conn_block_start+8+conn_len]))

  return {
    'num_points': num_points,
    'num_cells': num_cells,
    'points': points,
    'offsets': offsets,
    'types': cell_types,
    'connectivity': connectivity
  }

def print_vtu_info(filename):
  """Print coordinates and topology for a VTU file"""
  print(f"\n{'='*70}")
  print(f"File: {os.path.basename(filename)}")
  print(f"{'='*70}")
  
  data = read_vtu(filename)
  
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
    start_ptr = 0 if c_idx == 0 else data['offsets'][c_idx-1]
    end_ptr = data['offsets'][c_idx]
    point_indices = data['connectivity'][start_ptr:end_ptr]
    
    # Check for invalid indices
    invalid = [pid for pid in point_indices if pid < 0 or pid >= data['num_points']]
    if invalid:
      indices_str = f"{point_indices} ⚠ INVALID: {invalid}"
    else:
      indices_str = str(point_indices)
    
    print(f"{c_idx:<6} | {indices_str:<30} | {data['types'][c_idx]:<6}")
  
  print(f"{'='*70}\n")

"""
## Verification against the analytical fields set by the tests

The `.c` tests fill their fields with known analytical functions, so the written
file can be *verified*, not merely inspected. These routines print a compact,
deterministic report suitable for use as a Basilisk reference (`.ref`) file:
every quantity is either an integer, a rounded float, or a thresholded verdict,
so the output is stable across builds and optimisation levels.
"""

# VTK cell types we expect from the output routines
VTK_TYPE_NAMES = {3: 'LINE', 4: 'POLY_LINE', 7: 'POLYGON',
                  8: 'PIXEL', 9: 'QUAD', 11: 'VOXEL', 12: 'HEXAHEDRON'}

# struct format and item size for each VTK DataArray type
VTK_DTYPES = {'Float64': ('d', 8), 'Float32': ('f', 4), 'Int64': ('q', 8),
              'Int32': ('i', 4), 'Int8': ('b', 1), 'UInt8': ('B', 1)}


def read_vtu_arrays(filename):
  """Read every DataArray of a VTU file, driven by the offsets declared in its
  own XML header.

  The output routines do not emit the ``CellData``/``Points``/``Cells``
  sections in a fixed order, so nothing here may depend on that order: each
  array is located solely through its ``offset`` attribute.

  Returns a dict with 'num_points', 'num_cells', 'points' (list of xyz
  tuples), 'offsets', 'types', 'connectivity' and 'cell_data' (a dict mapping
  field name to a flat list of values).
  """
  import re

  with open(filename, 'rb') as f:
    content = f.read()

  xml_part = content.split(b'<AppendedData')[0].decode('utf-8')
  num_points = int(re.search(r'NumberOfPoints="(\d+)"', xml_part).group(1))
  num_cells = int(re.search(r'NumberOfCells="(\d+)"', xml_part).group(1))

  data_start = content.find(b'_', content.find(b'<AppendedData')) + 1

  def read_block(offset, dtype):
    """Read one appended block: a UInt64 byte count followed by its payload."""
    fmt, size = VTK_DTYPES[dtype]
    start = data_start + offset
    nbytes = struct.unpack('<Q', content[start:start + 8])[0]
    payload = content[start + 8:start + 8 + nbytes]
    return list(struct.unpack(f'<{nbytes // size}{fmt}', payload))

  # Walk every DataArray declared in the header, remembering which section it
  # belongs to so that Points and CellData with no Name are told apart.
  arrays, cell_data = {}, {}
  section = None
  for line in xml_part.splitlines():
    stripped = line.strip()
    for tag in ('Points', 'Cells', 'CellData', 'PointData'):
      if stripped.startswith(f'<{tag}'):
        section = tag
    if not stripped.startswith('<DataArray'):
      continue
    dtype = re.search(r'type="(\w+)"', stripped).group(1)
    offset = int(re.search(r'offset="(\d+)"', stripped).group(1))
    name_match = re.search(r'Name="([^"]+)"', stripped)
    ncomp_match = re.search(r'NumberOfComponents="(\d+)"', stripped)
    values = read_block(offset, dtype)
    if section == 'Points':
      arrays['points_flat'] = values
    elif section == 'Cells':
      arrays[name_match.group(1)] = values
    elif section == 'CellData':
      name = name_match.group(1)
      ncomp = int(ncomp_match.group(1)) if ncomp_match else 1
      cell_data[name] = values
      cell_data[name + '_ncomp'] = ncomp

  flat = arrays['points_flat']
  points = [tuple(flat[i * 3:i * 3 + 3]) for i in range(num_points)]

  return {
    'num_points': num_points,
    'num_cells': num_cells,
    'points': points,
    'offsets': arrays['offsets'],
    'types': arrays['types'],
    'connectivity': arrays['connectivity'],
    'cell_data': cell_data,
  }


def cell_vertices(data, c_idx):
  """Return the list of vertex coordinates of cell `c_idx`."""
  start = 0 if c_idx == 0 else data['offsets'][c_idx - 1]
  end = data['offsets'][c_idx]
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
  consistent = data['offsets'][-1] == len(data['connectivity'])
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
      # A quad is flat in one direction: in a 2D run that is z, in a
      # slice file it is the slice normal. Express the expected ordering
      # in the two remaining axes, so every quad is checked rather than
      # skipped. The convention is counter-clockwise from the low corner.
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
  """Verify the analytical fields written by `test_output_vtu[_box].c`.

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
    ncomp = cd.get('u.x_ncomp', 3)
    max_err = 0.0
    for c in range(data['num_cells']):
      verts = cell_vertices(data, c)
      centroid = [sum(v[d] for v in verts) / len(verts) for d in range(3)]
      written = cd['u.x'][c * ncomp:(c + 1) * ncomp]
      # All three components are compared, in 2D as well as in 3D: a 2D
      # run writes u.z as 0 and its cells are flat in z, so the centroid
      # is 0 too and the comparison stays meaningful.
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
      # Use the largest extent rather than the x one: cells of a slice
      # file are flat in the slice-normal direction, whose extent is 0.
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
    outside = 0
    for c in range(data['num_cells']):
      verts = cell_vertices(data, c)
      centroid = [sum(v[d] for v in verts) / len(verts) for d in range(3)]
      ndim = len(box[0])
      if any(not (box[0][d] <= centroid[d] <= box[1][d]) for d in range(ndim)):
        outside += 1
    print(f"cells inside box       "
          f"{'OK' if outside == 0 else f'FAIL ({outside} outside)'}")
    failures += outside > 0
  return failures


def check_cardioid(data, a=0.15):
  """Verify the interface written by `test_output_vtu_facets.c`.

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

  if 'kappa' in data['cell_data']:
    kappa = data['cell_data']['kappa']
    n_finite = sum(1 for k in kappa if math.isfinite(k))
    print(f"kappa finite           "
          f"{'OK' if n_finite == len(kappa) else f'FAIL ({n_finite}/{len(kappa)})'}")
    failures += n_finite != len(kappa)
  return failures


def check_plane(data, axis, value):
  """Verify that a slice file is genuinely planar.

  `output_slice_vtu()` is asked for the plane `axis = value`, so every vertex
  it writes must lie exactly on it, and the cells must be flat in that
  direction. This is the property that distinguishes a slice from a full
  domain dump, so it is checked separately.

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

  data = read_vtu_arrays(filename)
  failures = check_topology(data)
  if plane is not None:
    failures += check_plane(data, plane[0], plane[1])
  if cardioid:
    failures += check_cardioid(data)
  else:
    failures += check_analytical_fields(data, box=box)

  print(f"verdict                {'PASS' if failures == 0 else f'FAIL ({failures})'}")
  return 0 if failures == 0 else 1


if __name__ == "__main__":
  args = [a for a in sys.argv[1:] if not a.startswith('--')]
  flags = {a for a in sys.argv[1:] if a.startswith('--')}

  if not args:
    print("Usage: python3 test_output_vtu.py [--check] [--box=X0,Y0,X1,Y1] "
          "[--plane=AXIS=VALUE] [--cardioid] <file1.vtu> ...")
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
      print_vtu_info(filename)
    else:
      print(f"Error: {filename} not found")
