#!/usr/bin/env python3
"""
Checkers for the tests here, one subcommand each. The norms are computed in C;
what happens here is the thresholds and the structural checks.

  python3 test_spectra.py sample spectra_sample.asc
  python3 test_spectra.py modes  spectra_modes.asc
  python3 test_spectra.py amr    spectra_amr.asc
  python3 test_spectra.py ascii  spec.asc spec_u.asc
  python3 test_spectra.py hdf5   spectra_hdf5.asc

Each prints one line per case then PASS or FAIL, and exits non-zero on
failure. The constants mirror the .c sources and are kept in step by hand.
"""
import math
import re
import sys

L0 = 2.*3.14159265358979   # Basilisk's pi is a 15-digit macro, not math.pi
Z0 = -L0/2.
M = 32                     # the base grid every test starts from
TOL = 1e-12                # a value the lattice represents exactly: roundoff
LEAK_MIN = 1e-4            # a case sampled off the cell centres must show it


def nshells(m):
  """Mirrors nshells() in spectra_shell.h: the corner of the k plane."""
  return int(math.sqrt(2.*(m/2.)**2) + 0.5) + 1


def verdict(failures, message):
  if failures:
    print('FAIL: ' + '; '.join(failures))
    return 1
  print('PASS: ' + message)
  return 0


def rows_of(filename):
  """The non-comment, non-blank lines of a table, split into fields."""
  with open(filename) as f:
    for line in f:
      line = line.strip()
      if line and not line.startswith('#'):
        yield line.split()


# ---------------------------------------------------------------- sample ---
#
# `tag m holes err_c dz`. Both errors measure the distance between a lattice
# point and the cell centre that answered it, so a half-cell step bounds them;
# |d cos/dx| <= 1 carries that bound to `err_c`. `#if TREE` skips the tree
# cases on multigrid, so only the rows present are checked.

COARSE_MIN = 1e-3          # the coarse case must really be displaced
SAMPLE_REQUIRED = ('uniform', 'face')


def check_sample(filename):
  rows, base = {}, M
  with open(filename) as f:
    for line in f:
      w = line.split()
      if not w:
        continue
      if w[0] == '#':
        if len(w) > 4 and w[1] == 'L0':
          base = int(w[4])
        continue
      rows[w[0]] = (int(w[1]), int(w[2]), float(w[3]), float(w[4]))

  failures = [f'missing rows: {t}' for t in SAMPLE_REQUIRED if t not in rows]
  half = 0.5*L0/base
  for tag, (m, holes, err_c, dz) in sorted(rows.items()):
    print(f'{tag:<8} m {m:4d}  holes {holes}  err_c {err_c:.4g}  dz {dz:.4g}')
    if holes:
      failures.append(f'{tag}: {holes} unfilled lattice points')

  # a lattice point on a cell centre is answered exactly
  for tag in ('uniform', 'matched'):
    if tag in rows and (rows[tag][2] > TOL or rows[tag][3] > TOL):
      failures.append(f'{tag}: not exact (err_c {rows[tag][2]:.3g}, '
                      f'dz {rows[tag][3]:.3g})')

  # z = 0 is a face: the plane comes back exactly half a cell away
  if 'face' in rows and abs(rows['face'][3] - half) > TOL:
    failures.append(f'face: dz {rows["face"][3]:.6g}, expected the half cell '
                    f'{half:.6g}')

  # a coarser lattice is displaced by up to half a cell, and really is
  if 'coarse' in rows and not COARSE_MIN < rows['coarse'][2] <= half + TOL:
    failures.append(f'coarse: err_c {rows["coarse"][2]:.3g} outside '
                    f'({COARSE_MIN:g}, {half:.4g}]')

  # refining under the same lattice shrinks that displacement
  if 'coarse' in rows and 'refined' in rows:
    if not rows['refined'][2] < 0.5*rows['coarse'][2]:
      failures.append(f'refined: err_c {rows["refined"][2]:.3g} did not shrink '
                      f'below half of the coarse {rows["coarse"][2]:.3g}')

  return verdict(failures,
                 'sample_scalar_plane matches the cell-centre displacement.')


# ----------------------------------------------------------------- modes ---
#
# `tag kind m holes bpeak bexp E Eexp leak parseval sum`; `leak` and
# `parseval` are already residuals. `exact` rows use modes the lattice
# represents exactly, so their tolerance is a machine tolerance. A `leaky` row
# is sampled off the cell centres on purpose and must lose measurable energy
# -- only checking that it stays small would let the rule stop mattering.
# Parseval belongs to the transform, so it is held to TOL on every row.

def check_modes(filename):
  failures, seen = [], 0
  for w in rows_of(filename):
    tag, kind = w[0], w[1]
    m, holes, bpeak, bexp = (int(v) for v in w[2:6])
    E, Eexp, leak, parseval, total = (float(v) for v in w[6:11])
    seen += 1

    print(f'{tag:<11}{kind:<7} m {m:4d}  bin {bpeak:3d} (exp {bexp:3d})  '
          f'E {E:.12f} (exp {Eexp:.4f})  leak {leak:.2e}  '
          f'Parseval {parseval:.2e}  sum {total:.12f}')

    if holes:
      failures.append(f'{tag}: {holes} unfilled lattice points')
    if bpeak != bexp:
      failures.append(f'{tag}: peak in bin {bpeak}, expected {bexp}')
    if parseval > TOL:
      failures.append(f'{tag}: Parseval residual {parseval:.3g}')

    if kind == 'exact':
      if abs(E - Eexp) > TOL:
        failures.append(f'{tag}: E {E:.12g}, expected {Eexp:g}')
      if leak > TOL:
        failures.append(f'{tag}: {leak:.3g} of energy outside its bins')
    elif kind == 'leaky':
      if leak < LEAK_MIN:
        failures.append(f'{tag}: leak {leak:.3g} below {LEAK_MIN:g}, so '
                        f'sampling off the cell centres no longer shows')
    else:
      failures.append(f'{tag}: unknown kind {kind!r}')

  if not seen:
    failures.append('no rows read')
  return verdict(failures,
                 'single modes land in the right shell, and Parseval holds.')


# ------------------------------------------------------------------- amr ---
#
# `tag kind m holes lmin lmax E leak parseval`. Both rows carry the same mode,
# so E = 1/2 and no leak; what separates them is where the plane sits. The
# `exact` row lies at one level and must reproduce both; the `leaky` row
# straddles two and must lose measurable energy, or the `exact` row is no
# longer evidence for anything. `lmin`/`lmax` are checked too, so a refinement
# change that flattens the contrast fails rather than passing vacuously.

AMR_EXPECTED = 0.5


def check_amr(filename):
  failures, seen = [], 0
  for w in rows_of(filename):
    tag, kind = w[0], w[1]
    m, holes, lmin, lmax = (int(v) for v in w[2:6])
    E, leak, parseval = (float(v) for v in w[6:9])
    seen += 1

    print(f'{tag:<11}{kind:<7} m {m:4d}  levels {lmin}-{lmax}  '
          f'E {E:.12f}  leak {leak:.2e}  Parseval {parseval:.2e}')

    if holes:
      failures.append(f'{tag}: {holes} unfilled lattice points')
    if parseval > TOL:
      failures.append(f'{tag}: Parseval residual {parseval:.3g}')

    if kind == 'exact':
      if lmin != lmax:
        failures.append(f'{tag}: levels {lmin}-{lmax}, expected one level')
      if abs(E - AMR_EXPECTED) > TOL:
        failures.append(f'{tag}: E {E:.12g}, expected {AMR_EXPECTED}')
      if leak > TOL:
        failures.append(f'{tag}: {leak:.3g} of energy outside the bin')
    elif kind == 'leaky':
      if lmin == lmax:
        failures.append(f'{tag}: level {lmin} throughout, so the plane no '
                        f'longer straddles a refinement boundary')
      if leak < LEAK_MIN:
        failures.append(f'{tag}: leak {leak:.3g} below {LEAK_MIN:g}, so '
                        f'straddling the boundary no longer costs anything')
    else:
      failures.append(f'{tag}: unknown kind {kind!r}')

  if seen != 2:
    failures.append(f'{seen} rows, expected 2')
  return verdict(failures,
                 'the plane is exact at one level and leaks across two.')


# ----------------------------------------------------------------- ascii ---
#
# Structure rather than transform: block count, header metadata, row ordering,
# wavenumber columns, plane heights, and that each value is where it belongs.
# The fields are scaled by z, so a plane's expected spectrum follows from the
# height beside it -- which is what makes a transposed stride visible. That
# also sets the tolerance: z is written to nine digits and E derived from it,
# so energies are compared relatively at REL; a misplaced value is wrong by
# order one, not by a rounding.

NZ = 4
REL = 1e-7      # z is written to nine digits, and E is derived from it
META = 1e-5     # the block header writes L0 with %g

# (t, field names, {name: (bin, amplitude(z))}) per block of spec.asc.
BLOCKS = [
  (0.00, ['a', 'b'], {'a': (0, lambda z: (3. + z)**2),
                      'b': (5, lambda z: 0.5*(1. + z)**2)}),
  (0.25, ['a', 'b'], {'a': (0, lambda z: (5. + z)**2),
                      'b': (7, lambda z: 0.5*(2. + z)**2)}),
]
VECTOR = (0.50, ['v.x', 'v.y', 'v.z'], {'v.x': (5, lambda z: 0.5),
                                        'v.y': (0, lambda z: 4.),
                                        'v.z': (3, lambda z: 0.5)})

HEADER = re.compile(r'# Spectrum: t = (\S+), L0 = (\S+), m = (\d+), nz = (\d+),'
                    r' nk = (\d+), hmin = (\S+), hmax = (\S+)')


def read_blocks(filename):
  """Returns [(meta, field names, rows)], rows being lists of floats."""
  blocks, meta, names, rows = [], None, None, []
  with open(filename) as f:
    for line in f:
      line = line.rstrip('\n')
      m = HEADER.match(line)
      if m:
        if meta:
          blocks.append((meta, names, rows))
        t, l0, mm, nz, nk, hmin, hmax = m.groups()
        meta = dict(t=float(t), L0=float(l0), m=int(mm), nz=int(nz),
                    nk=int(nk), hmin=float(hmin), hmax=float(hmax))
        names, rows = None, []
      elif line.startswith('# [0]iz'):
        names = re.findall(r'\[\d+\]E\((.*?)\)', line)
      elif line.strip():
        rows.append([float(v) for v in line.split()])
  if meta:
    blocks.append((meta, names, rows))
  return blocks


def check_block(tag, meta, names, rows, expected, failures):
  t, exp_names, spec = expected
  nz, nk = meta['nz'], meta['nk']

  if abs(meta['t'] - t) > TOL:
    failures.append(f'{tag}: t = {meta["t"]}, expected {t}')
  if meta['m'] != M or abs(meta['L0'] - L0) > META*L0:
    failures.append(f'{tag}: m = {meta["m"]}, L0 = {meta["L0"]}')
  if nk != nshells(M):
    failures.append(f'{tag}: nk = {nk}, expected {nshells(M)}')
  if names != exp_names:
    failures.append(f'{tag}: columns {names}, expected {exp_names}')
    return
  if len(rows) != nz*nk:
    failures.append(f'{tag}: {len(rows)} rows, expected {nz}*{nk}')
    return
  if any(len(r) != 4 + len(names) for r in rows):
    failures.append(f'{tag}: rows are not {4 + len(names)} columns wide')
    return

  cell = L0/M
  heights = []
  for iz in range(nz):
    plane = rows[iz*nk:(iz + 1)*nk]

    # rows are ordered by plane and then by bin, and z is constant on a plane
    if any(int(r[0]) != iz for r in plane):
      failures.append(f'{tag}: plane {iz} rows do not all carry iz = {iz}')
    if any(int(r[2]) != b for b, r in enumerate(plane)):
      failures.append(f'{tag}: plane {iz} bins are not 0..{nk - 1} in order')
    if any(abs(r[3] - 2.*math.pi*r[2]/L0) > TOL for r in plane):
      failures.append(f'{tag}: plane {iz} kphys is not 2*pi*k/L0')

    z = plane[0][1]
    if any(r[1] != z for r in plane):
      failures.append(f'{tag}: plane {iz} z is not constant')
    heights.append(z)

    # snap_to_cell() puts every height on a cell centre inside [hmin, hmax]
    if abs((z - Z0)/cell - 0.5 - round((z - Z0)/cell - 0.5)) > 1e-6:
      failures.append(f'{tag}: z = {z:.6g} is not a cell centre')
    if not meta['hmin'] - cell <= z <= meta['hmax'] + cell:
      failures.append(f'{tag}: z = {z:.6g} outside '
                      f'[{meta["hmin"]:.6g}, {meta["hmax"]:.6g}]')

    # the expected spectrum follows from this plane's own height
    for k, name in enumerate(names):
      bexp, amp = spec[name]
      col = [r[4 + k] for r in plane]
      total, Eexp = sum(col), amp(z)
      if abs(col[bexp] - Eexp) > REL*Eexp:
        failures.append(f'{tag}: plane {iz} {name} bin {bexp} = '
                        f'{col[bexp]:.12g}, expected {Eexp:.12g}')
      if abs(total - Eexp) > REL*Eexp:
        failures.append(f'{tag}: plane {iz} {name} has {total - Eexp:.3g} '
                        f'outside bin {bexp}')

  if heights != sorted(heights) or len(set(heights)) != nz:
    failures.append(f'{tag}: heights {heights} are not increasing and distinct')
  print(f'{tag}: t {meta["t"]:.2f}  nz {nz}  nk {nk}  fields {",".join(names)}'
        f'  z [{heights[0]:.5f}, {heights[-1]:.5f}]')


def check_ascii(scalar_file, vector_file):
  failures = []
  blocks = read_blocks(scalar_file)
  if len(blocks) != len(BLOCKS):
    failures.append(f'{len(blocks)} blocks in {scalar_file}, '
                    f'expected {len(BLOCKS)}')
  for i, (block, expected) in enumerate(zip(blocks, BLOCKS)):
    check_block(f'block {i}', *block, expected, failures)

  vblocks = read_blocks(vector_file)
  if len(vblocks) != 1:
    failures.append(f'{len(vblocks)} blocks in {vector_file}, expected 1')
  else:
    meta, names, rows = vblocks[0]
    if meta['nz'] != 1:
      failures.append(f'vector: nz = {meta["nz"]}, expected 1')
    check_block('vector', meta, names, rows, VECTOR, failures)

  return verdict(failures, 'the ASCII blocks are laid out as documented.')


# ------------------------------------------------------------------ hdf5 ---
#
# The summary test_spectra_hdf5.c reads back out of the stored file:
#
#   attr <L0> <m> <nz> <nk>
#   shape <name> <dims...>
#   t <it> <value>
#   kphys <max residual against 2*pi*k/L0>
#   peak <field> <it> <iz> <bpeak> <Epeak> <sum> <z>
#
# The twelve spectra are all distinct -- `a`'s mode moves with the block, both
# amplitudes are scaled by z -- so the expected bin and energy follow from the
# block index and the stored height, and axes written in the wrong order fail
# on the values rather than only on the shape.

NT = 3
NK = nshells(M)
HDF5_REL = 1e-9

# the bin of field <f> in block <it>, and its energy at height z
HDF5_SPEC = {
  'a': (lambda it: 5 + it, lambda it, z: 0.5*(1. + z)**2),
  'b': (lambda it: 0,      lambda it, z: (3. + it + z)**2),
}


def check_hdf5(filename):
  shapes, times, peaks = {}, {}, []
  attr = dk = None
  for w in rows_of(filename):
    if w[0] == 'attr':
      attr = (float(w[1]), int(w[2]), int(w[3]), int(w[4]))
    elif w[0] == 'shape':
      shapes[w[1]] = tuple(int(v) for v in w[2:])
    elif w[0] == 't':
      times[int(w[1])] = float(w[2])
    elif w[0] == 'kphys':
      dk = float(w[1])
    elif w[0] == 'peak':
      peaks.append((w[1], int(w[2]), int(w[3]), int(w[4]),
                    float(w[5]), float(w[6]), float(w[7])))

  failures = []
  if attr != (L0, M, NZ, NK):
    failures.append(f'attributes {attr}, expected {(L0, M, NZ, NK)}')

  expected_shapes = {'t': (NT,), 'z': (NT, NZ), 'kphys': (NK,),
                     'E': (NT, NZ, NK), 'v.x': (NT, NZ, NK),
                     'v.y': (NT, NZ, NK), 'v.z': (NT, NZ, NK)}
  for name, exp in expected_shapes.items():
    if shapes.get(name) != exp:
      failures.append(f'/{name} has shape {shapes.get(name)}, expected {exp}')

  for it in range(NT):
    if abs(times.get(it, -1) - 0.1*it) > TOL:
      failures.append(f't[{it}] = {times.get(it)}, expected {0.1*it}')
  if dk is None or dk > TOL:
    failures.append(f'/kphys residual {dk}')
  if len(peaks) != 2*NT*NZ:
    failures.append(f'{len(peaks)} spectra summarised, expected {2*NT*NZ}')

  heights = {}
  for name, it, iz, bpeak, peak, total, z in peaks:
    bexp, Eexp = HDF5_SPEC[name][0](it), HDF5_SPEC[name][1](it, z)
    heights.setdefault(it, []).append(z)
    if bpeak != bexp:
      failures.append(f'{name}[{it}][{iz}]: peak in bin {bpeak}, '
                      f'expected {bexp}')
    if abs(peak - Eexp) > HDF5_REL*Eexp:
      failures.append(f'{name}[{it}][{iz}]: E {peak:.12g}, '
                      f'expected {Eexp:.12g}')
    if abs(total - Eexp) > HDF5_REL*Eexp:
      failures.append(f'{name}[{it}][{iz}]: {total - Eexp:.3g} outside '
                      f'bin {bexp}')

  # the zone widens between blocks, so the stored heights spread with it
  spans = [max(v) - min(v) for it, v in sorted(heights.items())]
  if spans != sorted(spans) or spans[0] == spans[-1]:
    failures.append(f'plane heights did not spread with the zone: {spans}')

  print(f'nt {NT}  nz {NZ}  nk {NK}  {len(peaks)} spectra  '
        f'z spans {", ".join(f"{s:.4f}" for s in spans)}')
  return verdict(failures,
                 'the HDF5 file stores the blocks on the axes it documents.')


CHECKS = {'sample': (check_sample, 1), 'modes': (check_modes, 1),
          'amr': (check_amr, 1), 'ascii': (check_ascii, 2),
          'hdf5': (check_hdf5, 1)}


def main():
  argv = sys.argv[1:]
  if not argv or argv[0] not in CHECKS or len(argv) - 1 != CHECKS[argv[0]][1]:
    print(__doc__.strip())
    return 1
  check, _ = CHECKS[argv[0]]
  return check(*argv[1:])


if __name__ == '__main__':
  sys.exit(main())
