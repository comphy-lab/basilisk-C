#!/usr/bin/env python3
"""
Checkers for the tests here, one subcommand each. The norms are computed in C;
what happens here is the thresholds and the structural checks.

  python3 test_spectra.py sample         spectra_sample.asc
  python3 test_spectra.py modes          spectra_modes.asc
  python3 test_spectra.py cross          spectra_cross.asc
  python3 test_spectra.py foliate_sample spectra_foliate_sample.asc
  python3 test_spectra.py foliate        spectra_foliate.asc spectra_foliate_cross.asc
  python3 test_spectra.py amr            spectra_amr.asc
  python3 test_spectra.py ascii          spec.asc spec_u.asc
  python3 test_spectra.py hdf5           spectra_hdf5.asc
  python3 test_spectra.py restrict_plane spectra_restrict_plane.asc
  python3 test_spectra.py restrict_rr    spectra_restrict_rr.asc

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


# ----------------------------------------------------------------- cross ---
#
# `tag kind m holes bexp E Eexp leak sum`. `exact` rows pin one bin the way
# `modes` does; `zero` rows check the *total* cross power, since two fields on
# disjoint Fourier modes must cross to zero everywhere, not just off the bin
# a peak-finder would have picked.

def check_cross(filename):
  failures, seen = [], 0
  for w in rows_of(filename):
    tag, kind = w[0], w[1]
    m, holes, bexp = (int(v) for v in w[2:5])
    E, Eexp, leak, total = (float(v) for v in w[5:9])
    seen += 1

    print(f'{tag:<11}{kind:<7} m {m:4d}  bin {bexp:3d}  E {E:.12f} '
          f'(exp {Eexp:.4f})  leak {leak:.2e}  sum {total:.12f}')

    if holes:
      failures.append(f'{tag}: {holes} unfilled lattice points')

    if kind == 'exact':
      if abs(E - Eexp) > TOL:
        failures.append(f'{tag}: E {E:.12g}, expected {Eexp:g}')
      if leak > TOL:
        failures.append(f'{tag}: {leak:.3g} of cross power outside bin {bexp}')
    elif kind == 'zero':
      if abs(total) > TOL:
        failures.append(f'{tag}: total cross power {total:.3g}, expected 0 '
                        f'(orthogonal modes)')
    else:
      failures.append(f'{tag}: unknown kind {kind!r}')

  if not seen:
    failures.append('no rows read')
  return verdict(failures,
                 'cross_shell_average matches the auto-spectrum on identical '
                 'fields and vanishes on orthogonal modes.')


# --------------------------------------------------------- foliate_sample ---
#
# `tag m nz holes err`. Each row's expected sum is analytic, exact to
# roundoff since every case lands on cell centres -- so all three rows share
# `sample`'s TOL rather than needing a case-specific tolerance.

FOLIATE_REQUIRED = ('zero_mean', 'offset_mean', 'wrong_mean')


def check_foliate_sample(filename):
  failures, seen = [], set()
  for w in rows_of(filename):
    tag = w[0]
    m, nz, holes = (int(v) for v in w[1:4])
    err = float(w[4])
    seen.add(tag)

    print(f'{tag:<12} m {m:4d}  nz {nz}  holes {holes}  err {err:.3g}')

    if holes:
      failures.append(f'{tag}: {holes} unfilled lattice points')
    if err > TOL:
      failures.append(f'{tag}: err {err:.3g}, expected <= {TOL:g}')

  failures += [f'missing row: {t}' for t in FOLIATE_REQUIRED if t not in seen]
  return verdict(failures,
                 'sample_scalar_stack_sum foliates to the analytic anomaly '
                 'sum, linear in the supplied per-slab means.')


# --------------------------------------------------------------- foliate ---
#
# `spectrum_scalar_foliated()`/`cross_spectrum_scalar_foliated()` reuse the
# writer and the transform, both already checked by `ascii` and `modes`; this
# checks only the foliation and its plumbing through them, against three
# blocks (`a`, `b`, `c`) in one ASCII file and one row in a separate cross
# file. `read_blocks()` is shared with `ascii`, defined further down.

FOLIATE_NZ = 4
FOLIATE_TAGS = ('a', 'b', 'c')
# tag -> {bin: expected E}; bins not listed must carry (near) no power
FOLIATE_EXPECTED = {'a': {}, 'b': {5: 0.5*FOLIATE_NZ**2},
                    'c': {0: float (FOLIATE_NZ**2), 5: 0.5*FOLIATE_NZ**2}}


def check_foliate(ascii_file, cross_file):
  failures = []
  blocks = read_blocks(ascii_file)
  if len(blocks) != len(FOLIATE_TAGS):
    failures.append(f'{len(blocks)} blocks in {ascii_file}, '
                    f'expected {len(FOLIATE_TAGS)}')

  for tag, (meta, names, rows) in zip(FOLIATE_TAGS, blocks):
    exp = FOLIATE_EXPECTED[tag]
    if meta['nz'] != 1:
      failures.append(f'{tag}: nz = {meta["nz"]}, expected 1 (one folded block)')
    if names != [tag]:
      failures.append(f'{tag}: columns {names}, expected [{tag!r}]')
      continue
    if len(rows) != meta['nk']:
      failures.append(f'{tag}: {len(rows)} rows, expected {meta["nk"]}')
      continue

    zmid = 0.5*(meta['hmin'] + meta['hmax'])
    z = rows[0][1]
    if abs(z - zmid) > TOL:
      failures.append(f'{tag}: stored z {z:.6g}, expected the midpoint '
                      f'{zmid:.6g}')

    E = {int(r[2]): r[4] for r in rows}
    total = sum(r[4] for r in rows)
    for b, Eexp in exp.items():
      if abs(E.get(b, 0.) - Eexp) > TOL:
        failures.append(f'{tag}: bin {b} = {E.get(b, 0.):.12g}, '
                        f'expected {Eexp:.12g}')
    leak = total - sum(exp.values())
    if abs(leak) > TOL:
      failures.append(f'{tag}: {leak:.3g} of power outside '
                      f'{sorted(exp) or "no bins (should be all-zero)"}')

    print(f'{tag}: nz {meta["nz"]}  z {z:.4f}  bins {sorted(exp) or "none"}  '
          f'total {total:.12f}')

  seen = 0
  for w in rows_of(cross_file):
    holes, bexp = int(w[0]), int(w[1])
    E, Eexp = float(w[2]), float(w[3])
    seen += 1
    print(f'cross_self: bin {bexp}  E {E:.12f} (exp {Eexp:.4f})  holes {holes}')
    if holes:
      failures.append(f'cross_self: {holes} unfilled lattice points')
    if abs(E - Eexp) > TOL:
      failures.append(f'cross_self: E {E:.12g}, expected {Eexp:g}')
  if not seen:
    failures.append(f'no rows read from {cross_file}')

  return verdict(failures,
                 'the foliated sum, its spectrum and its cross-spectrum all '
                 'match the analytic anomaly.')


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


# --------------------------------------------------------- restrict_plane ---
#
# `tag holes err_hvar err_zlin err_smooth err_curved err_radial dmin dmax
# filled`. `hvar` is a step at $x = 0$: no lattice, restricted or not, lands
# a cell centre exactly on it, so both strategies carry a real, unavoidable
# error there -- what matters is that `restrict` does not do *worse* than
# `region`. The other four fields are smooth enough that the coarse-cell
# average should recover them close to exactly, so `restrict` is held to a
# tight roundoff-scale tolerance on those, while `region` must show a real,
# nonzero bias (the vacuity guard: without it a passing `restrict` row would
# prove nothing).

RESTRICT_PLANE_REQUIRED = ('region', 'restrict')
RESTRICT_BIAS_MIN = 1e-3   # region's error on the smooth fields must be real
RESTRICT_SMOOTH_TOL = 2e-2 # restrict's residual on a coarser lattice, not 0


def check_restrict_plane(filename):
  failures, rows = [], {}
  for w in rows_of(filename):
    tag = w[0]
    holes = int(w[1])
    err = [float(v) for v in w[2:7]]
    dmin, dmax = float(w[7]), float(w[8])
    filled = int(w[9])
    rows[tag] = (holes, err, dmin, dmax, filled)
    print(f'{tag:<10} holes {holes:4d}  err {" ".join(f"{e:.2e}" for e in err)}'
          f'  D [{dmin:.4g}, {dmax:.4g}]  filled {filled}')
    if holes:
      failures.append(f'{tag}: {holes} unfilled lattice points')

  failures += [f'missing row: {t}' for t in RESTRICT_PLANE_REQUIRED
              if t not in rows]

  if 'region' in rows:
    holes, err, dmin, dmax, filled = rows['region']
    if dmin == dmax:
      failures.append('region: cells all the same size -- grid does not '
                      'straddle the refined band, test is vacuous')
    # smooth, curved, radial (indices 2-4): must show a real bias
    if min(err[2:5]) < RESTRICT_BIAS_MIN:
      failures.append(f'region: smallest smooth/curved/radial error '
                      f'{min(err[2:5]):.3g} below {RESTRICT_BIAS_MIN:g}, '
                      f'unrestricted sampling no longer shows a bias')

  if 'restrict' in rows:
    holes, err, dmin, dmax, filled = rows['restrict']
    if dmin != dmax:
      failures.append(f'restrict: cell sizes not uniform ({dmin} vs {dmax})')
    if max(err[2:5]) > RESTRICT_SMOOTH_TOL:
      failures.append(f'restrict: max smooth/curved/radial error '
                      f'{max(err[2:5]):.3g} above {RESTRICT_SMOOTH_TOL:g}')

  # hvar (index 0) and zlin (index 1): restrict must not be worse than region
  if 'region' in rows and 'restrict' in rows:
    for name, k in (('hvar', 0), ('zlin', 1)):
      rg, rs = rows['region'][1][k], rows['restrict'][1][k]
      if rs > rg + TOL:
        failures.append(f'{name}: restrict error {rs:.3g} worse than '
                        f'region {rg:.3g}')

  return verdict(failures,
                 'sample_scalar_plane_restrict() is uniform and at least as '
                 'accurate as the locate()-based sampler, which is biased.')


# ------------------------------------------------------------ restrict_rr ---
#
# `tag npe nz len m maxdiff`. Round-robin correctness is an equality, not a
# tolerance: every rank does the same FFT + shell-average work either way, so
# `spectrum_scalar_stack()`/`cross_spectrum_stack()` splitting it round-robin
# across ranks instead of stacking it on `pid() == 0` must reproduce the
# serial `E` to roundoff, at whatever `npe()` the test happened to build with.

RESTRICT_RR_REQUIRED = ('scalar_stack', 'cross_stack')


def check_restrict_rr(filename):
  failures, seen = [], set()
  for w in rows_of(filename):
    tag, npe_, nz, ln, m = w[0], int(w[1]), int(w[2]), int(w[3]), int(w[4])
    maxdiff = float(w[5])
    seen.add(tag)
    print(f'{tag:<14} npe {npe_}  nz {nz}  len {ln}  m {m}  '
          f'max|E_rr - E_ref| {maxdiff:.3e}')
    if maxdiff > TOL:
      failures.append(f'{tag}: max|E_rr - E_ref| {maxdiff:.3g} above {TOL:g}')

  failures += [f'missing row: {t}' for t in RESTRICT_RR_REQUIRED
              if t not in seen]
  return verdict(failures,
                 'round-robin FFT distribution matches the serial reference '
                 'exactly, independent of npe().')


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
#   format <v>
#   block <it> <group> <nz> <nk> <m> <t> <hmin> <hmax>
#   shape <it> <name> <dims...>
#   kphys <it> <nk> <max residual against 2*pi*k/L0>
#   peak <field> <it> <iz> <bpeak> <Epeak> <sum> <z>
#
# The spectra are all distinct -- `a`'s mode moves with the block, both
# amplitudes are scaled by z -- so the expected bin and energy follow from the
# block index and the stored height, and axes written in the wrong order fail
# on the values rather than only on the shape.
#
# One self-contained group per output, so the plane count differs per block
# (HDF5_NZS) with nothing shared for a later block to break. `it` is the
# position HDF5 listed the group at, ordered by name, so checking t against
# the block index is also checking that the zero-padded name sorts by time.

NT = 3
NK = nshells(M)
HDF5_REL = 1e-9
HDF5_FORMAT = 3
HDF5_NZS = (4, 2, 7)              # mirrors NZS in test_spectra_hdf5.c
HDF5_LZ = lambda it: 0.5 + 0.25*it   # the widening zone, as in the .c
HDF5_GROUP = lambda it: 't%09.4f' % (0.1*it)
HDF5_COMPS = ('v.x', 'v.y', 'v.z')

# the bin of field <f> in block <it>, and its energy at height z
HDF5_SPEC = {
  'a': (lambda it: 5 + it, lambda it, z: 0.5*(1. + z)**2),
  'b': (lambda it: 0,      lambda it, z: (3. + it + z)**2),
}


def check_hdf5(filename):
  shapes, blocks, kphys, peaks = {}, {}, {}, []
  fmt = None
  for w in rows_of(filename):
    if w[0] == 'format':
      fmt = int(w[1])
    elif w[0] == 'shape':
      shapes[(int(w[1]), w[2])] = tuple(int(v) for v in w[3:])
    elif w[0] == 'block':
      blocks[int(w[1])] = (w[2], int(w[3]), int(w[4]), int(w[5]),
                           float(w[6]), float(w[7]), float(w[8]))
    elif w[0] == 'kphys':
      kphys[int(w[1])] = (int(w[2]), float(w[3]))
    elif w[0] == 'peak':
      peaks.append((w[1], int(w[2]), int(w[3]), int(w[4]),
                    float(w[5]), float(w[6]), float(w[7])))

  failures = []
  if fmt != HDF5_FORMAT:
    failures.append(f'format {fmt}, expected {HDF5_FORMAT}')
  if len(blocks) != NT:
    failures.append(f'{len(blocks)} groups, expected {NT}')

  # each group stands alone: its own name, axes, metadata and plane count.
  # `it` is the position HDF5 listed it at, ordered by name, so matching t
  # to the block index also checks the name sorts by time.
  for it in range(NT):
    got = blocks.get(it)
    if got is None:
      failures.append(f'block {it}: missing')
      continue
    name, nz, nk, m, t, hmin, hmax = got
    nzexp, Lz = HDF5_NZS[it], HDF5_LZ(it)
    if name != HDF5_GROUP(it):
      failures.append(f'block {it}: group {name!r}, '
                      f'expected {HDF5_GROUP(it)!r} at this position')
    if (nz, nk, m) != (nzexp, NK, M):
      failures.append(f'block {it}: (nz, nk, m) = {(nz, nk, m)}, '
                      f'expected {(nzexp, NK, M)}')
    if abs(t - 0.1*it) > TOL:
      failures.append(f'block {it}: t = {t}, expected {0.1*it}')
    if abs(hmin + Lz) > TOL or abs(hmax - Lz) > TOL:
      failures.append(f'block {it}: bounds ({hmin:.6g}, {hmax:.6g}), '
                      f'expected ({-Lz:.6g}, {Lz:.6g})')

    # the axes are repeated per group, so every copy has to be right
    if kphys.get(it) is None or kphys[it][0] != NK or kphys[it][1] > TOL:
      failures.append(f'block {it}: kphys {kphys.get(it)}, '
                      f'expected ({NK}, residual <= {TOL:g})')
    if shapes.get((it, 'z')) != (nzexp,):
      failures.append(f'block {it}: z has shape {shapes.get((it, "z"))}, '
                      f'expected {(nzexp,)}')
    for f in list(HDF5_SPEC) + list(HDF5_COMPS):
      if shapes.get((it, f)) != (nzexp, NK):
        failures.append(f'block {it}: {f} has shape '
                        f'{shapes.get((it, f))}, expected {(nzexp, NK)}')

  if len(peaks) != len(HDF5_SPEC)*sum(HDF5_NZS):
    failures.append(f'{len(peaks)} spectra summarised, '
                    f'expected {len(HDF5_SPEC)*sum(HDF5_NZS)}')

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

  # the zone widens between blocks and each block's planes sit inside its own
  # bounds, so a group read at the wrong width shows up in the heights too
  cell = L0/M
  for it in range(NT):
    zs = sorted(set(heights.get(it, [])))
    Lz = HDF5_LZ(it)
    if len(zs) != HDF5_NZS[it]:
      failures.append(f'block {it}: {len(zs)} distinct heights, '
                      f'expected {HDF5_NZS[it]}')
    if zs and not -Lz - cell <= zs[0] <= zs[-1] <= Lz + cell:
      failures.append(f'block {it}: heights [{zs[0]:.6g}, {zs[-1]:.6g}] '
                      f'outside [{-Lz:.6g}, {Lz:.6g}]')

  print(f'format {fmt}  nt {NT}  nz {HDF5_NZS}  nk {NK}  '
        f'{len(peaks)} spectra')
  return verdict(failures,
                 'each output is a self-contained group, whatever its plane '
                 'count.')


CHECKS = {'sample': (check_sample, 1), 'modes': (check_modes, 1),
          'cross': (check_cross, 1),
          'foliate_sample': (check_foliate_sample, 1),
          'foliate': (check_foliate, 2),
          'amr': (check_amr, 1), 'ascii': (check_ascii, 2),
          'hdf5': (check_hdf5, 1),
          'restrict_plane': (check_restrict_plane, 1),
          'restrict_rr': (check_restrict_rr, 1)}


def main():
  argv = sys.argv[1:]
  if not argv or argv[0] not in CHECKS or len(argv) - 1 != CHECKS[argv[0]][1]:
    print(__doc__.strip())
    return 1
  check, _ = CHECKS[argv[0]]
  return check(*argv[1:])


if __name__ == '__main__':
  sys.exit(main())
