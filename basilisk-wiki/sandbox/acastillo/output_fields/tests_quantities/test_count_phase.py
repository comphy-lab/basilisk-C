#!/usr/bin/env python3
"""
Checker for `count_phase()`, shared by `test_count_phase.c` (`above = 1`) and
`test_count_phase_bubbles.c` (`above = 0`). Both use the same three spheres and
the same velocity field, so both expect the same numbers.

Row layout: `i t j volume bx by bz ux uy uz area` in 3D, dropping bz and uz in
2D; the dimension is inferred from the row width.

Closed form for sphere k, with u = U0 + Omega*z_hat x r:

  volume = 4/3*pi*R^3      b = c          area = 4*pi*R^2
  <u> = (U0x - Omega*cy, U0y + Omega*cx, U0z)

the average being exact because a sphere is symmetric about its centre and u is
linear. In 2D the same centres and radii give three disjoint discs, with
volume = pi*R^2 and area = 2*pi*R (a perimeter).

Regions are matched to spheres by volume, since `tag()` numbers them in
traversal order; the radii are distinct to make that unambiguous.

Volume and centre of mass converge fast; the interfacial area is the weak one,
being a sum over reconstructed VOF planes, hence its looser tolerance.

The constants below mirror the .c sources and must be kept in step by hand.

Usage:
  python3 test_count_phase.py [--droplets <file>] [--bubbles <file>]

Exits non-zero (and prints FAIL) when any region deviates by more than its
tolerance, or when the region count is wrong.
"""
import math
import sys

# Mirrors the .c sources
R_SPH = [0.18, 0.14, 0.10]
C_SPH = [(-0.25, -0.25, -0.25),
         (0.25, 0.25, -0.25),
         (0.20, -0.28, 0.28)]
U0 = (0.1, 0.2, 0.3)
OMEGA = 1.5

VOLUME_TOL = 0.01   # relative
AREA_TOL = 0.03     # relative
POS_TOL = 0.005     # absolute, in units of L0
VEL_TOL = 0.01      # absolute


def read_regions(filename):
  """Read the last time step's rows, skipping comment (`#`) and blank lines."""
  rows = []
  with open(filename) as f:
    for line in f:
      line = line.strip()
      if not line or line.startswith('#'):
        continue
      values = [float(v) for v in line.split()]
      if len(values) not in (9, 11):
        raise ValueError(f'unexpected row width {len(values)} '
                         f'(expected 9 in 2D or 11 in 3D)')
      rows.append(values)
  if not rows:
    return rows
  last = rows[-1][0]
  return [r for r in rows if r[0] == last]


def check_regions(filename, label):
  """Returns 0 on success, 1 on failure."""
  rows = read_regions(filename)
  if len(rows) != len(R_SPH):
    print(f'FAIL: {len(rows)} regions in {filename}, expected {len(R_SPH)}')
    return 1

  dim = 3 if len(rows[0]) == 11 else 2

  # Match regions to spheres by volume: tag() numbers them in traversal order.
  rows.sort(key=lambda r: r[3])
  order = sorted(range(len(R_SPH)), key=lambda k: R_SPH[k])

  worst = {'volume': 0., 'area': 0., 'b': 0., 'u': 0.}
  for row, k in zip(rows, order):
    volume, area = row[3], row[-1]
    b = row[4:4 + dim]
    u = row[4 + dim:4 + 2 * dim]
    r, c = R_SPH[k], C_SPH[k][:dim]

    if dim == 3:
      volume_exact = 4. / 3. * math.pi * r ** 3
      area_exact = 4. * math.pi * r ** 2
    else:
      volume_exact = math.pi * r ** 2
      area_exact = 2. * math.pi * r
    u_exact = (U0[0] - OMEGA * c[1], U0[1] + OMEGA * c[0], U0[2])[:dim]

    worst['volume'] = max(worst['volume'],
                          abs(volume - volume_exact) / volume_exact)
    worst['area'] = max(worst['area'], abs(area - area_exact) / area_exact)
    worst['b'] = max(worst['b'], max(abs(p - q) for p, q in zip(b, c)))
    worst['u'] = max(worst['u'], max(abs(p - q) for p, q in zip(u, u_exact)))

  print(f'Processed {len(rows)} regions in {dim}D ({label}). '
        f'Max relative error: volume {worst["volume"]:.3g}, '
        f'area {worst["area"]:.3g}. '
        f'Max absolute error: b {worst["b"]:.3g}, u {worst["u"]:.3g}')

  failures = []
  if worst['volume'] > VOLUME_TOL:
    failures.append('volume')
  if worst['area'] > AREA_TOL:
    failures.append('area')
  if worst['b'] > POS_TOL:
    failures.append('centre of mass')
  if worst['u'] > VEL_TOL:
    failures.append('velocity')

  if failures:
    print(f'FAIL: count_phase mismatch on {", ".join(failures)}.')
    return 1
  print('PASS: count_phase matches the sphere closed form.')
  return 0


CHECKS = {
  '--droplets': lambda f: check_regions(f, 'above = 1'),
  '--bubbles': lambda f: check_regions(f, 'above = 0'),
}

USAGE = ('Usage: python3 test_count_phase.py [--droplets <file>] '
         '[--bubbles <file>]')


def main():
  argv = sys.argv[1:]
  jobs = []  # (check_fn, filename), in the order given on the command line
  i = 0
  while i < len(argv):
    flag = argv[i]
    if flag not in CHECKS or i + 1 >= len(argv):
      print(USAGE)
      return 1
    jobs.append((CHECKS[flag], argv[i + 1]))
    i += 2

  if not jobs:
    print(USAGE)
    return 1

  status = 0
  for check, filename in jobs:
    print(f'--- {filename}')
    status |= check(filename)
  return status


if __name__ == '__main__':
  sys.exit(main())
