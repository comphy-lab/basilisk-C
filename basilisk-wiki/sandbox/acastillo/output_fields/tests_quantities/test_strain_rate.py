#!/usr/bin/env python3
"""
Checker for `strain_rate_sq()`, shared by `test_strain_rate_affine.c` and
`test_strain_rate_smooth.c`. Each writes one row per resolution; the norms are
computed there (with MPI reductions), the pass/fail criteria live here.

--affine (`N err_rel`): u_i = A_ij x_j with tr(A) = 0 has a constant gradient,
  so centred differences reproduce S^2 = S_ij S_ij exactly. The residual is
  floating-point cancellation in (u[1] - u[-1])/(2*Delta); an algebraic mistake
  instead gives an O(1) relative error at every resolution. Threshold
  AFFINE_TOL.

--smooth (`N err_max err_l2`): second-order differences on a smooth field, so
  both norms must fall as Delta^2. Order is log2(err(N)/err(2N)) between
  successive rows, checked against 2 within ORDER_TOL -- loose, since the
  location of the max norm moves between grids.

Usage:
  python3 test_strain_rate.py [--affine <file>] [--smooth <file>]

Exits non-zero (and prints FAIL) if the affine error exceeds AFFINE_TOL or any
observed order falls below 2 - ORDER_TOL.
"""
import math
import sys

AFFINE_TOL = 1e-10
ORDER_TOL = 0.25
ORDER_EXPECTED = 2.0


def read_table(filename, width):
  """Read the error table, skipping comment (`#`) and blank lines."""
  rows = []
  with open(filename) as f:
    for line in f:
      line = line.strip()
      if not line or line.startswith('#'):
        continue
      values = line.split()
      if len(values) != width:
        raise ValueError(f'unexpected row width {len(values)} '
                         f'(expected {width})')
      rows.append([int(values[0])] + [float(v) for v in values[1:]])
  return rows


def check_affine(filename):
  """Verify exactness on the affine field. Returns 0 on success, 1 on failure."""
  rows = read_table(filename, 2)
  if not rows:
    print(f'FAIL: no data rows in {filename}')
    return 1

  worst = max(row[1] for row in rows)
  shown = ', '.join(f'{row[1]:.3g}' for row in rows)
  print(f'Processed {len(rows)} resolutions. Relative error: {shown}')
  if worst > AFFINE_TOL:
    print('FAIL: strain_rate_sq is not exact on an affine field.')
    return 1
  print('PASS: strain_rate_sq exact to roundoff on an affine field.')
  return 0


def check_smooth(filename):
  """Verify second-order convergence. Returns 0 on success, 1 on failure."""
  rows = read_table(filename, 3)
  if len(rows) < 2:
    print(f'FAIL: need at least two resolutions in {filename}, '
          f'got {len(rows)}')
    return 1

  print(f'Processed {len(rows)} resolutions.')
  status = 0
  for name, col in (('max', 1), ('L2', 2)):
    series = [row[col] for row in rows]
    if any(e <= 0. for e in series[1:]):
      print(f'FAIL: vanishing {name}-norm error, cannot measure an order')
      status = 1
      continue
    orders = [math.log2(c / f) for c, f in zip(series, series[1:])]
    shown_err = ', '.join(f'{e:.3g}' for e in series)
    shown_ord = ', '.join(f'{o:.2f}' for o in orders)
    if min(orders) < ORDER_EXPECTED - ORDER_TOL:
      print(f'FAIL: {name}-norm error {shown_err}, order {shown_ord} '
            f'(expected {ORDER_EXPECTED:g})')
      status = 1
    else:
      print(f'PASS: {name}-norm error {shown_err}, order {shown_ord}.')
  return status


CHECKS = {'--affine': check_affine, '--smooth': check_smooth}

USAGE = ('Usage: python3 test_strain_rate.py [--affine <file>] '
         '[--smooth <file>]')


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
