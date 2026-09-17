#!/usr/bin/env python3
"""
Checker for `interfacial_power.h`, shared by
`test_interfacial_power_static.c` and `test_interfacial_power_routes.c`. Each
writes one row per resolution; the norms are computed there (with MPI
reductions), the pass/fail criteria live here.

--static (`N err_uniform err_radial err_field err_linear`): a drop of radius R
  with a constant potential sigma*kappa, for which Psi_sigma has a closed form.
  err_uniform is the null of a uniform velocity, err_radial the non-zero answer
  of a radial one; both are second-order convergent once the fraction comes
  from Vofi rather than fraction(). Order is log2(err(N)/err(2N)) between
  successive rows, checked against 2 within ORDER_TOL -- loose, because the
  interface samples the grid differently at each resolution.

  err_field (the output field integrates to the return value) and err_linear
  (Psi_sigma is linear in u) are identities, not approximations, and are held
  to ROUNDOFF_TOL. Their values are reported as a bound rather than printed:
  they sit at the roundoff of a cancelling sum, so the digits depend on the
  order of the MPI reduction and would make the reference rank-count specific.
  The convergent columns carry no such cancellation and are stable to the
  three digits shown -- verified identical at 1 and 4 ranks.

--energy (`N err_area resid`): the identity Psi_sigma = -d(sigma A)/dt on a
  drop relaxing from a strongly deformed shape to its equilibrium circle.
  err_area compares the final area against the closed-form equal-area circle
  and is the guard that equilibrium was reached -- without it the residual
  compares against an arbitrary endpoint, which is what makes the *undamped*
  oscillating drop useless for this. Held to AREA_TOL. resid is the relative
  gap between the integrated power and the surface energy released; it must
  keep converging, at ORDER_ENERGY rather than the 2 asserted elsewhere, for
  the reason given at that constant.

--amr (`mode i err_routes psi`): the same comparison on a tree, twice. mode 0
  keeps the interface inside a band refined to maxlevel and must agree to
  ROUNDOFF_TOL, as on a uniform grid; mode 1 puts a level jump through the drop
  and must *disagree* by at least DIFFER_MIN, since iforce.h's prolongation
  swap applies to route 1 and not to route 2. The second half is what stops the
  first from passing vacuously.

--routes (`i err_routes psi`): the two routes of the header integrate the same
  u.f_sigma, one reading back the applied acceleration and one recomputing the
  potential. They agree to roundoff on a uniform grid, where iforce.h's
  prolongation swap is a no-op. Held to ROUNDOFF_TOL, normalised by |psi|, and
  psi itself must not be zero or the agreement is vacuous.

--force (`N err_field err_psi psi`): interfacial_power_curvature() is built on
  surface_tension_force(), so u.F_sigma from the standalone vector must match
  the density it returns (err_field, normalised by max|d|) and the density
  must integrate to the same Psi_sigma (err_psi, normalised by |psi|) -- both
  exact identities, held to ROUNDOFF_TOL. psi itself must not be zero or the
  comparison is vacuous.

Usage:
  python3 test_interfacial_power.py [--static <file>] [--routes <file>]
                                    [--force <file>]

Exits non-zero (and prints FAIL) if any identity exceeds ROUNDOFF_TOL, if any
observed order falls below 2 - ORDER_TOL, or if the routes test never sees a
non-zero Psi_sigma.
"""
import math
import sys

ROUNDOFF_TOL = 1e-12
ORDER_TOL = 0.25
ORDER_EXPECTED = 2.0

# Psi_sigma of the radial case is O(1), so a null that has not converged shows
# up here well before the order check would catch it.
NULL_TOL = 1e-3

# The relaxed drop must sit on its equal-area circle this closely, or it has
# not reached equilibrium and the endpoint of the energy comparison is
# arbitrary. Measured ~3e-5.
AREA_TOL = 1e-3

# The energy residual converges, but not cleanly at second order: measured
# 0.128, 0.0331, 0.0130 at D/Delta = 12.8, 25.6, 51.2, i.e. orders 1.95 then
# 1.35. Three points cannot separate "order ~1.4" from "approaching a floor
# near 1%", so the criterion asserts only what is supported -- that the
# residual keeps falling at first order or better.
ORDER_ENERGY = 1.0

# Below this the routes agreement is not being measured, only |psi| ~ 0.
PSI_MIN = 1e-6

# A level jump through the interface must move the two routes apart by at
# least this much, or the comparison is not sensitive to the swap it is meant
# to detect. Measured margin is far larger; this only has to be well clear of
# ROUNDOFF_TOL.
DIFFER_MIN = 1e-6


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


def check_order(name, series, expected=None, tol=None):
  """Verify the convergence order of one column. 0 on success, 1 on failure."""
  expected = ORDER_EXPECTED if expected is None else expected
  tol = ORDER_TOL if tol is None else tol
  if any(e <= 0. for e in series[1:]):
    print(f'FAIL: vanishing {name} error, cannot measure an order')
    return 1
  orders = [math.log2(c / f) for c, f in zip(series, series[1:])]
  shown_err = ', '.join(f'{e:.3g}' for e in series)
  shown_ord = ', '.join(f'{o:.2f}' for o in orders)
  if min(orders) < expected - tol:
    print(f'FAIL: {name} error {shown_err}, order {shown_ord} '
          f'(expected {expected:g})')
    return 1
  print(f'PASS: {name} error {shown_err}, order {shown_ord}.')
  return 0


def check_static(filename):
  """Verify the closed forms and the two identities. 0 on success, 1 on failure."""
  rows = read_table(filename, 5)
  if len(rows) < 2:
    print(f'FAIL: need at least two resolutions in {filename}, '
          f'got {len(rows)}')
    return 1

  print(f'Processed {len(rows)} resolutions.')
  status = 0

  worst = max(row[1] for row in rows)
  shown = ', '.join(f'{row[1]:.3g}' for row in rows)
  if worst > NULL_TOL:
    print(f'FAIL: uniform velocity is not a null: {shown}')
    status = 1
  else:
    print(f'PASS: uniform velocity gives zero power: {shown}')

  status |= check_order('radial', [row[2] for row in rows])

  for name, col, claim in (
      ('field', 3, 'the output field integrates to the return value'),
      ('linear', 4, 'Psi_sigma is linear in u')):
    worst = max(row[col] for row in rows)
    if worst > ROUNDOFF_TOL:
      shown = ', '.join(f'{row[col]:.3g}' for row in rows)
      print(f'FAIL: {claim} -- violated: {shown}')
      status = 1
    else:
      print(f'PASS: {claim}, below {ROUNDOFF_TOL:g}.')
  return status


def check_routes(filename):
  """Verify the two routes agree. Returns 0 on success, 1 on failure."""
  rows = read_table(filename, 3)
  if not rows:
    print(f'FAIL: no data rows in {filename}')
    return 1

  print(f'Processed {len(rows)} timesteps.')
  status = 0

  biggest = max(abs(row[2]) for row in rows)
  if biggest < PSI_MIN:
    print(f'FAIL: Psi_sigma never exceeds {PSI_MIN:g} '
          f'(max {biggest:.3g}); the agreement is vacuous')
    return 1

  worst = max(row[1] for row in rows)
  if worst > ROUNDOFF_TOL:
    print(f'FAIL: the two routes disagree: worst {worst:.3g} relative')
    status = 1
  else:
    print(f'PASS: the two routes agree below {ROUNDOFF_TOL:g} relative, '
          f'over |Psi_sigma| up to {biggest:.3g}.')
  return status


def check_force(filename):
  """Verify surface_tension_force() reproduces interfacial_power_curvature()'s
  decomposition, pointwise and integrated. Returns 0 on success, 1 on failure."""
  rows = read_table(filename, 4)
  if not rows:
    print(f'FAIL: no data rows in {filename}')
    return 1

  status = 0
  biggest = max(abs(row[3]) for row in rows)
  if biggest < PSI_MIN:
    print(f'FAIL: Psi_sigma never exceeds {PSI_MIN:g} '
          f'(max {biggest:.3g}); the comparison is vacuous')
    return 1

  for name, col in (('field', 1), ('psi', 2)):
    worst = max(row[col] for row in rows)
    if worst > ROUNDOFF_TOL:
      shown = ', '.join(f'{row[col]:.3g}' for row in rows)
      print(f'FAIL: {name} identity violated: {shown}')
      status = 1
    else:
      print(f'PASS: {name} identity holds below {ROUNDOFF_TOL:g}.')
  return status


def check_amr(filename):
  """Verify the tree behaviour on both sides. Returns 0 on success, 1 on failure."""
  rows = read_table(filename, 4)
  if not rows:
    print(f'FAIL: no data rows in {filename}')
    return 1

  status = 0
  for mode, name, expect_agreement in ((0, 'contained', True),
                                       (1, 'straddling', False)):
    sel = [row for row in rows if row[0] == mode]
    if not sel:
      print(f'FAIL: no rows for the {name} case')
      status = 1
      continue

    biggest = max(abs(row[3]) for row in sel)
    if biggest < PSI_MIN:
      print(f'FAIL: {name}: Psi_sigma never exceeds {PSI_MIN:g} '
            f'(max {biggest:.3g}); the comparison is vacuous')
      status = 1
      continue

    worst = max(row[2] for row in sel)
    if expect_agreement:
      if worst > ROUNDOFF_TOL:
        print(f'FAIL: {name}: the routes should agree on one resolution, '
              f'worst {worst:.3g} relative')
        status = 1
      else:
        print(f'PASS: {name}: the routes agree below {ROUNDOFF_TOL:g} '
              f'relative over {len(sel)} timesteps.')
    else:
      if worst < DIFFER_MIN:
        print(f'FAIL: {name}: a level jump through the interface left the '
              f'routes agreeing to {worst:.3g}; the test cannot see the '
              f'prolongation swap it is meant to detect')
        status = 1
      else:
        print(f'PASS: {name}: a level jump separates the routes by at least '
              f'{DIFFER_MIN:g} relative (worst {orderof(worst)}).')
  return status


def orderof(value):
  """Order of magnitude, so the reference does not depend on the rank count."""
  return f'1e{math.floor(math.log10(value)):d}' if value > 0. else '0'


def check_energy(filename):
  """Verify the energy identity on a relaxing drop. 0 on success, 1 on failure."""
  rows = read_table(filename, 3)
  if len(rows) < 2:
    print(f'FAIL: need at least two resolutions in {filename}, '
          f'got {len(rows)}')
    return 1

  print(f'Processed {len(rows)} resolutions.')
  status = 0

  worst = max(row[1] for row in rows)
  if worst > AREA_TOL:
    shown = ', '.join(f'{row[1]:.3g}' for row in rows)
    print(f'FAIL: the drop has not reached its equilibrium circle: {shown}; '
          f'the energy comparison has no well-defined endpoint')
    status = 1
  else:
    print(f'PASS: relaxed onto the equal-area circle, within {AREA_TOL:g}.')

  series = [row[2] for row in rows]
  if any(e <= 0. for e in series[1:]):
    print('FAIL: vanishing energy residual, cannot measure an order')
    return status | 1

  orders = [math.log2(c / f) for c, f in zip(series, series[1:])]

  # Two significant digits, and the order as a bound rather than a value: the
  # residual is a small difference of two much larger integrals accumulated
  # over thousands of timesteps, so its third digit and the orders derived
  # from it move with the domain decomposition. Measured at 1 and 4 ranks:
  # 0.12785/0.12750, 0.033126/0.032775, 0.0129702/0.0129740.
  shown = ', '.join(f'{e:.2g}' for e in series)
  if min(orders) < ORDER_ENERGY - ORDER_TOL:
    shown_ord = ', '.join(f'{o:.2f}' for o in orders)
    print(f'FAIL: energy residual {shown}, order {shown_ord} '
          f'(expected at least {ORDER_ENERGY:g})')
    status = 1
  else:
    print(f'PASS: energy residual {shown}, converging at order '
          f'{ORDER_ENERGY:g} or better.')
  return status


CHECKS = {'--static': check_static, '--routes': check_routes,
          '--force': check_force,
          '--amr': check_amr, '--energy': check_energy}

USAGE = ('Usage: python3 test_interfacial_power.py [--static <file>] '
         '[--routes <file>] [--force <file>] [--amr <file>] [--energy <file>]')


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
