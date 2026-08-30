#!/usr/bin/env python3
"""
Checker for the `profile_*_foreach_region()` ASCII outputs.

Verifies profile `.asc` files against the closed form of the two-phase
shear flow set up by `tests_profiles/test_profiles.c`:

  f = 1 for y < 0 (phase 1), f = 0 for y > 0 (phase 2)
  rho = rho1*f + rho2*(1-f), mu = mu1*f + mu2*(1-f)
  U_i(x,y) = alpha_shear*x + beta_i*y
  V_i(x,y) = -alpha_shear*y

with beta_i = tau/mu_i so the shear stress mu_i*beta_i is continuous across
the interface, even though beta_i (and dU/dy) jumps. Since the interface is
horizontal, every profile plane lies entirely within one phase, so each
height gets an exact, phase-uniform closed-form target.

Usage:
  python3 test_profiles.py --scalar profiles_scalar.asc
                            --product profiles_product.asc
                            --dissipation profiles_dissipation.asc
                            --scalar3D profiles_scalar3D.asc
                            --product3D profiles_product3D.asc
                            --dissipation3D profiles_dissipation3D.asc

Any subset of the flags may be given. Exits non-zero (and prints FAIL) when
any profile height's values deviate from the closed form by more than the
tolerance.
"""
import sys

RHO1, RHO2 = 3., 1.
MU1, MU2 = 3., 1.
TAU = 1.
ALPHA_SHEAR = 0.5
BETA1, BETA2 = TAU / MU1, TAU / MU2
TOL = 0.05


def read_profile(filename):
  """Read a `profile_*_foreach_region()` .asc file.

  Returns a list of float rows, skipping comment (`#`) and blank lines.
  """
  rows = []
  with open(filename) as f:
    for line in f:
      line = line.strip()
      if not line or line.startswith('#'):
        continue
      rows.append([float(v) for v in line.split()])
  return rows


def check_scalar(filename):
  """Verify `profile_foreach_region()` on the two-phase field.

  Row layout (2D): iprof, y, delta, mean(u.x), mean(u.x^2),
                    mean(u.y), mean(u.y^2), mean(rho), mean(rho^2)

  <U>(y)    = alpha_shear*<x> + beta_i*y = beta_i*y   (domain centered at
                                                        X0=-L0/2, so <x>=0)
  <U^2>(y)  = alpha_shear^2*<x^2> + beta_i^2*y^2  (cross term vanishes: <x>=0)
  <V>(y)    = -alpha_shear*y
  <rho>(y)  = rho1  (y<0)  or  rho2  (y>0)

  Returns 0 on success, 1 on failure.
  """
  rows = read_profile(filename)
  if not rows:
    print(f"FAIL: no data rows in {filename}")
    return 1

  L0 = 1.0  # domain size used by test_profiles.c
  max_err = 0.
  for row in rows:
    if len(row) != 9:
      print(f"FAIL: unexpected row width {len(row)} (expected 9)")
      return 1
    _, hprof, _, U_mean, U_sq, V_mean, V_sq, rho_mean, rho_sq = row

    phase1 = hprof < 0
    beta = BETA1 if phase1 else BETA2
    rho_expected = RHO1 if phase1 else RHO2
    U_sq_expected = ALPHA_SHEAR**2 * L0**2 / 12. + beta**2 * hprof**2

    errs = {
      '<U>': abs(U_mean - beta * hprof),
      '<V>': abs(V_mean - (-ALPHA_SHEAR * hprof)),
      '<rho>': abs(rho_mean - rho_expected),
      '<U^2>': abs(U_sq - U_sq_expected),
    }
    max_err = max(max_err, max(errs.values()))

  print(f"Processed {len(rows)} profile heights. Max error: {max_err:.3g}")
  if max_err > TOL:
    print("FAIL: profile_foreach_region mismatch on two-phase field.")
    return 1
  print("PASS: profile_foreach_region matches two-phase closed form.")
  return 0


def check_product(filename):
  """Verify `profile_product_foreach_region()` on the two-phase field.

  Row layout (2D): iprof, y, delta, mean(u.x*u.y), mean(u.y*u.y)

  Since <x> = 0 (domain centered at X0=-L0/2) and V does not depend on x,
  the cross terms vanish and the product profile reduces to a clean
  per-phase closed form:

    <U*V>(y) = alpha_shear*<x>*(-alpha_shear*y) + beta_i*y*(-alpha_shear*y)
             = -alpha_shear*beta_i*y^2
    <V*V>(y) = alpha_shear^2*y^2

  Returns 0 on success, 1 on failure.
  """
  rows = read_profile(filename)
  if not rows:
    print(f"FAIL: no data rows in {filename}")
    return 1

  max_err = 0.
  for row in rows:
    if len(row) != 5:
      print(f"FAIL: unexpected row width {len(row)} (expected 5)")
      return 1
    _, hprof, _, UV, VV = row

    phase1 = hprof < 0
    beta = BETA1 if phase1 else BETA2
    UV_expected = -ALPHA_SHEAR * beta * hprof**2
    VV_expected = ALPHA_SHEAR**2 * hprof**2

    errs = {'<U*V>': abs(UV - UV_expected), '<V*V>': abs(VV - VV_expected)}
    max_err = max(max_err, max(errs.values()))

  print(f"Processed {len(rows)} profile heights. Max error: {max_err:.3g}")
  if max_err > TOL:
    print("FAIL: profile_product_foreach_region mismatch on two-phase field.")
    return 1
  print("PASS: profile_product_foreach_region matches two-phase closed form.")
  return 0


def check_dissipation(filename):
  """Verify `profile_dissipation_foreach_region()` on the two-phase field.

  Row layout (2D): iprof, y, delta, grad[4], grad_sq[4], strain_sq[4]
    grad      = dudx, dvdx, dudy, dvdy
    grad_sq   = dudx^2, dvdx^2, dudy^2, dvdy^2
    strain_sq = nu*S^2, Sxx^2, Sxy^2, Syy^2

  Returns 0 on success, 1 on failure.
  """
  rows = read_profile(filename)
  if not rows:
    print(f"FAIL: no data rows in {filename}")
    return 1

  max_err = 0.
  for row in rows:
    if len(row) != 15:
      print(f"FAIL: unexpected row width {len(row)} (expected 15)")
      return 1
    _, hprof, _, dudx, dvdx, dudy, dvdy, \
        dudx_sq, dvdx_sq, dudy_sq, dvdy_sq, \
        nuS2, Sxx_sq, Sxy_sq, Syy_sq = row

    phase1 = hprof < 0
    beta = BETA1 if phase1 else BETA2
    nu = (MU1 if phase1 else MU2) / (RHO1 if phase1 else RHO2)
    S2 = 2. * ALPHA_SHEAR**2 + beta**2 / 2.

    errs = {
      'dudx': abs(dudx - ALPHA_SHEAR),
      'dvdx': abs(dvdx - 0.),
      'dudy': abs(dudy - beta),
      'dvdy': abs(dvdy - (-ALPHA_SHEAR)),
      'nu*S^2': abs(nuS2 - nu * S2),
      'Sxx^2': abs(Sxx_sq - ALPHA_SHEAR**2),
      'Sxy^2': abs(Sxy_sq - (beta / 2.)**2),
      'Syy^2': abs(Syy_sq - ALPHA_SHEAR**2),
    }
    max_err = max(max_err, max(errs.values()))

  print(f"Processed {len(rows)} profile heights. Max error: {max_err:.3g}")
  if max_err > TOL:
    print("FAIL: profile_dissipation_foreach_region mismatch on two-phase field.")
    return 1
  print("PASS: profile_dissipation_foreach_region matches two-phase closed form.")
  return 0


def check_scalar3D(filename):
  """Verify `profile_foreach_region()` on the 3D two-phase field.

  Row layout (3D): iprof, z, delta, mean(u.x), mean(u.x^2),
                    mean(u.y), mean(u.y^2), mean(u.z), mean(u.z^2),
                    mean(rho), mean(rho^2)

  Field (see test_profiles3D.c): U_i = alpha_shear*x + beta_i*z, V_i = 0,
  W_i = -alpha_shear*z. Domain centered at X0=Y0=-L0/2, so <x>=<y>=0:

    <U>(z)    = beta_i*z
    <U^2>(z)  = alpha_shear^2*<x^2> + beta_i^2*z^2   (cross term vanishes)
    <V>(z)    = 0
    <W>(z)    = -alpha_shear*z
    <rho>(z)  = rho1  (z<0)  or  rho2  (z>0)

  Returns 0 on success, 1 on failure.
  """
  rows = read_profile(filename)
  if not rows:
    print(f"FAIL: no data rows in {filename}")
    return 1

  L0 = 1.0  # domain size used by test_profiles3D.c
  max_err = 0.
  for row in rows:
    if len(row) != 11:
      print(f"FAIL: unexpected row width {len(row)} (expected 11)")
      return 1
    _, hprof, _, U_mean, U_sq, V_mean, V_sq, W_mean, W_sq, rho_mean, rho_sq = row

    phase1 = hprof < 0
    beta = BETA1 if phase1 else BETA2
    rho_expected = RHO1 if phase1 else RHO2
    U_sq_expected = ALPHA_SHEAR**2 * L0**2 / 12. + beta**2 * hprof**2

    errs = {
      '<U>': abs(U_mean - beta * hprof),
      '<V>': abs(V_mean - 0.),
      '<W>': abs(W_mean - (-ALPHA_SHEAR * hprof)),
      '<rho>': abs(rho_mean - rho_expected),
      '<U^2>': abs(U_sq - U_sq_expected),
    }
    max_err = max(max_err, max(errs.values()))

  print(f"Processed {len(rows)} profile heights. Max error: {max_err:.3g}")
  if max_err > TOL:
    print("FAIL: profile_foreach_region mismatch on 3D two-phase field.")
    return 1
  print("PASS: profile_foreach_region matches 3D two-phase closed form.")
  return 0


def check_product3D(filename):
  """Verify `profile_product_foreach_region()` on the 3D two-phase field.

  Row layout (3D): iprof, z, delta, mean(u.x*u.z), mean(u.z*u.z)

  Since <x> = 0 and W does not depend on x, the cross term vanishes:

    <U*W>(z) = -alpha_shear*beta_i*z^2
    <W*W>(z) = alpha_shear^2*z^2

  Returns 0 on success, 1 on failure.
  """
  rows = read_profile(filename)
  if not rows:
    print(f"FAIL: no data rows in {filename}")
    return 1

  max_err = 0.
  for row in rows:
    if len(row) != 5:
      print(f"FAIL: unexpected row width {len(row)} (expected 5)")
      return 1
    _, hprof, _, UW, WW = row

    phase1 = hprof < 0
    beta = BETA1 if phase1 else BETA2
    UW_expected = -ALPHA_SHEAR * beta * hprof**2
    WW_expected = ALPHA_SHEAR**2 * hprof**2

    errs = {'<U*W>': abs(UW - UW_expected), '<W*W>': abs(WW - WW_expected)}
    max_err = max(max_err, max(errs.values()))

  print(f"Processed {len(rows)} profile heights. Max error: {max_err:.3g}")
  if max_err > TOL:
    print("FAIL: profile_product_foreach_region mismatch on 3D two-phase field.")
    return 1
  print("PASS: profile_product_foreach_region matches 3D two-phase closed form.")
  return 0


def check_dissipation3D(filename):
  """Verify `profile_dissipation_foreach_region()` on the 3D two-phase field.

  Row layout (3D): iprof, z, delta, grad[9], grad_sq[9], strain_sq[7]
    grad      = dudx, dvdx, dudy, dvdy, dwdx, dwdy, dudz, dvdz, dwdz
    grad_sq   = same order, squared
    strain_sq = nu*S^2, Sxx^2, Sxy^2, Syy^2, Szz^2, Sxz^2, Syz^2

  Field: U_i = alpha_shear*x + beta_i*z, V_i = 0, W_i = -alpha_shear*z, so
  dudx=alpha_shear, dudz=beta_i, dwdz=-alpha_shear, everything else 0.
  Strain: Sxx=alpha_shear, Syy=0, Szz=-alpha_shear, Sxz=beta_i/2, Sxy=Syz=0.
  S^2 = Sxx^2 + Syy^2 + Szz^2 + 2*Sxy^2 + 2*Sxz^2 + 2*Syz^2
      = 2*alpha_shear^2 + beta_i^2/2

  Returns 0 on success, 1 on failure.
  """
  rows = read_profile(filename)
  if not rows:
    print(f"FAIL: no data rows in {filename}")
    return 1

  max_err = 0.
  for row in rows:
    if len(row) != 28:
      print(f"FAIL: unexpected row width {len(row)} (expected 28)")
      return 1
    _, hprof, _, dudx, dvdx, dudy, dvdy, dwdx, dwdy, dudz, dvdz, dwdz, \
        dudx_sq, dvdx_sq, dudy_sq, dvdy_sq, dwdx_sq, dwdy_sq, dudz_sq, dvdz_sq, dwdz_sq, \
        nuS2, Sxx_sq, Sxy_sq, Syy_sq, Szz_sq, Sxz_sq, Syz_sq = row

    phase1 = hprof < 0
    beta = BETA1 if phase1 else BETA2
    nu = (MU1 if phase1 else MU2) / (RHO1 if phase1 else RHO2)
    S2 = 2. * ALPHA_SHEAR**2 + beta**2 / 2.

    errs = {
      'dudx': abs(dudx - ALPHA_SHEAR),
      'dvdx': abs(dvdx - 0.),
      'dudy': abs(dudy - 0.),
      'dvdy': abs(dvdy - 0.),
      'dwdx': abs(dwdx - 0.),
      'dwdy': abs(dwdy - 0.),
      'dudz': abs(dudz - beta),
      'dvdz': abs(dvdz - 0.),
      'dwdz': abs(dwdz - (-ALPHA_SHEAR)),
      'nu*S^2': abs(nuS2 - nu * S2),
      'Sxx^2': abs(Sxx_sq - ALPHA_SHEAR**2),
      'Sxy^2': abs(Sxy_sq - 0.),
      'Syy^2': abs(Syy_sq - 0.),
      'Szz^2': abs(Szz_sq - ALPHA_SHEAR**2),
      'Sxz^2': abs(Sxz_sq - (beta / 2.)**2),
      'Syz^2': abs(Syz_sq - 0.),
    }
    max_err = max(max_err, max(errs.values()))

  print(f"Processed {len(rows)} profile heights. Max error: {max_err:.3g}")
  if max_err > TOL:
    print("FAIL: profile_dissipation_foreach_region mismatch on 3D two-phase field.")
    return 1
  print("PASS: profile_dissipation_foreach_region matches 3D two-phase closed form.")
  return 0


CHECKS = {
  '--scalar': check_scalar,
  '--product': check_product,
  '--dissipation': check_dissipation,
  '--scalar3D': check_scalar3D,
  '--product3D': check_product3D,
  '--dissipation3D': check_dissipation3D,
}


if __name__ == "__main__":
  argv = sys.argv[1:]
  jobs = []  # (check_fn, filename) pairs, in the order given on the command line
  i = 0
  while i < len(argv):
    flag = argv[i]
    if flag not in CHECKS or i + 1 >= len(argv):
      print("Usage: python3 test_profiles.py [--scalar <file>] "
            "[--product <file>] [--dissipation <file>] "
            "[--scalar3D <file>] [--product3D <file>] [--dissipation3D <file>]")
      sys.exit(1)
    jobs.append((CHECKS[flag], argv[i + 1]))
    i += 2

  if not jobs:
    print("Usage: python3 test_profiles.py [--scalar <file>] "
          "[--product <file>] [--dissipation <file>] "
          "[--scalar3D <file>] [--product3D <file>] [--dissipation3D <file>]")
    sys.exit(1)

  status = 0
  for check, filename in jobs:
    print(f"--- {filename}")
    status |= check(filename)
  sys.exit(status)
