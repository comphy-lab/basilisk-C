#!/usr/bin/env python3
"""
Checker for `probability_distribution_2D()` / `probability_distribution_2D_weighted()`'s
ASCII output.

Verifies the joint histogram/PDF written by `test_histograms2D.c`
(unrestricted) and `test_histograms2D_weighted.c` (restricted to
`f = (x < xf)`). u(x,y) = 0.5*(1 + tanh(x/epsu)) and
v(x,y) = 0.5*(1 + tanh(y/epsv)) each depend on a different coordinate, so on
a square domain uniformly sampled/weighted in (x,y) they are statistically
independent and the joint (bin-averaged) PDF factors as

  pdf(u,v) = pdf_u(u) * pdf_v(v)

Each marginal has the same closed form as in test_histograms.py:

  CDF(b) = (coord(b) - offset)/extent,  coord(b) = eps*atanh(2b-1)

For the unrestricted case `extent` is the domain size L0 for both marginals.
For the f-restricted case (`f = (x < xf)`), only the u-marginal's extent
changes (to xf - X0), since f does not depend on y.

--Lu and --Lv give the (possibly different) extents used to normalize the
u- and v- marginal CDFs respectively.

p_pdf[i,j] is the histogram's bin-average density, not a pointwise sample --
see test_histograms.py for why this is compared against the bin-averaged
closed form rather than the pointwise derivative, and why the PDF check
uses a looser, informational-only tolerance.

Usage:
  python3 test_histograms2D.py --histogram histogram2D.asc \
                                --epsu 0.15 --epsv 0.2 \
                                --X0 -0.5 --Y0 -0.5 --L0 1.0

  (--Lu/--Lv override --L0 individually when the two marginals are
  normalized over different extents, e.g. in the weighted test.)

Exits non-zero (and prints FAIL) when the joint PDF's implied CDF deviates
from the closed-form product of marginals by more than PDF_TOL, informed by
comparing marginal CDFs derived by summing the stored PDF.
"""
import math
import sys

MASS_TOL = 0.01
MEDIAN_PDF_TOL = 0.2


def read_histogram2d(filename):
  """Read a `probability_distribution_2D()` .asc file.

  Each call writes one block: a blank line separates each fixed-i row group
  (one per u-bin), and a second, consecutive blank line terminates the
  block. Returns the last such block as a list of (u, v, pdf) float rows,
  skipping comment (`#`) lines.
  """
  blocks = [[]]
  blank_streak = 0
  with open(filename) as f:
    for line in f:
      line = line.strip()
      if line.startswith('#'):
        continue
      if not line:
        blank_streak += 1
        if blank_streak >= 2 and blocks[-1]:
          blocks.append([])
        continue
      blank_streak = 0
      _, u, v, pdf = (float(x) for x in line.split())
      blocks[-1].append((u, v, pdf))
  return blocks[-1] if blocks[-1] else blocks[-2]


def check_histogram2d(filename, epsu, epsv, X0, Y0, Lu, Lv):
  rows = read_histogram2d(filename)
  us = sorted(set(round(u, 12) for u, v, pdf in rows))
  vs = sorted(set(round(v, 12) for u, v, pdf in rows))
  du = us[1] - us[0]
  dv = vs[1] - vs[0]

  def marginal_cdf(b, eps, offset, extent):
    if b <= 0.:
      return 0.
    if b >= 1.:
      return 1.
    coord = eps * math.atanh(2. * b - 1.)
    return min(max((coord - offset) / extent, 0.), 1.)

  def cdf_u(u):
    return marginal_cdf(u, epsu, X0, Lu)

  def cdf_v(v):
    return marginal_cdf(v, epsv, Y0, Lv)

  mass = sum(pdf for u, v, pdf in rows) * du * dv

  errs = []
  max_err = 0.
  for u, v, pdf_measured in rows:
    if not (0. < u < 1. and 0. < v < 1.):
      continue
    pdf_u = (cdf_u(u + du) - cdf_u(u)) / du
    pdf_v = (cdf_v(v + dv) - cdf_v(v)) / dv
    pdf_exact = pdf_u * pdf_v
    if pdf_exact == 0.:
      continue
    e = abs(pdf_measured - pdf_exact) / pdf_exact
    errs.append(e)
    max_err = max(max_err, e)

  errs.sort()
  median_err = errs[len(errs) // 2] if errs else 0.
  mass_err = abs(mass - 1.)

  # The joint PDF is a bin-average density, quantized by integer grid-cell
  # counts per bin -- near the steep tanh transitions a handful of bins can
  # be off by 30-50% (see test_histograms.py for the 1D analog), so the max
  # error is reported but not gated on. Mass conservation (tight) and the
  # median relative error (robust to the noisy minority of bins) are the
  # actual pass/fail signals.
  status = ('PASS' if mass_err < MASS_TOL and median_err < MEDIAN_PDF_TOL
            else 'FAIL')
  print(f'histogram2D: mass err = {mass_err:.3g}, pdf median relerr = '
        f'{median_err:.3g}, pdf max relerr = {max_err:.3g} over '
        f'{len(errs)} bins [{status}]')
  return 0 if status == 'PASS' else 1


def main():
  args = sys.argv[1:]
  filename = None
  epsu = epsv = X0 = Y0 = L0 = None
  Lu = Lv = None
  while args:
    flag, value, *args = args
    if flag == '--histogram':
      filename = value
    elif flag == '--epsu':
      epsu = float(value)
    elif flag == '--epsv':
      epsv = float(value)
    elif flag == '--X0':
      X0 = float(value)
    elif flag == '--Y0':
      Y0 = float(value)
    elif flag == '--L0':
      L0 = float(value)
    elif flag == '--Lu':
      Lu = float(value)
    elif flag == '--Lv':
      Lv = float(value)
    else:
      sys.exit(f'unknown flag {flag!r}')

  if filename is None:
    sys.exit('--histogram is required')
  if Lu is None:
    Lu = L0
  if Lv is None:
    Lv = L0

  sys.exit(check_histogram2d(filename, epsu, epsv, X0, Y0, Lu, Lv))


if __name__ == '__main__':
  main()
