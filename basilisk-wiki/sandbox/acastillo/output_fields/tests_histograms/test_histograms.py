#!/usr/bin/env python3
"""
Checker for `probability_distribution_1D()` / `probability_distribution_1D_weighted()`'s
ASCII output.

Verifies the histogram/PDF/CDF of `c(x,y) = 0.5*(1 + tanh(y/eps))` written by
`tests_histograms/test_histograms.c` (unrestricted) and
`tests_histograms/test_histograms_weighted.c` (restricted to `f = (y < yf)`).
c is uniform in x and monotonically increasing in y, so on any y-range
[Y0, Y0+L0] where it is uniformly sampled/weighted, the CDF has a closed
form:

  CDF(c) = (y(c) - Y0)/L0,  y(c) = eps*atanh(2c-1)

For the unrestricted case L0 is the domain size; for the f-restricted case
it is instead the restricted region's extent (yf - Y0), since c(y) is
monotonic so restricting to y < yf is equivalent to renormalizing the same
CDF shape over that sub-range.

p_range[ii] is the *left* edge of bin ii, and p_sum[ii] is checked directly
against CDF(p_range[ii]) -- this is the primary correctness check, and is
tight (order db) since the CDF is exact at bin edges up to normalization
roundoff.

p_pdf[ii] is the histogram's bin-average density, not a pointwise sample,
so it is compared against the bin-averaged closed form
(CDF(right edge) - CDF(left edge)) / db rather than against the pointwise
derivative eps/(4*L0*c*(1-c)) -- the two only agree in the continuum limit.
Even bin-averaged, the PDF stays noticeably noisier than the CDF: each bin's
count is an integer number of grid rows, so its density is quantized in
steps of 1/(vol_cells*db) per row, and near c=0.5 (where dc/dy is largest)
a single row's width in c can exceed a bin's width, making that quantization
coarse. The PDF check therefore uses a looser tolerance than the CDF check,
and is informational rather than the primary pass/fail signal.

Usage:
  python3 test_histograms.py --histogram histogram.asc \
                              --eps 0.05 --Y0 -0.5 --L0 1.0

Exits non-zero (and prints FAIL) when the stored CDF deviates from the
closed form by more than CDF_TOL, away from the c=0,1 tails where atanh
diverges. The PDF error is reported but does not affect the exit code.
"""
import math
import sys

CDF_TOL = 0.01
PDF_TOL = 0.3


def read_histogram(filename):
  """Read a `get_distribution_function()` .asc file.

  Returns the last block as a list of (c, pdf, cdf) float rows, skipping
  comment (`#`) and blank lines.
  """
  blocks = [[]]
  with open(filename) as f:
    for line in f:
      line = line.strip()
      if line.startswith('#'):
        continue
      if not line:
        if blocks[-1]:
          blocks.append([])
        continue
      _, c, pdf, cdf = (float(x) for x in line.split())
      blocks[-1].append((c, pdf, cdf))
  return blocks[-1] if blocks[-1] else blocks[-2]


def check_histogram(filename, eps, Y0, L0):
  rows = read_histogram(filename)
  db = rows[1][0] - rows[0][0]

  def cdf(c):
    if c <= 0.:
      return 0.
    if c >= 1.:
      return 1.
    y = eps * math.atanh(2. * c - 1.)
    return min(max((y - Y0) / L0, 0.), 1.)

  cdf_err = 0.
  pdf_err = 0.
  for c, pdf, cdf_measured in rows:
    if not 0. < c < 1.:
      continue
    cdf_exact = cdf(c)
    pdf_exact = (cdf(c + db) - cdf(c)) / db
    # Skip the saturated tails: a handful of bins there are dominated by
    # discretization/roundoff, not a real error.
    if cdf_exact in (0., 1.) or pdf_exact == 0.:
      continue
    cdf_err = max(cdf_err, abs(cdf_measured - cdf_exact))
    pdf_err = max(pdf_err, abs(pdf - pdf_exact) / pdf_exact)
  status = 'PASS' if cdf_err < CDF_TOL else 'FAIL'
  pdf_status = 'ok' if pdf_err < PDF_TOL else 'noisy'
  print(f'histogram: cdf max err = {cdf_err:.3g} [{status}], '
        f'pdf max relerr = {pdf_err:.3g} [{pdf_status}]')
  return 0 if status == 'PASS' else 1


def main():
  args = sys.argv[1:]
  filename = None
  eps = Y0 = L0 = None
  while args:
    flag, value, *args = args
    if flag == '--histogram':
      filename = value
    elif flag == '--eps':
      eps = float(value)
    elif flag == '--Y0':
      Y0 = float(value)
    elif flag == '--L0':
      L0 = float(value)
    else:
      sys.exit(f'unknown flag {flag!r}')

  if filename is None:
    sys.exit('--histogram is required')

  sys.exit(check_histogram(filename, eps, Y0, L0))


if __name__ == '__main__':
  main()
