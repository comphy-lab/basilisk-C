#!/usr/bin/env python3
"""
Checker for `profiles_slab.h`, shared by `test_profiles_slab.c` and
`test_profiles_slab_bias.c`. Each writes the same profile twice, once through
`profiles.h` and once through `profiles_slab.h`; this compares the two files
column by column.

--pair <a> <b>: the two must agree to roundoff. They accumulate in different
  orders -- one sums over sampled planes, the other over cells in grid order --
  so the comparison is against PAIR_TOL rather than bit-for-bit, even though on
  the present fixtures they mostly come out identical. Deviations are measured
  against the largest value in the file; see compare().

--differ <a> <b>: the same comparison where the header says the two are *not*
  interchangeable, because the sampled region spans more than one resolution.
  Passes when the deviation exceeds DIFFER_MIN, i.e. when the precondition is
  demonstrably load-bearing. If this ever fails, `profiles_slab.h` has become
  more general than its documentation claims and the header should say so.

Deviations are reported as an order of magnitude rather than a value: MPI
reductions sum in rank order, so the last digits move with the rank count and
a literal value would make the frozen reference depend on how many ranks the
test happened to run under.

Usage:
  python3 test_profiles_slab.py [--pair <a> <b>] [--differ <a> <b>] ...

Exits non-zero (and prints FAIL) if a --pair disagrees beyond PAIR_TOL or a
--differ agrees more closely than DIFFER_MIN.
"""
import math
import sys

PAIR_TOL = 1e-10
DIFFER_MIN = 1e-6


def read_profile(filename):
  """Read a profile file as a list of rows of floats, ignoring comments."""
  rows = []
  with open(filename) as f:
    for line in f:
      line = line.strip()
      if not line or line.startswith('#'):
        continue
      rows.append([float(v) for v in line.split()])
  return rows


def compare(file_a, file_b):
  """Max deviation between two profiles, relative to the file's overall scale.

  Returns (deviation, message); deviation is None if the files cannot be
  compared at all, in which case the message says why.
  """
  a, b = read_profile(file_a), read_profile(file_b)
  if not a or not b:
    return None, f'no data rows in {file_a if not a else file_b}'
  if len(a) != len(b):
    return None, f'row counts differ: {len(a)} vs {len(b)}'
  if len(a[0]) != len(b[0]):
    return None, f'column counts differ: {len(a[0])} vs {len(b[0])}'

  # Every deviation is measured against the largest value in the file, not
  # against the entry it came from. These profiles carry columns that are
  # analytically zero and hold nothing but roundoff; per-entry normalisation
  # makes such a column its own yardstick and turns 1e-17 against 1e-17 into
  # O(1). Against the file's own scale, a genuine error in even a small column
  # still lands far above PAIR_TOL.
  biggest = max(abs(v) for row in a for v in row)
  if biggest == 0.:
    biggest = 1.

  worst = 0.
  for row_a, row_b in zip(a, b):
    for va, vb in zip(row_a, row_b):
      dev = abs(va - vb)/biggest
      if dev > worst:
        worst = dev
  return worst, ''


def magnitude(x):
  """Order of magnitude of x, as a string, or '0' when x vanishes."""
  if x == 0.:
    return '0'
  return f'1e{round(math.log10(x))}'


def check_pair(file_a, file_b):
  """The two must agree. Returns 0 on success, 1 on failure."""
  dev, why = compare(file_a, file_b)
  if dev is None:
    print(f'FAIL: {file_a} vs {file_b}: {why}')
    return 1
  if dev > PAIR_TOL:
    print(f'FAIL: {file_a} vs {file_b} disagree, '
          f'max relative deviation ~{magnitude(dev)}')
    return 1
  print(f'PASS: {file_a} vs {file_b} agree to ~{magnitude(dev)}')
  return 0


def check_differ(file_a, file_b):
  """The two must disagree. Returns 0 on success, 1 on failure."""
  dev, why = compare(file_a, file_b)
  if dev is None:
    print(f'FAIL: {file_a} vs {file_b}: {why}')
    return 1
  if dev < DIFFER_MIN:
    print(f'FAIL: {file_a} vs {file_b} agree to ~{magnitude(dev)}; the '
          f'documented precondition looks unnecessary')
    return 1
  print(f'PASS: {file_a} vs {file_b} differ by ~{magnitude(dev)}, '
        f'as the precondition requires')
  return 0


CHECKS = {'--pair': check_pair, '--differ': check_differ}

USAGE = ('Usage: python3 test_profiles_slab.py [--pair <a> <b>] '
         '[--differ <a> <b>] ...')


def main():
  argv = sys.argv[1:]
  jobs = []  # (check_fn, file_a, file_b), in command-line order
  i = 0
  while i < len(argv):
    flag = argv[i]
    if flag not in CHECKS or i + 2 >= len(argv):
      print(USAGE)
      return 1
    jobs.append((CHECKS[flag], argv[i + 1], argv[i + 2]))
    i += 3

  if not jobs:
    print(USAGE)
    return 1

  status = 0
  for check, file_a, file_b in jobs:
    status |= check(file_a, file_b)
  return status


if __name__ == '__main__':
  sys.exit(main())
