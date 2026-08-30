#!/usr/bin/env python3
"""
Run a group of profile tests and summarise the results.

  ./run_tests.py profiles   # the profiles/ tests
  ./run_tests.py all        # every group
  ./run_tests.py            # lists the groups

This is only a driver: each test is still run by `make <test>.tst`, which
remains the single source of truth for how a test is compiled, executed and
compared against its reference file. Run `load_gcc` first, as usual.

Each test in a group may be run against several grids (`-grid=...`, passed
through `EXTRA_CFLAGS`): a plain test name runs with the default grid, a
tuple `(test, grid)` runs the same source with `-grid=<grid>` and compares
against `<test>.<grid>.ref` instead of `<test>.ref` when that file exists --
`runtest` only knows about `<test>.ref`, so the variant's reference is
copied into place before each run. The grid name alone selects the
dimension too (`multigrid3D`/`octree` are 3D, `multigrid`/`quadtree` are
2D), so no separate dimension flag or `.3D.c` symlink is needed.
"""
import os
import shutil
import subprocess
import sys

# Each entry is either a bare test name (default grid) or a (test, grid)
# pair. `test_profiles.c` (2D, unity weight), `test_profiles_weighted.c`
# (2D, rho weight -- coincides with the unweighted closed form since rho is
# constant on every sampled plane, but still exercises the weighted code
# path) and `test_profiles3D.c` (3D) are separate sources; multigrid/
# multigrid3D have no AMR near the interface, so they are expected to
# resolve the dissipation profile less well than quadtree/octree -- that is
# captured in each grid's own frozen .ref, not papered over.
# `test_average_bias.c` and `test_profile_bias.c` branch on
# `#if dimension == 3` internally and check their own tolerance via the
# process exit code, so they need no .ref at all -- in 2D or in 3D.
GROUPS = {
  'profiles': ([('test_profiles', 'multigrid'),
                ('test_profiles', 'quadtree'),
                ('test_profiles3D', 'multigrid3D'),
                ('test_profiles3D', 'octree'),
                ('test_profiles_weighted', 'multigrid'),
                ('test_profiles_weighted', 'quadtree'),
                'test_average_bias',
                ('test_average_bias', 'octree'),
                'test_profile_bias',
                ('test_profile_bias', 'octree')], False),
}

NP = os.environ.get('NP', '4')  # MPI ranks, when a group needs them


def run_one(entry, mpi):
  """Run a single test (optionally under a given grid) and return True on pass."""
  test, grid = entry if isinstance(entry, tuple) else (entry, None)
  ref_variant = f'{test}.{grid}.ref' if grid else f'{test}.ref'

  # Remove the stamp so the test runs again instead of being reported up to
  # date, and the object so that a switch between serial/MPI or grid recompiles.
  subprocess.run(['rm', '-rf', test, f'{test}.tst', f'{test}.s'], check=False)

  # `runtest` always diffs against `<test>.ref`; put the variant's reference
  # there for the duration of this run, when there is one.
  if grid and os.path.exists(ref_variant):
    shutil.copyfile(ref_variant, f'{test}.ref')

  env = dict(os.environ)
  if mpi:
    env['CC'] = f'mpicc -D_MPI={NP}'
  if grid:
    env['EXTRA_CFLAGS'] = f'-grid={grid} ' + env.get('EXTRA_CFLAGS', '')

  subprocess.run(['make', f'{test}.tst'], env=env,
                 stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

  label = f'{test} [grid={grid}]' if grid else test
  if mpi:
    label += f' [mpi -np {NP}]'
  if os.path.exists(os.path.join(test, 'fail')):
    verdict = 'FAIL'
  elif os.path.exists(os.path.join(test, 'pass')):
    verdict = 'pass'
  else:
    verdict = 'ERROR (did not build or run)'
  print(f'  {label:<44} {verdict}')
  return verdict == 'pass'


def main():
  if len(sys.argv) != 2 or sys.argv[1] in ('-h', '--help'):
    print(__doc__.strip())
    print('\ngroups: ' + ' | '.join(list(GROUPS) + ['all']))
    return 2

  requested = sys.argv[1]
  if requested != 'all' and requested not in GROUPS:
    print(f'unknown group {requested!r}; '
          f'choose from {", ".join(list(GROUPS) + ["all"])}', file=sys.stderr)
    return 2

  # The environment comes from load_gcc; do not try to reproduce it here.
  if not os.environ.get('BASILISK'):
    print('BASILISK is not set -- run load_gcc first.', file=sys.stderr)
    return 2

  os.chdir(os.path.dirname(os.path.abspath(__file__)))

  names = list(GROUPS) if requested == 'all' else [requested]
  results = []
  for name in names:
    entries, mpi = GROUPS[name]
    print(f'== {name}')
    for entry in entries:
      results.append(run_one(entry, mpi))

  failed = results.count(False)
  print()
  if failed:
    print(f'{failed} of {len(results)} test(s) failed '
          f'(see <test>/fail for the diff against the reference)')
  else:
    print(f'all {len(results)} test(s) passed')
  return 1 if failed else 0


if __name__ == '__main__':
  sys.exit(main())
