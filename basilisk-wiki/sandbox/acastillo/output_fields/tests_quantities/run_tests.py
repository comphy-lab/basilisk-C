#!/usr/bin/env python3
"""
Run a group of tests and summarise the results.

  ./run_tests.py strain_rate   # the strain_rate.h tests
  ./run_tests.py all           # every group
  ./run_tests.py               # lists the groups

Only a driver: each test is still run by `make <test>.tst`, which remains the
single source of truth for how it is compiled, executed and compared against
its reference. Run `load_gcc` first.

A bare test name uses the default grid; a `(test, grid)` pair passes
`-grid=<grid>` via EXTRA_CFLAGS and compares against `<test>.<grid>.ref`, which
is copied into place first since `runtest` only knows `<test>.ref`. The grid
name selects the dimension too (`multigrid3D`/`octree` are 3D).
"""
import os
import shutil
import subprocess
import sys

GROUPS = {
  'strain_rate': ([('test_strain_rate_affine', 'multigrid'),
                   ('test_strain_rate_affine', 'multigrid3D'),
                   ('test_strain_rate_smooth', 'multigrid'),
                   ('test_strain_rate_smooth', 'multigrid3D')], True),
  # count_phase.h resolves the interface with AMR, so only the tree grids are
  # registered: on multigrid the `#if TREE` refinement is skipped and the
  # coarse base level would freeze a failure with nothing to learn from it.
  'count_phase': ([('test_count_phase', 'quadtree'),
                   ('test_count_phase', 'octree'),
                   ('test_count_phase_bubbles', 'quadtree'),
                   ('test_count_phase_bubbles', 'octree')], True),
}

NP = os.environ.get('NP', '4')  # MPI ranks, when a group needs them

# `multigrid3D` decomposes the domain into octants, so it insists on 8^i ranks
# and refuses the default 4.
GRID_NP = {'multigrid3D': '8'}


def run_one(entry, mpi):
  """Run a single test (optionally under a given grid) and return True on pass."""
  test, grid = entry if isinstance(entry, tuple) else (entry, None)
  ref_variant = f'{test}.{grid}.ref' if grid else f'{test}.ref'

  # Remove the stamp so the test runs again instead of being reported up to
  # date, and the object so a switch between serial/MPI or grid recompiles.
  subprocess.run(['rm', '-rf', test, f'{test}.tst', f'{test}.s'], check=False)

  if grid and os.path.exists(ref_variant):
    shutil.copyfile(ref_variant, f'{test}.ref')

  np = GRID_NP.get(grid, NP)
  env = dict(os.environ)
  if mpi:
    env['CC'] = f'mpicc -D_MPI={np}'
  if grid:
    env['EXTRA_CFLAGS'] = f'-grid={grid} ' + env.get('EXTRA_CFLAGS', '')

  subprocess.run(['make', f'{test}.tst'], env=env,
                 stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

  label = f'{test} [grid={grid}]' if grid else test
  if mpi:
    label += f' [mpi -np {np}]'
  if os.path.exists(os.path.join(test, 'fail')):
    verdict = 'FAIL'
  elif os.path.exists(os.path.join(test, 'pass')):
    verdict = 'pass'
  else:
    verdict = 'ERROR (did not build or run)'
  print(f'  {label:<48} {verdict}')
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
