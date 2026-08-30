#!/usr/bin/env python3
"""
Run the spectra tests and summarise the results.

  ./run_tests.py spectra   # every spectra.h test
  ./run_tests.py all
  ./run_tests.py           # lists the groups

Only a driver: each test is still run by `make <test>.tst`, which remains the
single source of truth for how it is compiled, executed and compared against
its reference. Run `load_gcc` first.

A `(test, grid)` pair passes `-grid=<grid>` via EXTRA_CFLAGS and compares
against `<test>.<grid>.ref`, copied into place first since `runtest` only
knows `<test>.ref`. Everything here is 3D -- `spectra_sample.h` refuses to
build otherwise -- so the grids are `octree` and `multigrid3D`.

`test_spectra_hdf5` needs hdf5.h, which `spectra_output.h` detects at compile
time; without it the test builds into a different program with a log no single
reference can cover, so it is skipped. The probe uses `mpicc`, since a
parallel HDF5 includes mpi.h from hdf5.h.
"""
import os
import shutil
import subprocess
import sys

GROUPS = {
  # test_spectra_amr refines, so it is registered on the tree grid only: on
  # multigrid the `#if TREE` guard leaves it with nothing to compare.
  'spectra': ([('test_spectra_sample', 'octree'),
               ('test_spectra_sample', 'multigrid3D'),
               ('test_spectra_modes', 'octree'),
               ('test_spectra_modes', 'multigrid3D'),
               ('test_spectra_amr', 'octree'),
               ('test_spectra_ascii', 'octree'),
               ('test_spectra_ascii', 'multigrid3D'),
               ('test_spectra_hdf5', 'octree'),
               ('test_spectra_hdf5', 'multigrid3D')], True),
}

NP = os.environ.get('NP', '4')  # MPI ranks, when a group needs them

# `multigrid3D` decomposes the domain into octants, so it insists on 8^i ranks
# and refuses the default 4.
GRID_NP = {'multigrid3D': '8'}

NEEDS_HDF5 = ('test_spectra_hdf5',)


def have_hdf5(cc='mpicc'):
  """Whether `#include <hdf5.h>` resolves for the compiler the tests use."""
  try:
    return subprocess.run([cc, '-E', '-x', 'c', '-'],
                          input='#include <hdf5.h>\n', text=True,
                          stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL).returncode == 0
  except OSError:
    return False


def run_one(entry, mpi):
  """Run a single test under a given grid and return True on pass."""
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
  hdf5 = have_hdf5()

  names = list(GROUPS) if requested == 'all' else [requested]
  results, skipped = [], 0
  for name in names:
    entries, mpi = GROUPS[name]
    print(f'== {name}')
    for entry in entries:
      test = entry[0] if isinstance(entry, tuple) else entry
      if test in NEEDS_HDF5 and not hdf5:
        print(f'  {test:<48} skipped (no hdf5.h)')
        skipped += 1
        continue
      results.append(run_one(entry, mpi))

  failed = results.count(False)
  print()
  tail = f' ({skipped} skipped)' if skipped else ''
  if failed:
    print(f'{failed} of {len(results)} test(s) failed{tail} '
          f'(see <test>/fail for the diff against the reference)')
  else:
    print(f'all {len(results)} test(s) passed{tail}')
  return 1 if failed else 0


if __name__ == '__main__':
  sys.exit(main())
