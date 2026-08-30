#!/usr/bin/env python3
"""
Run a group of output tests and summarise the results.

  ./run_tests.py vtu        # the VTU tests
  ./run_tests.py all        # every group
  ./run_tests.py            # lists the groups

This is only a driver: each test is still run by `make <test>.tst`, which
remains the single source of truth for how a test is compiled, executed and
compared against its reference file. Run `load_gcc` first, as usual.
"""
import os
import subprocess
import sys

# Each group lists its tests and whether they have to run under MPI. The xmf
# and vtkhdf tests do: the HDF5 installation they need is an MPI build, so
# hdf5.h pulls in mpi.h and a serial compilation cannot succeed. That is a
# property of the test, not something worth asking for on the command line.
GROUPS = {
  'dump':     (['test_dump'], False),
  'vtu':      (['test_output_vtu',
                'test_output_vtu_box',
                'test_output_vtu_facets'], False),
  'xmf':      (['test_output_xmf',
                'test_output_xmf2',
                'test_output_xmf_box',
                'test_output_xmf_facets',
                'test_output_xmf_facets_list',
                'test_output_xmf_facets_box'], True),
  'vtkhdf':   (['test_output_vtkhdf',
                'test_output_vtkhdf2',
                'test_output_vtkhdf_box',
                'test_output_vtkhdf_facets',
                'test_output_vtkhdf_facets_box'], True),
  'catalyst': (['test_output_catalyst'], False),
}

NP = os.environ.get('NP', '4')  # MPI ranks, when a group needs them


def run_one(test, mpi):
  """Run a single test through make and return True when it passed."""
  # Remove the stamp so the test runs again instead of being reported up to
  # date, and the object so that a switch between serial and MPI recompiles.
  subprocess.run(['rm', '-rf', test, f'{test}.tst', f'{test}.s'], check=False)

  env = dict(os.environ)
  if mpi:
    env['CC'] = f'mpicc -D_MPI={NP}'

  subprocess.run(['make', f'{test}.tst'], env=env,
                 stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

  label = f'{test} [mpi -np {NP}]' if mpi else test
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
    tests, mpi = GROUPS[name]
    print(f'== {name}')
    for test in tests:
      results.append(run_one(test, mpi))

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
