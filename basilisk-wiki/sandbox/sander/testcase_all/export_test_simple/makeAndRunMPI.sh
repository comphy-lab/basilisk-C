
set -x
make clean && make export_test_simple_mpi && mpirun -np 8 --use-hwthread-cpus export_test_simple_mpi