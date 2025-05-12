minlevel=$1
maxlevel=$2
for np in 1 2 4 8 16 32; do
    mpirun -np $np ./mpi-laplacian $maxlevel $minlevel
done > res-$minlevel-$maxlevel
