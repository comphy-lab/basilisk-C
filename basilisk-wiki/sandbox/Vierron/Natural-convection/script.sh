~~~~~~~~~~~~~~~~
#!/bin/sh
#SBATCH --exclusive
#SBATCH --ntasks=16   
#SBATCH --nodes=1 
#SBATCH --partition data                 
echo "Running on: $SLURM_NODELIST"
cd $PWD
export OMP_NUM_THREADS=1
module purge
module load intel/compiler
module load intel/mpi
module load slurm
unset I_MPI_PMI_LIBRARY
export I_MPI_EXTRA_FILESYSTEM=on
export I_MPI_EXTRA_FILESYSTEM_LIST=gpfs

mpirun ./cav
~~~~~~~~~~~~~~~~