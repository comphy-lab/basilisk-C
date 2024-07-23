#!/bin/bash
#SBATCH -J basilisk
##SBATCH -d singleton
#SBATCH --nodes=1
##SBATCH --ntasks=2
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --time=10:00
#SBATCH --output basilisk.output
#SBATCH --exclusive

LEVEL=12

module purge

# module load gcc/4.9.2 
# module load intel/16.0.1
# module load bullxmpi/MPI3.gcc.4.9-beta
# mpicc -Wall -std=c99 -O2 -D_MPI=1 -DTRACE=2 _atomisation.c -o atomisation -lm

module load bullxmpi
module load intel

srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS ./atomisation $LEVEL \
     2> log-$LEVEL-$SLURM_NTASKS > out-$LEVEL-$SLURM_NTASKS
