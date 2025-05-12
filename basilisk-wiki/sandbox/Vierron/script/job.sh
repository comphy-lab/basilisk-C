
#SBATCH -A pt
#SBATCH -p epic,epic2
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexandre.vierron@univ-amu.fr
#SBATCH -J test
#SBATCH -o ./%x-%A_%a.out
#SBATCH -e ./%x-%A_%a.err

. /usr/share/Modules/init/bash
module purge
module load openmpi-gcc/4.0.4-pmix_v2

export LD_LIBRARY_PATH=/home/vierron/local/lib:$LD_LIBRARY_PATH
export PATH=$PATH:/home/vierron/local/bin

#mpicc -Wall -O2 -std=c99 -D_XOPEN_SOURCE=700 _rb.c -o RBnf -I/home/vierron/ -L/home/vierron/gl -lfb_osmesa -L/home/vierron/local/lib -lglutils -lGLU -lOSMesa -lm

srun --mpi=pmix_v2 ./RBnf
