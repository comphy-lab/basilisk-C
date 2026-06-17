#!/bin/bash
#SBATCH --job-name=compile      # nom du job
#SBATCH --ntasks=1
#SBATCH --time=00:05:00            # Temps d’exécution maximum demande (HH:MM:SS)
#SBATCH --output=out_compile%j.out  # Nom du fichier de sortie
#SBATCH --error=err_compile%j.out   # Nom du fichier d'erreur (ici commun avec la sortie)
#SBATCH --partition=compil
#SBATCH -A rkm@cpu

# echo des commandes
set -x

export HOME=/linkhome/rech/genige01/uji61qv

# nettoyage des modules charges en interactif et herites par defaut
module purge

module load arch/a100 # specific module for A100 GPUs

# chargement des modules
module load gcc/14.2.0
# NOTE: Basilisk needs to be compiled with the same version of GCC
module load cuda/12.6.3 # same version used to compile src/grid/cuda/libbuda.a
# NOTE: unlike other HPC, loading the module loads both libcuda.so and libnvrtc.so

LIBGPU=$HOME/basilisk/src/grid/gpu # needed when using -grid=cuda/multigrid

NAME=turbulence
CASE=cuda

# compilation
CFLAGS='-DTRACE=2 -DSINGLE_PRECISION -DBENCHMARK=1 '
CFLAGS+='-autolink -disable-dimensions -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2'
INCLUDE="-I${LIBGPU}"
LIBS="-lm -lnvrtc -lcuda"
mkdir -p ${NAME}.${CASE}
ln -sf ${NAME}.c ${NAME}.${CASE}.c
$BASILISK/qcc $CFLAGS $INCLUDE -grid=${CASE}/multigrid -o $NAME.${CASE}/$NAME.${CASE} $NAME.${CASE}.c $LIBS
