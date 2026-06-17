#!/bin/bash
# name of the job
#OAR -n compil
# ressources, max time
#OAR -l /nodes=1/gpu=1,walltime=00:10:00
# select GPU model
#OAR -p gpumodel='V100'
# faster queue time, short jobs
#OAR -t devel
# stdout
#OAR --stdout out_compil
# stderr
#OAR --stderr err_compil
# accouting, project
#OAR --project pr-data-ocean

# print all commands
set -x

BASILISK=/home/jacqhugo/basilisk/src
WORK=/bettik/PROJECTS/pr-data-ocean/jacqhugo/
cd $WORK/bench_cleanr

# loading guix env
source /applis/site/guix-start.sh # start guix session
refresh_guix basilisk_profile     # load profile

# source nvidia env
source /applis/environments/cuda_env.sh 12.6

# what is the gpu config
nvidia-smi >>bigfoot_nvidia_hardware.md

NAME=turbulence

# -> Libs folders
LIBGPU=$HOME/basilisk/src/grid/gpu              # liberrors.a
LIBBUDA=$HOME/basilisk/src/grid/cuda            # libbuda.a, needed when using -grid=cuda/multigrid
LIB_NVRTC=/applis/bigfoot/cuda/cuda-12.6/lib64/ # libnvrtc.so
# We also need libcuda.so, but this is loaded only when using a gpu ! so compiling on frontal node is not possible ... unless you are in interactive mode

# location of cuda.h and nvrtc.h : $CUDADIR/include

# -> Flags
CFLAGS='-autolink -disable-dimensions  -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 '
CFLAGS+='-DTRACE=2 -DBENCHMARK=1'

# -> Compile
INCLUDE="-I${LIBBUDA} -I${LIB_NVRTC} -I${LIBGPU}"
LIBS="-L${LIBBUDA} -L${LIB_NVRTC} -L${LIBGPU} -lm"
mkdir -p turbulence.cuda
ln -sf turbulence.c turbulence.cuda.c
$BASILISK/qcc $CFLAGS $INCLUDE -grid=cuda/multigrid -o $NAME.cuda/$NAME.cuda $NAME.cuda.c $LIBS
