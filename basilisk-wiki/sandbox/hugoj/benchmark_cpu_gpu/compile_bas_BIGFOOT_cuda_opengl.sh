#!/bin/bash

#OAR -n compile
#OAR -l /nodes=1/gpu=1,walltime=00:05:00
#OAR -p gpumodel='V100'
#OAR --stdout out_compile
#OAR --stderr err_compile
#OAR --project pr-data-ocean
#OAR -t devel

# print all commands
set -x


export BASILISK=/home/jacqhugo/basilisk/src
cd $BASILISK

# start guix
source /applis/site/guix-start.sh
# source my guix profile
refresh_guix basilisk_profile

# source nvidia env
source /applis/environments/cuda_env.sh 12.6

# what is the gpu config
nvidia-smi

CFLAGS=-I$CUDADIR/include/

# qcc, bview
ln -sf config.gcc config # fixme: use proper bigfoot config
make qcc

# gpu (opengl)
cd $BASILISK/grid/gpu
CFLAGS=$CFLAGS make libgpu.a liberrors.a

# gpu (cuda)
cd $BASILISK/grid/cuda
CFLAGS=$CFLAGS make libbuda.a

# gpu (HIP) --> to do later on dedicated nodes
#cd $BASILISK/grid/hip
#make
