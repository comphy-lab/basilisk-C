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

CFLAGS=-I$CUDADIR/include/ # cuda.h and nvrtc.h

# qcc, bview
ln -sf config.gcc config # fixme: use proper bigfoot config
make

# gpu (opengl)
# I need glfw.h that is exposed with the glfw package of the basilisk_profile
cd $BASILISK/grid/gpu
CFLAGS=$CFLAGS make

# gpu (cuda)
# I cuda.h and nvrtc.h that are exposed with the cuda_env
cd $BASILISK/grid/cuda
CFLAGS=$CFLAGS make

# gpu (HIP) --> to do later on dedicated nodes
#cd $BASILISK/grid/hip
#make
