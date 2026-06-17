#!/bin/bash

#OAR -n turbbench
#OAR -l /nodes=1/gpu=1,walltime=01:00:00
#OAR -p gpumodel='V100'
# If you are not in production, use OAR -t devel
#OAR --stdout out_turbbench
#OAR --stderr err_turbbench
#OAR --project pr-data-ocean

# print all commands
set -x

BASILISK=/home/jacqhugo/basilisk/src
WORK=/bettik/PROJECTS/pr-data-ocean/jacqhugo/
cd $WORK/benchmark_bas

# loading guix env
source /applis/site/guix-start.sh # start guix session
refresh_guix basilisk_profile     # load profile

# source nvidia env
source /applis/environments/cuda_env.sh 12.6

# what is the gpu config
nvidia-smi >bigfoot_nvidia_hardware.md

SCRIPT_DIR=/home/jacqhugo/work/benchmark_bas
FILE_PATH="${SCRIPT_DIR}/result_bench_bigfoot_case1_cuda_"

list_res=(128 256 512 1024 2048)
for i in "${list_res[@]}"; do
  echo $i
  #strace -e trace=open,openat,read,write,exit_group ./turbulence.cuda/turbulence.cuda 2>&1 | cat # debug
  turbulence.cuda/turbulence.cuda $i >"${FILE_PATH}${i}" 2>&1
done
