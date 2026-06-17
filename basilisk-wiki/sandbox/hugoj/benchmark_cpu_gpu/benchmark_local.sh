#!/usr/bin/env bash

# Monitor with htop and nvtop

set -x

TODO='all'
NTHREADS=8
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

echo $SCRIPT_DIR
NAME='results_bench_local'

ln -sf turbulence_benchmark.c turbulence.c

##########################
# Hardware
rm "${NAME}_hardware"
echo "CPU" >"${NAME}_hardware"
cat /proc/cpuinfo | grep -i 'model name' >>"${NAME}_hardware"
echo "GPU" >>"${NAME}_hardware"
__NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia glxinfo -B >>"${NAME}_hardware"

##########################
# Actual benchmark

CASE='turbulence'
echo "Case is: ${CASE}.c"

FILE=${SCRIPT_DIR}/"${NAME}_${CASE}"
mv $FILE "${FILE}_old"
touch $FILE

list_res=(128 256) # 128 256 512 1024

if [ $TODO = 'cpu' ] || [ $TODO = 'all' ]; then
  echo 'cpu' >>$FILE
  cd $BASILISK/examples
  echo 'compile once on CPU' # in single precision
  CFLAGS='-DTRACE=2 -DSINGLE_PRECISION -DBENCHMARK=1' make turbulence.tst
  for i in "${list_res[@]}"; do
    echo $i >>$FILE
    OMP_NUM_THREADS=$NTHREADS turbulence/turbulence $i >>"${FILE}" 2>&1
  done
fi

list_res=(128 256 512 1024 2048 4096) # 128 256 512 1024
list_backend=('gpu' 'cuda' 'hip')

if [ $TODO = 'allgpu' ] || [ $TODO = 'all' ]; then
  for backend in "${list_backend[@]}"; do
    echo $backend >>$FILE
    cd $BASILISK/examples
    __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia OMP_NUM_THREADS=1 CFLAGS='-DTRACE=2 -DBENCHMARK=1' make turbulence.gpu.tst # compile first
    for i in "${list_res[@]}"; do
      echo $i >>$FILE
      __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia OMP_NUM_THREADS=1 turbulence.gpu/turbulence.gpu $i >>"${FILE}" 2>&1
    done
  done
fi

# if [ $TODO = 'gpu' ] || [ $TODO = 'all' ] || [ $TODO = 'allgpu1' ]; then
#   echo '-> GPU'
#   FILE_PATH="${SCRIPT_DIR}/${name}_turbulence_gpu_"
#   cd $BASILISK/examples
#   __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia OMP_NUM_THREADS=1 CFLAGS='-DTRACE=2 -DBENCHMARK=1' make turbulence.gpu.tst # compile first
#   for i in "${list_res[@]}"; do
#     echo $i
#     __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia OMP_NUM_THREADS=1 turbulence.gpu/turbulence.gpu $i >"${FILE_PATH}${i}" 2>&1
#   done
# fi
# if [ $TODO = 'cuda' ] || [ $TODO = 'all' ] || [ $TODO = 'allgpu' ]; then
#   echo '-> CUDA'
#   FILE_PATH="${SCRIPT_DIR}/${name}_turbulence_cuda_"
#   cd $BASILISK/examples
#   __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia OMP_NUM_THREADS=1 CFLAGS='-DTRACE=2 -DBENCHMARK=1' make turbulence.cuda.tst # compile first
#   for i in "${list_res[@]}"; do
#     echo $i
#     __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia OMP_NUM_THREADS=1 turbulence.cuda/turbulence.cuda $i >"${FILE_PATH}${i}" 2>&1
#   done
# fi
# if [ $TODO = 'hip' ] || [ $TODO = 'all' ] || [ $TODO = 'allgpu' ]; then
#   echo '-> HIP'
#   FILE_PATH="${SCRIPT_DIR}/${name}_turbulence_hip_"
#   cd $BASILISK/examples
#   __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia OMP_NUM_THREADS=1 CFLAGS='-DTRACE=2 -DBENCHMARK=1' make turbulence.hip.tst # compile first
#   for i in "${list_res[@]}"; do
#     echo $i
#     __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia OMP_NUM_THREADS=1 turbulence.hip/turbulence.hip $i >"${FILE_PATH}${i}" 2>&1
#   done
# fi

echo 'CASE 2'
list_res=(128, 256, 512, 1024, 2048)

echo 'TODO AFTER CASE 1 HAS BEEN DOWN FOR ALL HARDWARE'
