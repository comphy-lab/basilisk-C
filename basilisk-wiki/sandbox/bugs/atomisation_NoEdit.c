I am trying to simulate the atomisation case in the folder /example. I do not edit any code in the atomisation.c, but restarting the simulation is not possile. The below is the procedure I compile with MPI.
  
1. qcc -source -D_MPI=1 -grid=quadtree atomisation.c

2. mpicc -Wall -O2 -std=c99 _atomisation.c -o run -L$BASILISK/gl -lglutils -lfb_osmesa -lOSMesa -lGLU -lm

and the errors are printed on the screen:
[xiaopang-Inspiron-7466:128945] *** Process received signal ***
[xiaopang-Inspiron-7466:128945] Signal: Segmentation fault (11)
[xiaopang-Inspiron-7466:128945] Signal code: Address not mapped (1)
[xiaopang-Inspiron-7466:128945] Failing at address: 0x55c7ed1871c0
[xiaopang-Inspiron-7466:128944] *** Process received signal ***
[xiaopang-Inspiron-7466:128944] Signal: Segmentation fault (11)
[xiaopang-Inspiron-7466:128944] Signal code: Address not mapped (1)
[xiaopang-Inspiron-7466:128944] Failing at address: 0x56a5edab6b40
[xiaopang-Inspiron-7466:128945] [ 0] /lib/x86_64-linux-gnu/libpthread.so.0(+0x143c0)[0x7f5ccf0083c0]
[xiaopang-Inspiron-7466:128945] [ 1] ./run(+0x3d95b)[0x556beb6f695b]
[xiaopang-Inspiron-7466:128945] [ 2] ./run(+0x3f84c)[0x556beb6f884c]
[xiaopang-Inspiron-7466:128944] [ 0] [xiaopang-Inspiron-7466:128945] [ 3] ./run(+0x49e4a)[0x556beb702e4a]
[xiaopang-Inspiron-7466:128945] [ 4] ./run(+0x6088a)[0x556beb71988a]
[xiaopang-Inspiron-7466:128945] [ 5] ./run(+0x1aefa)[0x556beb6d3efa]
[xiaopang-Inspiron-7466:128945] [ 6] ./run(+0x4a3ad)[0x556beb7033ad]
[xiaopang-Inspiron-7466:128945] [ 7] /lib/x86_64-linux-gnu/libpthread.so.0(+0x143c0)[0x7f9798d953c0]
[xiaopang-Inspiron-7466:128944] [ 1] ./run(+0x3d95b)[0x5649ebd4495b]
[xiaopang-Inspiron-7466:128944] [ 2] ./run(+0x3f84c)[0x5649ebd4684c]
./run(+0x7198)[0x556beb6c0198]
[xiaopang-Inspiron-7466:128945] [ 8] [xiaopang-Inspiron-7466:128944] [ 3] ./run(+0x49e4a)[0x5649ebd50e4a]
[xiaopang-Inspiron-7466:128944] [ 4] ./run(+0x6088a)[0x5649ebd6788a]
[xiaopang-Inspiron-7466:128944] [ 5] ./run(+0x1aefa)[0x5649ebd21efa]
[xiaopang-Inspiron-7466:128944] [ 6] ./run(+0x4a3ad)[0x5649ebd513ad]
[xiaopang-Inspiron-7466:128944] [ 7] ./run(+0x7198)[0x5649ebd0e198]
[xiaopang-Inspiron-7466:128944] [ 8] /lib/x86_64-linux-gnu/libc.so.6(__libc_start_main+0xf3)[0x7f5ccee260b3]
[xiaopang-Inspiron-7466:128945] [ 9] ./run(+0x720e)[0x556beb6c020e]
[xiaopang-Inspiron-7466:128945] *** End of error message ***
/lib/x86_64-linux-gnu/libc.so.6(__libc_start_main+0xf3)[0x7f9798bb30b3]
[xiaopang-Inspiron-7466:128944] [ 9] ./run(+0x720e)[0x5649ebd0e20e]
[xiaopang-Inspiron-7466:128944] *** End of error message ***
--------------------------------------------------------------------------
Primary job  terminated normally, but 1 process returned
a non-zero exit code. Per user-direction, the job has been aborted.
--------------------------------------------------------------------------
--------------------------------------------------------------------------
mpirun noticed that process rank 0 with PID 0 on node xiaopang-Inspiron-7466 exited on signal 11 (Segmentation fault).

  
  Note that when compiling without mpicc, i.e., running the simulation with one core, it is OK.