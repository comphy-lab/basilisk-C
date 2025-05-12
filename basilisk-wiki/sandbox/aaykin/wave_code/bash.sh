#define number of threads
export OMP_NUM_THREADS=11
qcc -fopenmp -autolink -DgnuX=1 -O2 -g -Wall -Wno-unused-function -pipe -D_FORTIFY_SOURCE=2 wave.c -o wave -lm
./wave
#erase the executable
rm wave
