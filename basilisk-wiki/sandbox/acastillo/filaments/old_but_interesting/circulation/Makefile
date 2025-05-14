include $(BASILISK)/Makefile.defs
export BASILISK_INCLUDE_PATH=/home/castillo/Workspace/wiki/sandbox/Antoonvh/

CC99 = gcc -std=c99
CFLAGS += -Wall -std=c99 -O3
OPENGLIBS=-lglutils -lfb_osmesa -lGLU -lOSMesa -lgsl -lgslcblas

base:
	export OMP_NUM_THREADS=1;
	export KMP_NUM_THREADS=1;
	qcc $(QCFLAGS) -source helix2_base.c ;
	mpicc $(QCFLAGS) $(CFLAGS) _helix2_base.c -I$(BASILISK) -L$(BASILISK)/gl $(OPENGLIBS) -o helix2_base/helix2_base -lm ;
	cd helix2_base ; mpirun -np 8 ./helix2_base > out 2> log ; cd ..
