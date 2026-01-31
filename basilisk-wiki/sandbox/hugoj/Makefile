CFLAGS += -O2 -I$(DOCUMENT_ROOT)/sandbox -Wdimensions

all: ml_breaking_simple.tst plot

plot: ml_mpi/plots ml_openmp/plots ml_breaking_simple/plots

ml_mpi.c: ml_breaking_simple.c
	ln -sf ml_breaking_simple.c ml_mpi.c

ml_openmp.c: ml_breaking_simple.c
	ln -sf ml_breaking_simple.c ml_openmp.c

ml_mpi.tst: CC = mpicc -D_MPI=16

ml_openmp.tst: CFLAGS += -fopenmp

ml_breaking_simple.tst: ml_mpi.tst ml_openmp.tst

include $(BASILISK)/Makefile.defs
