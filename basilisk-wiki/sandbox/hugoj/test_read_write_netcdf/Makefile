CFLAGS += -O2 -I$(DOCUMENT_ROOT)/sandbox 
CFLAGS += -disable-dimensions # skip dimensions check for now

all: myexample.tst

myexample.tst: CC = mpicc -D_MPI=16

_myexample.c: CFLAGS += -D_MPI=1
# Note: using -D_MPI=1 make the code crash, even without netcdf_bas.h !

include $(BASILISK)/Makefile.defs
