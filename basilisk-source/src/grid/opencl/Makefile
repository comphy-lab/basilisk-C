.PHONY: all clean

CFLAGS += -I$(BASILISK) -I.. -DSINGLE_PRECISION -O2

all: # libocl.a

# these are not Basilisk programs
EXCLUDE = opencl.c

libocl.a: opencl.o
	ar cr libocl.a opencl.o

opencl.o: ../gpu/backend.h ../externals.h opencl.c

include ../../Makefile.defs
