.PHONY: all clean

CFLAGS += -I$(BASILISK) -I.. -DSINGLE_PRECISION -O2

all: # libbuda.a

# these are not Basilisk programs
EXCLUDE = cuda.c

libbuda.a: cuda.o
	ar cr libbuda.a cuda.o

cuda.o: ../gpu/backend.h ../externals.h cuda.c

include ../../Makefile.defs
