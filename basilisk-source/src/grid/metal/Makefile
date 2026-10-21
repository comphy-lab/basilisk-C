.PHONY: all clean

# metal.c is Objective-C (the Metal API is not available in plain C)
CFLAGS += -I$(BASILISK) -I.. -DSINGLE_PRECISION -O2 -x objective-c -fobjc-exceptions

all: # libmetal.a

# these are not Basilisk programs
EXCLUDE = metal.c

libmetal.a: metal.o
	ar cr libmetal.a metal.o

metal.o: ../gpu/backend.h ../externals.h metal.c

include ../../Makefile.defs
