.PHONY: all clean

# Change -D__HIP_PLATFORM_NVIDIA__ to -D__HIP_PLATFORM_AMD__ to use ROCm on AMD

CFLAGS += -I$(BASILISK) -I.. -DSINGLE_PRECISION -O2 -D__HIP_PLATFORM_NVIDIA__

all: libhip.a

libhip.a: hip.o
	ar cr libhip.a hip.o

cuda.o: ../gpu/backend.h hip.c

include ../../Makefile.defs
