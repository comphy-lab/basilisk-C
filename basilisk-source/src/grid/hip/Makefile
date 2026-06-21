.PHONY: all clean

CFLAGS += -I$(BASILISK) -I.. -DSINGLE_PRECISION -O2

all: libhipcuda.a libhipamd.a

libhipcuda.a: hipcuda.o
	ar cr libhipcuda.a hipcuda.o
hipcuda.o: ../gpu/backend.h hipcuda.c
hipcuda.o: CFLAGS += -D__HIP_PLATFORM_NVIDIA__
hipcuda.c: hip.c
	ln -s -f hip.c hipcuda.c

libhipamd.a: hipamd.o
	ar cr libhipamd.a hipamd.o
hipamd.o: ../gpu/backend.h hipamd.c
hipamd.o: CFLAGS += -D__HIP_PLATFORM_AMD__
hipamd.c: hip.c
	ln -s -f hip.c hipamd.c

include ../../Makefile.defs
