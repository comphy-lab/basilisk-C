.PHONY: all clean

CFLAGS += -I$(BASILISK) -I.. -DSINGLE_PRECISION -Wno-implicit-int -O2

all: # libhipcuda.a libhipamd.a

# these are not Basilisk programs
EXCLUDE = hip.c

libhipcuda.a: hipcuda.o
	ar cr libhipcuda.a hipcuda.o
hipcuda.o: ../gpu/backend.h ../externals.h hipcuda.c
hipcuda.o: CFLAGS += -D__HIP_PLATFORM_NVIDIA__
hipcuda.c: hip.c
	ln -s -f hip.c hipcuda.c

libhipamd.a: hipamd.o
	ar cr libhipamd.a hipamd.o
hipamd.o: ../gpu/backend.h ../externals.h hipamd.c
hipamd.o: CFLAGS += -D__HIP_PLATFORM_AMD__
hipamd.c: hip.c
	ln -s -f hip.c hipamd.c

include ../../Makefile.defs
