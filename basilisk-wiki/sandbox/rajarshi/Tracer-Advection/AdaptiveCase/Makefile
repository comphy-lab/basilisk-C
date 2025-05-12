# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3

include $(BASILISK)/Makefile.defs

advection-weno.c: advection.c
	ln -s advection.c advection-weno.c
advection-weno.tst: CFLAGS += -DWENO=1

advection.tst: advection-weno.tst
