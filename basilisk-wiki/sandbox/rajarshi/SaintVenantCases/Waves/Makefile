# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3

include $(BASILISK)/Makefile.defs

implicit_weno.c: implicit.c
	ln -s implicit.c implicit_weno.c
implicit_weno.tst: CFLAGS += -DWENO=1

implicit.tst: implicit_weno.tst

