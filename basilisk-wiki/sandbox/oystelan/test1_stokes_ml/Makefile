# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -fopenmp -Wdimensions -DMTRACE=3 -I$(DOCUMENT_ROOT)/sandbox

include $(BASILISK)/Makefile.defs

stokes-new.c: stokes.c
	ln -sf stokes.c stokes-new.c
stokes-new.s: CFLAGS += -Dswitch_solver=1
stokes-new.tst: CFLAGS += -Dswitch_solver=1

stokes.tst: stokes-new.tst
