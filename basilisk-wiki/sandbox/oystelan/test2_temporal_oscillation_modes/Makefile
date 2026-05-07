# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -fopenmp -Wdimensions -DMTRACE=3 -I$(DOCUMENT_ROOT)/sandbox

include $(BASILISK)/Makefile.defs

instability-new.c: instability.c
	ln -sf instability.c instability-new.c
instability-new.s: CFLAGS += -Dswitch_solver=1
instability-new.tst: CFLAGS += -Dswitch_solver=1

instability.tst: instability-new.tst
