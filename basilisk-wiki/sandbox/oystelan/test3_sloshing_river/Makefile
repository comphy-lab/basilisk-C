# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -fopenmp -Wdimensions -DMTRACE=3 -I$(DOCUMENT_ROOT)/sandbox

include $(BASILISK)/Makefile.defs

river-new.c: river.c
	ln -sf river.c river-new.c
river-new.s: CFLAGS += -Dswitch_solver=1
river-new.tst: CFLAGS += -Dswitch_solver=1

river.tst: river-new.tst


