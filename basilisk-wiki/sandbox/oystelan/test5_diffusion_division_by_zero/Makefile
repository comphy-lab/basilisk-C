# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -fopenmp -Wdimensions -DMTRACE=3 -I$(DOCUMENT_ROOT)/sandbox

include $(BASILISK)/Makefile.defs

division_by_zero-new.c: division_by_zero.c
	ln -sf division_by_zero.c division_by_zero-new.c
division_by_zero-new.s: CFLAGS += -Dswitch_solver=1
division_by_zero-new.tst: CFLAGS += -Dswitch_solver=1

division_by_zero.tst: division_by_zero-new.tst

#clean:
#	rm -R division_by_zero \
#	      division_by_zero-new
