# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2

include $(BASILISK)/Makefile.defs

wave.tst: CFLAGS += -fopenmp

non-coalescence.tst: coalescence.tst

dune-unsteady.c: dune.c
	ln -s $< $@
dune-unsteady.tst: dune.s
dune-unsteady.tst: CFLAGS += -DUNSTEADY=1

dune.tst: dune-unsteady.tst
