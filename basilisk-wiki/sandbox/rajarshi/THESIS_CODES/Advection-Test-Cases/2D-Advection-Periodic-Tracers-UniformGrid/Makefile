# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3

include $(BASILISK)/Makefile.defs

2D-PeriodicAdvection-UniformGrid-ver-weno.c: 2D-PeriodicAdvection-UniformGrid-ver-bcg.c
	ln -s 2D-PeriodicAdvection-UniformGrid-ver-bcg.c 2D-PeriodicAdvection-UniformGrid-ver-weno.c
2D-PeriodicAdvection-UniformGrid-ver-weno.tst: CFLAGS += -DWENO=1

2D-PeriodicAdvection-UniformGrid-ver-bcg.tst: 2D-PeriodicAdvection-UniformGrid-ver-weno.tst
