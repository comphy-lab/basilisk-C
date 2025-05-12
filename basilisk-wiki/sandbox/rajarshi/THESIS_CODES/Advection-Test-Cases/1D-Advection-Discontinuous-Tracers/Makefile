# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3

include $(BASILISK)/Makefile.defs

discontinuousadvection1D-Minmod.c: discontinuousadvection1D-NoLimiter.c
	ln -s discontinuousadvection1D-NoLimiter.c discontinuousadvection1D-Minmod.c
discontinuousadvection1D-Minmod.tst: CFLAGS += -DMINMOD=1

discontinuousadvection1D-NoLimiter.tst: discontinuousadvection1D-Minmod.tst

discontinuousadvection1D-Superbee.c: discontinuousadvection1D-NoLimiter.c
	ln -s discontinuousadvection1D-NoLimiter.c discontinuousadvection1D-Superbee.c
discontinuousadvection1D-Superbee.tst: CFLAGS += -DSUPERBEE=1

discontinuousadvection1D-NoLimiter.tst: discontinuousadvection1D-Superbee.tst

discontinuousadvection1D-weno.c: discontinuousadvection1D-NoLimiter.c
	ln -s discontinuousadvection1D-NoLimiter.c discontinuousadvection1D-weno.c
discontinuousadvection1D-weno.tst: CFLAGS += -DWENO=1 -DLIMITED=1 -DORDER=3

discontinuousadvection1D-NoLimiter.tst: discontinuousadvection1D-weno.tst
