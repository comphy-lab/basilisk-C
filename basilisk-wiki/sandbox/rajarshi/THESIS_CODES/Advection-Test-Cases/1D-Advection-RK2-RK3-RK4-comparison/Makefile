# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3

include $(BASILISK)/Makefile.defs

1D-PeriodicFuncAdvection_TimeComparison-Ver-RK4.tst: CFLAGS += -DWENO=1

1D-PeriodicFuncAdvection_TimeComparison-Ver-RK2.c: 1D-PeriodicFuncAdvection_TimeComparison-Ver-RK4.c
	ln -s 1D-PeriodicFuncAdvection_TimeComparison-Ver-RK4.c 1D-PeriodicFuncAdvection_TimeComparison-Ver-RK2.c
1D-PeriodicFuncAdvection_TimeComparison-Ver-RK2.tst: CFLAGS += -DORDER=2

1D-PeriodicFuncAdvection_TimeComparison-Ver-SSPRK3.c: 1D-PeriodicFuncAdvection_TimeComparison-Ver-RK4.c
	ln -s 1D-PeriodicFuncAdvection_TimeComparison-Ver-RK4.c 1D-PeriodicFuncAdvection_TimeComparison-Ver-SSPRK3.c
1D-PeriodicFuncAdvection_TimeComparison-Ver-SSPRK3.tst: CFLAGS += -DORDER=3

1D-PeriodicFuncAdvection_TimeComparison-Ver-RK4.tst: CFLAGS += -DORDER=4

1D-PeriodicFuncAdvection_TimeComparison-Ver-RK4.tst: 1D-PeriodicFuncAdvection_TimeComparison-Ver-RK2.tst
1D-PeriodicFuncAdvection_TimeComparison-Ver-RK4.tst: 1D-PeriodicFuncAdvection_TimeComparison-Ver-SSPRK3.tst

