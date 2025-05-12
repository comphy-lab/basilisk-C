# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3

include $(BASILISK)/Makefile.defs

Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-weno.c: Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-bcg.c
	ln -s Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-bcg.c Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-weno.c
Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-weno.tst: CFLAGS += -DWENO=1 -DLIMITED=1

Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-bcg.tst: Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-weno.tst

Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-O5.c: Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-bcg.c
	ln -s Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-bcg.c Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-O5.c
Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-O5.tst: CFLAGS += -DWENO=1

Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-bcg.tst: Test_Case-Advection1DPeriodic-bcg_O5_weno-ver-O5.tst



