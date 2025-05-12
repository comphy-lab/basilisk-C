# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3

include $(BASILISK)/Makefile.defs

Test_Case-Prolongation-NoLimiting-ver_higherNL.c: Test_Case-Prolongation-NoLimiting-ver_bilinear.c
	ln -s Test_Case-Prolongation-NoLimiting-ver_bilinear.c Test_Case-Prolongation-NoLimiting-ver_higherNL.c
Test_Case-Prolongation-NoLimiting-ver_higherNL.tst: CFLAGS += -DHIGHER=1

Test_Case-Prolongation-NoLimiting-ver_higherL.c: Test_Case-Prolongation-NoLimiting-ver_bilinear.c
	ln -s Test_Case-Prolongation-NoLimiting-ver_bilinear.c Test_Case-Prolongation-NoLimiting-ver_higherL.c
Test_Case-Prolongation-NoLimiting-ver_higherL.tst: CFLAGS += -DHIGHER=1 -DLIMITED=1

Test_Case-Prolongation-NoLimiting-ver_bilinear.tst: Test_Case-Prolongation-NoLimiting-ver_higherNL.tst
Test_Case-Prolongation-NoLimiting-ver_bilinear.tst: Test_Case-Prolongation-NoLimiting-ver_higherL.tst
