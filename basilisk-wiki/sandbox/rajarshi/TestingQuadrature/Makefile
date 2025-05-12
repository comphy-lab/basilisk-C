# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3

include $(BASILISK)/Makefile.defs

TestingQuadratureYes.c: TestingQuadrature.c
	ln -s TestingQuadrature.c TestingQuadratureYes.c
TestingQuadratureYes.tst: CFLAGS += -DQuad=1

TestingQuadrature.tst: TestingQuadratureYes.tst
