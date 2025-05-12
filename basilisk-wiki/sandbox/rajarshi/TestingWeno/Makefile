/**
#Makefile
*/

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3

include $(BASILISK)/Makefile.defs

TestingWenoScheme-quadrature.c: TestingWenoScheme.c
	ln -s TestingWenoScheme.c TestingWenoScheme-quadrature.c
TestingWenoScheme-quadrature.tst: CFLAGS += -DQuad=1

TestingWenoScheme.tst:  TestingWenoScheme-quadrature.tst