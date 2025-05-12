# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3

include $(BASILISK)/Makefile.defs

Gaussian-three.c: Gaussian-two.c
	ln -s Gaussian-two.c Gaussian-three.c
Gaussian-three.tst: CFLAGS += -DTHREEPOINT=1

Gaussian-two.tst: Gaussian-three.tst


