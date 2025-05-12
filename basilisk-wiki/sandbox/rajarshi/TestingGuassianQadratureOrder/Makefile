# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3

include $(BASILISK)/Makefile.defs

gaussian-three.c: gaussian.c
	ln -s gaussian.c gaussian-three.c
gaussian-three.tst: CFLAGS += -DTHREEPOINT=1

gaussian.tst: gaussian-three.tst
