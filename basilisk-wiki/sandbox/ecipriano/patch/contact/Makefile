include $(BASILISK)/Makefile.defs

CFLAGS += -O2

sessile-gerris-small.c: sessile-gerris.c
	ln -s sessile-gerris.c sessile-gerris-small.c
sessile-gerris-small.s: CFLAGS += -DSMALL_ANGLES=1
sessile-gerris-small.tst: CFLAGS += -DSMALL_ANGLES=1

sessile-gerris.tst: sessile-gerris-small.tst
