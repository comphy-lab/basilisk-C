include $(BASILISK)/Makefile.defs

#BASILISK_INCLUDES ?= -I../src
#BASILISK_LIBS ?= 

#CFLAGS += -O2 -disable-dimensions $(BASILISK_INCLUDES) $(BASILISK_LIBS)
CFLAGS += -O2 -disable-dimensions -I../src

poisson-gfm.tst: CFLAGS += -DGFM=1
capwave.tst: CFLAGS += -DGFM=1
spurious.tst: CFLAGS += -DGFM=1
oscillation.tst: CFLAGS += -DGFM=1

pinning-axi.c: pinning.c
	ln -s pinning.c pinning-axi.c
pinning-axi.s: CFLAGS += -DAXI=1
pinning-axi.tst: CFLAGS += -DAXI=1

pinning.tst: pinning-axi.tst