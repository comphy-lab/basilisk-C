all: embed-p-color-2D-tests

include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3 -DDISPLAY=1

# The following tests have special dependencies/compilation requirements

couette-color.tst: CFLAGS += -grid=multigrid

couette-color-adapt.tst: CFLAGS += -grid=quadtree
couette-color-adapt.c: couette-color.c
	ln -s couette-color.c couette-color-adapt.c

cylinder-confined-settling1.tst: CFLAGS += -DDENSITY=1
cylinder-confined-settling1.c: cylinder-confined-settling.c
	ln -s cylinder-confined-settling.c cylinder-confined-settling1.c

cylinder-confined-drifting1.tst: CFLAGS += -DDENSITY=1
cylinder-confined-drifting1.c: cylinder-confined-drifting.c
	ln -s cylinder-confined-drifting.c cylinder-confined-drifting1.c

embed-p-color-2D-tests: couette-color.tst couette-color-adapt.tst cylinder-confined-settling.tst cylinder-confined-settling1.tst cylinder-confined-drifting.tst cylinder-confined-drifting1.tst
