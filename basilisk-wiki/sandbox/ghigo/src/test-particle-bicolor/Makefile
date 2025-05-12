all: embed-p-bicolor-2D-tests

include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3 -DDISPLAY=1

# The following tests have special dependencies/compilation requirements

dump-restore.tst: CC = mpicc -D_MPI=4

couette-bicolor.tst: CFLAGS += -grid=multigrid

couette-bicolor-adapt.tst: CFLAGS += -grid=quadtree
couette-bicolor-adapt.c: couette-bicolor.c
	ln -s couette-bicolor.c couette-bicolor-adapt.c

bicylinder-steady.tst: CFLAGS += -grid=multigrid

bicylinder-steady-moving.tst: CFLAGS += -grid=multigrid

bicylinder-buoyant.tst: CFLAGS += -grid=multigrid

bicylinder-buoyant-adapt.tst: CFLAGS += -grid=quadtree
bicylinder-buoyant-adapt.c: bicylinder-buoyant.c
	ln -s bicylinder-buoyant.c bicylinder-buoyant-adapt.c

bicylinder-buoyant1.tst: CFLAGS += -grid=multigrid

bicylinder-buoyant1-adapt.tst: CFLAGS += -grid=quadtree
bicylinder-buoyant1-adapt.c: bicylinder-buoyant1.c
	ln -s bicylinder-buoyant1.c bicylinder-buoyant1-adapt.c

taylor-green-bicylinder.tst: CFLAGS += -grid=multigrid

embed-p-color-2D-tests: couette-bicolor.tst couette-bicolor-adapt.tst bicylinder-steady.tst bicylinder-steady-moving.tst bicylinder-buoyant.tst bicylinder-buoyant-adapt.tst bicylinder-buoyant1.tst bicylinder-buoyant1-adapt.tst taylor-green-bicylinder.tst
