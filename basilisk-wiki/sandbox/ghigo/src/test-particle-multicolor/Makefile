all: embed-p-multicolor-2D-tests

include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3 -DDISPLAY=1

# The following tests have special dependencies/compilation requirements

dump-restore.tst: CC = mpicc -D_MPI=4

couette-multicolor.tst: CFLAGS += -grid=multigrid

couette-multicolor-adapt.tst: CFLAGS += -grid=quadtree
couette-multicolor-adapt.c: couette-multicolor.c
	ln -s couette-multicolor.c couette-multicolor-adapt.c

tricylinder-steady.tst: CFLAGS += -grid=multigrid

tricylinder-steady-adapt.tst: CFLAGS += -grid=quadtree
tricylinder-steady-adapt.c: tricylinder-steady.c
	ln -s tricylinder-steady.c tricylinder-steady-adapt.c

tricylinder-steady-moving.tst: CFLAGS += -grid=multigrid

tricylinder-steady-moving-adapt.tst: CFLAGS += -grid=quadtree
tricylinder-steady-moving-adapt.c: tricylinder-steady-moving.c
	ln -s tricylinder-steady-moving.c tricylinder-steady-moving-adapt.c

tricylinder-buoyant.tst: CFLAGS += -grid=multigrid

tricylinder-buoyant-adapt.tst: CFLAGS += -grid=quadtree
tricylinder-buoyant-adapt.c: tricylinder-buoyant.c
	ln -s tricylinder-buoyant.c tricylinder-buoyant-adapt.c

tricylinder-buoyant1.tst: CFLAGS += -grid=multigrid

tricylinder-buoyant1-adapt.tst: CFLAGS += -grid=quadtree
tricylinder-buoyant1-adapt.c: tricylinder-buoyant1.c
	ln -s tricylinder-buoyant1.c tricylinder-buoyant1-adapt.c

taylor-green-quadricylinder.tst: CFLAGS += -grid=multigrid

cylinder-confined-settling1.tst: CFLAGS += -DDENSITY=1
cylinder-confined-settling1.c: cylinder-confined-settling.c
	ln -s cylinder-confined-settling.c cylinder-confined-settling1.c

# bicylinders-pure-wake-interaction.tst: CC = mpicc -D_MPI=6
bicylinders-pure-wake-interaction.tst: CFLAGS += -DADAPTOMEGA=1

embed-p-color-2D-tests: couette-multicolor.tst couette-multicolor-adapt.tst tricylinder-steady.tst tricylinder-steady-moving.tst tricylinder-buoyant.tst tricylinder-buoyant-adapt.tst tricylinder-buoyant1.tst tricylinder-buoyant1-adapt.tst taylor-green-quadricylinder.tst cylinder-confined-settling.tst cylinder-confined-settling1.tst bicylinders-pure-wake-interaction.tst bicylinders-drafting-kissing-tumbling.tst cylinder-impacting-wall.tst
