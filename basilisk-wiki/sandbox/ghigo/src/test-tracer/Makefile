all:	embed-tracer-2D-tests \
	embed-tracer-3D-tests

include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3 -DDISPLAY=1

# Minimal tests

embed-tracer-minimal-tests: 	cylinder-tracer.tst cylinder-tracer-moving.tst \
				sphere-tracer.tst sphere-tracer-moving.tst

# The following tests have special dependencies/compilation requirements

## 2D

cylinder-tracer.tst: CFLAGS += -grid=multigrid

cylinder-tracer-neumann.tst: CFLAGS += -DNEUMANN=1 -grid=multigrid
cylinder-tracer-neumann.c: cylinder-tracer.c
	ln -s cylinder-tracer.c cylinder-tracer-neumann.c

cylinder-tracer-moving.tst: CFLAGS += -grid=multigrid

cylinder-tracer-moving-neumann.tst: CFLAGS += -DNEUMANN=1 -grid=multigrid
cylinder-tracer-moving-neumann.c: cylinder-tracer-moving.c
	ln -s cylinder-tracer-moving.c cylinder-tracer-moving-neumann.c

embed-tracer-2D-tests: cylinder-tracer.tst cylinder-tracer-neumann.tst cylinder-tracer-moving.tst cylinder-tracer-moving-neumann.tst

## 3D

sphere-tracer-neumann.tst: CFLAGS += -DNEUMANN=1
sphere-tracer-neumann.c: sphere-tracer.c
	ln -s sphere-tracer.c sphere-tracer-neumann.c

sphere-tracer-moving-neumann.tst: CFLAGS += -DNEUMANN=1
sphere-tracer-moving-neumann.c: sphere-tracer-moving.c
	ln -s sphere-tracer-moving.c sphere-tracer-moving-neumann.c

embed-tracer-3D-tests: sphere-tracer.tst sphere-tracer-neumann.tst sphere-tracer-moving.tst sphere-tracer-moving-neumann.tst
