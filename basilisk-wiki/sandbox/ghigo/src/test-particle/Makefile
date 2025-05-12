all: 	embed-p-2D-tests \
	embed-p-3D-tests

include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3 -DDISPLAY=1

# Check

# embed-p-check-tests: CFLAGS += -catch -DTRASH=1
embed-p-check-tests: 	embed-p-2D-tests \
			sphere-advection-pressure-gradient.tst

# Minimal tests

embed-p-minimal-tests: 	cylinder-buoyant.tst cylinder-buoyant1.tst \
			cylinder-advection-pressure-gradient.tst sphere-advection-pressure-gradient.tst \
			sphere-settling.tst sphere-settling1.tst sphere-settling2.tst sphere-settling3.tst

# The following tests have special dependencies/compilation requirements

## 2D

cylinder-buoyant.tst: CFLAGS += -grid=multigrid

cylinder-buoyant-adapt.tst: CFLAGS += -grid=quadtree
cylinder-buoyant-adapt.c: cylinder-buoyant.c
	ln -s cylinder-buoyant.c cylinder-buoyant-adapt.c

cylinder-buoyant1.tst: CFLAGS += -grid=multigrid

cylinder-buoyant1-adapt.tst: CFLAGS += -grid=quadtree
cylinder-buoyant1-adapt.c: cylinder-buoyant1.c
	ln -s cylinder-buoyant1.c cylinder-buoyant1-adapt.c

cylinder-advection-pressure-gradient.tst: CFLAGS += -grid=multigrid

cylinder-advection-pressure-gradient-adapt.tst: CFLAGS += -grid=quadtree
cylinder-advection-pressure-gradient-adapt.c: cylinder-advection-pressure-gradient.c
	ln -s cylinder-advection-pressure-gradient.c cylinder-advection-pressure-gradient-adapt.c

embed-p-2D-tests: cylinder-buoyant.tst cylinder-buoyant-adapt.tst cylinder-buoyant1.tst cylinder-buoyant1-adapt.tst cylinder-advection-pressure-gradient.tst cylinder-advection-pressure-gradient-adapt.tst cylinder-impacting-wall.tst

## 3D

sphere-advection-pressure-gradient.tst: CC = mpicc -D_MPI=4

sphere-settling.tst: CC = mpicc -D_MPI=4

sphere-settling1.tst: CC = mpicc -D_MPI=4
sphere-settling1.tst: CFLAGS += -DRE=1
sphere-settling1.c: sphere-settling.c
	ln -s sphere-settling.c sphere-settling1.c

sphere-settling2.tst: CC = mpicc -D_MPI=4
sphere-settling2.tst: CFLAGS += -DRE=2
sphere-settling2.c: sphere-settling.c
	ln -s sphere-settling.c sphere-settling2.c

sphere-settling3.tst: CC = mpicc -D_MPI=4
sphere-settling3.tst: CFLAGS += -DRE=3
sphere-settling3.c: sphere-settling.c
	ln -s sphere-settling.c sphere-settling3.c

sphere-settling-pbc.tst: CC = mpicc -D_MPI=4
sphere-settling-pbc.tst: CFLAGS += -DPNEUMANN=1
sphere-settling-pbc.c: sphere-settling.c
	ln -s sphere-settling.c sphere-settling-pbc.c

sphere-settling1-pbc.tst: CC = mpicc -D_MPI=4
sphere-settling1-pbc.tst: CFLAGS += -DRE=1 -DPNEUMANN=1
sphere-settling1-pbc.c: sphere-settling.c
	ln -s sphere-settling.c sphere-settling1-pbc.c

sphere-settling2-pbc.tst: CC = mpicc -D_MPI=4
sphere-settling2-pbc.tst: CFLAGS += -DRE=2 -DPNEUMANN=1
sphere-settling2-pbc.c: sphere-settling.c
	ln -s sphere-settling.c sphere-settling2-pbc.c

sphere-settling3-pbc.tst: CC = mpicc -D_MPI=4
sphere-settling3-pbc.tst: CFLAGS += -DRE=3 -DPNEUMANN=1
sphere-settling3-pbc.c: sphere-settling.c
	ln -s sphere-settling.c sphere-settling3-pbc.c

sphere-settling-large-domain.tst: CC = mpicc -D_MPI=4

sphere-settling-large-domain1.tst: CC = mpicc -D_MPI=4
sphere-settling-large-domain1.tst: CFLAGS += -DRE=1
sphere-settling-large-domain1.c: sphere-settling-large-domain.c
	ln -s sphere-settling-large-domain.c sphere-settling-large-domain1.c

sphere-settling-large-domain2.tst: CC = mpicc -D_MPI=4
sphere-settling-large-domain2.tst: CFLAGS += -DRE=2
sphere-settling-large-domain2.c: sphere-settling-large-domain.c
	ln -s sphere-settling-large-domain.c sphere-settling-large-domain2.c

embed-p-3D-tests: sphere-advection-pressure-gradient.tst sphere-settling.tst sphere-settling1.tst sphere-settling2.tst sphere-settling3.tst sphere-settling-large-domain.tst sphere-settling-large-domain1.tst sphere-settling-large-domain2.tst
