all:	embed-heat-fixed-neumann-tests \
	embed-heat-fixed-dirichlet-tests \
	embed-heat-shrinking-neumann-tests \
	embed-heat-shrinking-dirichlet-tests \
	embed-heat-expanding-neumann-tests \
	embed-heat-expanding-dirichlet-tests

include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3 -DDISPLAY=1

# Check 

embed-heat-check-tests:	embed-heat-fixed-neumann-tests \
			embed-heat-fixed-dirichlet-tests \
			embed-heat-shrinking-neumann-tests \
			embed-heat-shrinking-dirichlet-tests \
			embed-heat-expanding-neumann-tests \
			embed-heat-expanding-dirichlet-tests

# Minimal tests

embed-heat-minimal-tests: 	neumann3D.tst \
				dirichlet3D.tst \
				neumann3D-expanding.tst \
				dirichlet3D-expanding.tst \
				neumann3D-shrinking.tst \
				dirichlet3D-shrinking.tst

# The following tests have special dependencies/compilation requirements

## 3D Neumann

neumann3D.tst: CC = mpicc -D_MPI=8
neumann3D.tst: CFLAGS += -grid=multigrid3D

neumann3D-adapt.tst: CC = mpicc -D_MPI=4
neumann3D-adapt.tst: CFLAGS += -grid=octree
neumann3D-adapt.c: neumann3D.c
	ln -s neumann3D.c neumann3D-adapt.c

neumann3D1.tst: CC = mpicc -D_MPI=8
neumann3D1.tst: CFLAGS += -DGEOM=1 -grid=multigrid3D
neumann3D1.c: neumann3D.c
	ln -s neumann3D.c neumann3D1.c

neumann3D1-adapt.tst: CC = mpicc -D_MPI=4
neumann3D1-adapt.tst: CFLAGS += -DGEOM=1 -grid=octree
neumann3D1-adapt.c: neumann3D.c
	ln -s neumann3D.c neumann3D1-adapt.c

embed-heat-fixed-neumann-tests: neumann3D.tst neumann3D-adapt.tst neumann3D1.tst neumann3D1-adapt.tst

## 3D Dirichlet

dirichlet3D.tst: CC = mpicc -D_MPI=8
dirichlet3D.tst: CFLAGS += -DDIRICHLET=1 -grid=multigrid3D
dirichlet3D.c: neumann3D.c
	ln -s neumann3D.c dirichlet3D.c

dirichlet3D-adapt.tst: CC = mpicc -D_MPI=4
dirichlet3D-adapt.tst: CFLAGS += -DDIRICHLET=1 -grid=octree
dirichlet3D-adapt.c: neumann3D.c
	ln -s neumann3D.c dirichlet3D-adapt.c

dirichlet3D1.tst: CC = mpicc -D_MPI=8
dirichlet3D1.tst: CFLAGS += -DGEOM=1 -DDIRICHLET=1 -grid=multigrid3D
dirichlet3D1.c: neumann3D.c
	ln -s neumann3D.c dirichlet3D1.c

dirichlet3D1-adapt.tst: CC = mpicc -D_MPI=4
dirichlet3D1-adapt.tst: CFLAGS += -DGEOM=1 -DDIRICHLET=1 -grid=octree
dirichlet3D1-adapt.c: neumann3D.c
	ln -s neumann3D.c dirichlet3D1-adapt.c

embed-heat-fixed-dirichlet-tests: dirichlet3D.tst dirichlet3D-adapt.tst dirichlet3D1.tst dirichlet3D1-adapt.tst

## 3D Neumann

neumann3D-shrinking.tst: CC = mpicc -D_MPI=8
neumann3D-shrinking.tst: CFLAGS += -DSHRINKING=1 -grid=multigrid3D
neumann3D-shrinking.c: neumann3D.c
	ln -s neumann3D.c neumann3D-shrinking.c

neumann3D-shrinking-adapt.tst: CC = mpicc -D_MPI=4
neumann3D-shrinking-adapt.tst: CFLAGS += -DSHRINKING=1 -grid=octree
neumann3D-shrinking-adapt.c: neumann3D.c
	ln -s neumann3D.c neumann3D-shrinking-adapt.c

neumann3D1-shrinking.tst: CC = mpicc -D_MPI=8
neumann3D1-shrinking.tst: CFLAGS += -DGEOM=1 -DSHRINKING=1 -grid=multigrid3D
neumann3D1-shrinking.c: neumann3D.c
	ln -s neumann3D.c neumann3D1-shrinking.c

neumann3D1-shrinking-adapt.tst: CC = mpicc -D_MPI=4
neumann3D1-shrinking-adapt.tst: CFLAGS += -DGEOM=1 -DSHRINKING=1 -grid=octree
neumann3D1-shrinking-adapt.c: neumann3D.c
	ln -s neumann3D.c neumann3D1-shrinking-adapt.c

neumann3D1-shrinking-ex1.tst: CC = mpicc -D_MPI=8
neumann3D1-shrinking-ex1.tst: CFLAGS += -DGEOM=1 -DSHRINKING=1 -DEXTRAPOLATE=1 -grid=multigrid3D
neumann3D1-shrinking-ex1.c: neumann3D.c
	ln -s neumann3D.c neumann3D1-shrinking-ex1.c

embed-heat-shrinking-neumann-tests: neumann3D-shrinking.tst neumann3D-shrinking-adapt.tst neumann3D1-shrinking.tst neumann3D1-shrinking-ex1.tst neumann3D1-shrinking-adapt.tst

## 3D Dirichlet

dirichlet3D-shrinking.tst: CC = mpicc -D_MPI=8
dirichlet3D-shrinking.tst: CFLAGS += -DDIRICHLET=1 -DSHRINKING=1 -grid=multigrid3D
dirichlet3D-shrinking.c: neumann3D.c
	ln -s neumann3D.c dirichlet3D-shrinking.c

dirichlet3D-shrinking-adapt.tst: CC = mpicc -D_MPI=4
dirichlet3D-shrinking-adapt.tst: CFLAGS += -DDIRICHLET=1 -DSHRINKING=1 -grid=octree
dirichlet3D-shrinking-adapt.c: neumann3D.c
	ln -s neumann3D.c dirichlet3D-shrinking-adapt.c

dirichlet3D1-shrinking.tst: CC = mpicc -D_MPI=8
dirichlet3D1-shrinking.tst: CFLAGS += -DGEOM=1 -DDIRICHLET=1 -DSHRINKING=1 -grid=multigrid3D
dirichlet3D1-shrinking.c: neumann3D.c
	ln -s neumann3D.c dirichlet3D1-shrinking.c

dirichlet3D1-shrinking-adapt.tst: CC = mpicc -D_MPI=4
dirichlet3D1-shrinking-adapt.tst: CFLAGS += -DGEOM=1 -DDIRICHLET=1 -DSHRINKING=1 -grid=octree
dirichlet3D1-shrinking-adapt.c: neumann3D.c
	ln -s neumann3D.c dirichlet3D1-shrinking-adapt.c

embed-heat-shrinking-dirichlet-tests: dirichlet3D-shrinking.tst dirichlet3D-shrinking-adapt.tst dirichlet3D1-shrinking.tst dirichlet3D1-shrinking-adapt.tst

## 3D Neumann

neumann3D-expanding.tst: CC = mpicc -D_MPI=8
neumann3D-expanding.tst: CFLAGS += -DEXPANDING=1 -grid=multigrid3D
neumann3D-expanding.c: neumann3D.c
	ln -s neumann3D.c neumann3D-expanding.c

neumann3D-expanding-ex1.tst: CC = mpicc -D_MPI=8
neumann3D-expanding-ex1.tst: CFLAGS += -DEXPANDING=1 -DEXTRAPOLATE=1 -grid=multigrid3D
neumann3D-expanding-ex1.c: neumann3D.c
	ln -s neumann3D.c neumann3D-expanding-ex1.c

neumann3D-expanding-adapt.tst: CC = mpicc -D_MPI=4
neumann3D-expanding-adapt.tst: CFLAGS += -DEXPANDING=1 -grid=octree
neumann3D-expanding-adapt.c: neumann3D.c
	ln -s neumann3D.c neumann3D-expanding-adapt.c

neumann3D1-expanding.tst: CC = mpicc -D_MPI=8
neumann3D1-expanding.tst: CFLAGS += -DGEOM=1 -DEXPANDING=1 -grid=multigrid3D
neumann3D1-expanding.c: neumann3D.c
	ln -s neumann3D.c neumann3D1-expanding.c

neumann3D1-expanding-adapt.tst: CC = mpicc -D_MPI=4
neumann3D1-expanding-adapt.tst: CFLAGS += -DGEOM=1 -DEXPANDING=1 -grid=octree
neumann3D1-expanding-adapt.c: neumann3D.c
	ln -s neumann3D.c neumann3D1-expanding-adapt.c

embed-heat-expanding-neumann-tests: neumann3D-expanding.tst neumann3D-expanding-ex1.tst neumann3D-expanding-adapt.tst neumann3D1-expanding.tst neumann3D1-expanding-adapt.tst

## 3D Dirichlet

dirichlet3D-expanding.tst: CC = mpicc -D_MPI=8
dirichlet3D-expanding.tst: CFLAGS += -DDIRICHLET=1 -DEXPANDING=1 -grid=multigrid3D
dirichlet3D-expanding.c: neumann3D.c
	ln -s neumann3D.c dirichlet3D-expanding.c

dirichlet3D-expanding-adapt.tst: CC = mpicc -D_MPI=4
dirichlet3D-expanding-adapt.tst: CFLAGS += -DDIRICHLET=1 -DEXPANDING=1 -grid=octree
dirichlet3D-expanding-adapt.c: neumann3D.c
	ln -s neumann3D.c dirichlet3D-expanding-adapt.c

dirichlet3D1-expanding.tst: CC = mpicc -D_MPI=8
dirichlet3D1-expanding.tst: CFLAGS += -DGEOM=1 -DDIRICHLET=1 -DEXPANDING=1 -grid=multigrid3D
dirichlet3D1-expanding.c: neumann3D.c
	ln -s neumann3D.c dirichlet3D1-expanding.c

dirichlet3D1-expanding-adapt.tst: CC = mpicc -D_MPI=4
dirichlet3D1-expanding-adapt.tst: CFLAGS += -DGEOM=1 -DDIRICHLET=1 -DEXPANDING=1 -grid=octree
dirichlet3D1-expanding-adapt.c: neumann3D.c
	ln -s neumann3D.c dirichlet3D1-expanding-adapt.c

embed-heat-expanding-dirichlet-tests: dirichlet3D-expanding.tst dirichlet3D-expanding-adapt.tst dirichlet3D1-expanding.tst dirichlet3D1-expanding-adapt.tst

