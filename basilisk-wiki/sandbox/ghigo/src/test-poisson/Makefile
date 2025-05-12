all: 	embed-poisson-2D-neumann-tests \
	embed-poisson-2D-dirichlet-tests \
	embed-poisson-3D-neumann-tests \
	embed-poisson-3D-dirichlet-tests \
	embed-extrapolation-tests \
	embed-boundary-condition.tst

include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3

# Check

embed-poisson-check-tests: CFLAGS += -catch -DTRASH=1
embed-poisson-check-tests: 	embed-poisson-2D-neumann-tests \
			  	embed-poisson-2D-dirichlet-tests \
				embed-poisson-3D-neumann-tests \
				embed-poisson-3D-dirichlet-tests \
				embed-extrapolation-tests \
				embed-boundary-condition.tst

# Minimal tests

embed-poisson-minimal-tests:	neumann.tst \
				neumann3D.tst \
				dirichlet.tst \
				dirichlet3D.tst \
				extrapolation.tst \
				extrapolation3D.tst \
				embed-boundary-condition.tst

# The following tests have special dependencies/compilation requirements

## 2D Dirichlet

dirichlet1.tst: CFLAGS += -DDIRICHLET=1 -DGEOM=1 -grid=multigrid
dirichlet1.c: neumann.c
	ln -s neumann.c dirichlet1.c

dirichlet1-o2.tst: CFLAGS += -DDIRICHLET=1 -DGEOM=1 -grid=multigrid -DORDER2=1
dirichlet1-o2.c: neumann.c
	ln -s neumann.c dirichlet1-o2.c

dirichlet1-adapt.tst: CFLAGS += -DDIRICHLET=1 -DGEOM=1 -grid=quadtree
dirichlet1-adapt.c: neumann.c
	ln -s neumann.c dirichlet1-adapt.c

dirichlet2.tst: CFLAGS += -DDIRICHLET=1 -DGEOM=2 -grid=multigrid
dirichlet2.c: neumann.c
	ln -s neumann.c dirichlet2.c

dirichlet2-o2.tst: CFLAGS += -DDIRICHLET=1 -DGEOM=2 -grid=multigrid -DORDER2=1
dirichlet2-o2.c: neumann.c
	ln -s neumann.c dirichlet2-o2.c

dirichlet2-adapt.tst: CFLAGS += -DDIRICHLET=1 -DGEOM=2 -grid=quadtree
dirichlet2-adapt.c: neumann.c
	ln -s neumann.c dirichlet2-adapt.c

dirichlet3.tst: CFLAGS += -DDIRICHLET=1 -DGEOM=3 -grid=multigrid
dirichlet3.c: neumann.c
	ln -s neumann.c dirichlet3.c

dirichlet3-o2.tst: CFLAGS += -DDIRICHLET=1 -DGEOM=3 -grid=multigrid -DORDER2=1
dirichlet3-o2.c: neumann.c
	ln -s neumann.c dirichlet3-o2.c

dirichlet3-adapt.tst: CFLAGS += -DDIRICHLET=1 -DGEOM=3 -grid=quadtree
dirichlet3-adapt.c: neumann.c
	ln -s neumann.c dirichlet3-adapt.c

dirichlet.tst: CFLAGS += -DDIRICHLET=1 -grid=multigrid
dirichlet.c: neumann.c
	ln -s neumann.c dirichlet.c

dirichlet-o2.tst: CFLAGS += -DDIRICHLET=1 -grid=multigrid -DORDER2=1
dirichlet-o2.c: neumann.c
	ln -s neumann.c dirichlet-o2.c

dirichlet-adapt.tst: CFLAGS += -DDIRICHLET=1 -grid=quadtree
dirichlet-adapt.c: neumann.c
	ln -s neumann.c dirichlet-adapt.c

dirichlet-adapt-ng.tst: CFLAGS += -DDIRICHLET=1 -grid=quadtree -DNGRAD=1
dirichlet-adapt-ng.c: neumann.c
	ln -s neumann.c dirichlet-adapt-ng.c

embed-poisson-2D-dirichlet-tests: dirichlet.tst dirichlet-o2.tst dirichlet-adapt.tst dirichlet-adapt-ng.tst dirichlet1.tst dirichlet1-o2.tst dirichlet1-adapt.tst dirichlet2.tst dirichlet2-o2.tst dirichlet2-adapt.tst dirichlet3.tst dirichlet3-o2.tst dirichlet3-adapt.tst

## 2D Neumann

neumann1.tst: CFLAGS += -DGEOM=1 -grid=multigrid
neumann1.c: neumann.c
	ln -s neumann.c neumann1.c

neumann1-o2.tst: CFLAGS += -DGEOM=1 -grid=multigrid -DORDER2=1
neumann1-o2.c: neumann.c
	ln -s neumann.c neumann1-o2.c

neumann1-adapt.tst: CFLAGS += -DGEOM=1 -grid=quadtree
neumann1-adapt.c: neumann.c
	ln -s neumann.c neumann1-adapt.c

neumann2.tst: CFLAGS += -DGEOM=2 -grid=multigrid
neumann2.c: neumann.c
	ln -s neumann.c neumann2.c

neumann2-o2.tst: CFLAGS += -DGEOM=2 -grid=multigrid -DORDER2=1
neumann2-o2.c: neumann.c
	ln -s neumann.c neumann2-o2.c

neumann2-adapt.tst: CFLAGS += -DGEOM=2 -grid=quadtree
neumann2-adapt.c: neumann.c
	ln -s neumann.c neumann2-adapt.c

neumann3.tst: CFLAGS += -DGEOM=3 -grid=multigrid
neumann3.c: neumann.c
	ln -s neumann.c neumann3.c

neumann3-o2.tst: CFLAGS += -DGEOM=3 -grid=multigrid -DORDER2=1
neumann3-o2.c: neumann.c
	ln -s neumann.c neumann3-o2.c

neumann3-adapt.tst: CFLAGS += -DGEOM=3 -grid=quadtree
neumann3-adapt.c: neumann.c
	ln -s neumann.c neumann3-adapt.c

neumann.tst: CFLAGS += -grid=multigrid

neumann-o2.tst: CFLAGS += -grid=multigrid -DORDER2=1
neumann-o2.c: neumann.c
	ln -s neumann.c neumann-o2.c

neumann-adapt.tst: CFLAGS += -grid=quadtree
neumann-adapt.c: neumann.c
	ln -s neumann.c neumann-adapt.c

embed-poisson-2D-neumann-tests: neumann.tst neumann-o2.tst neumann-adapt.tst neumann1.tst neumann1-o2.tst neumann1-adapt.tst neumann2.tst neumann2-o2.tst neumann2-adapt.tst neumann3.tst neumann3-o2.tst neumann3-adapt.tst

## 3D Dirichlet

dirichlet3D.tst: CFLAGS += -DDIRICHLET=1 -grid=multigrid3D
dirichlet3D.c: neumann3D.c
	ln -s neumann3D.c dirichlet3D.c

dirichlet3D-adapt.tst: CFLAGS += -DDIRICHLET=1 -grid=octree
dirichlet3D-adapt.c: neumann3D.c
	ln -s neumann3D.c dirichlet3D-adapt.c

dirichlet3D1.tst: CFLAGS += -DGEOM=1 -DDIRICHLET=1 -grid=multigrid3D
dirichlet3D1.c: neumann3D.c
	ln -s neumann3D.c dirichlet3D1.c

dirichlet3D1-adapt.tst: CFLAGS += -DGEOM=1 -DDIRICHLET=1 -grid=octree
dirichlet3D1-adapt.c: neumann3D.c
	ln -s neumann3D.c dirichlet3D1-adapt.c

embed-poisson-3D-dirichlet-tests: dirichlet3D.tst dirichlet3D-adapt.tst dirichlet3D1.tst dirichlet3D1-adapt.tst

## 3D Neumann

neumann3D.tst: CFLAGS += -grid=multigrid3D

neumann3D-adapt.tst: CFLAGS += -grid=octree
neumann3D-adapt.c: neumann3D.c
	ln -s neumann3D.c neumann3D-adapt.c

neumann3D1.tst: CFLAGS += -DGEOM=1 -grid=multigrid3D
neumann3D1.c: neumann3D.c	
	ln -s neumann3D.c neumann3D1.c

neumann3D1-adapt.tst: CFLAGS += -DGEOM=1 -grid=octree
neumann3D1-adapt.c: neumann3D.c
	ln -s neumann3D.c neumann3D1-adapt.c

embed-poisson-3D-neumann-tests: neumann3D.tst neumann3D-adapt.tst neumann3D1.tst neumann3D1-adapt.tst

## 2D extrapolation

extrapolation.tst: CFLAGS += -grid=multigrid

extrapolation-adapt.tst: CFLAGS += -grid=quadtree
extrapolation-adapt.c: extrapolation.c
	ln -s extrapolation.c extrapolation-adapt.c

extrapolation-ls.tst: CFLAGS += -DLS=1 -grid=multigrid
extrapolation-ls.c: extrapolation.c
	ln -s extrapolation.c extrapolation-ls.c

extrapolation-ls-adapt.tst: CFLAGS += -DLS=1 -grid=quadtree
extrapolation-ls-adapt.c: extrapolation.c
	ln -s extrapolation.c extrapolation-ls-adapt.c

extrapolation1.tst: CFLAGS += -DGEOM=1 -grid=multigrid
extrapolation1.c: extrapolation.c
	ln -s extrapolation.c extrapolation1.c

extrapolation1-adapt.tst: CFLAGS += -DGEOM=1 -grid=quadtree
extrapolation1-adapt.c: extrapolation.c
	ln -s extrapolation.c extrapolation1-adapt.c

extrapolation1-ls.tst: CFLAGS += -DLS=1 -DGEOM=1 -grid=multigrid
extrapolation1-ls.c: extrapolation.c
	ln -s extrapolation.c extrapolation1-ls.c

extrapolation1-ls-adapt.tst: CFLAGS += -DLS=1 -DGEOM=1 -grid=quadtree
extrapolation1-ls-adapt.c: extrapolation.c
	ln -s extrapolation.c extrapolation1-ls-adapt.c

extrapolation2.tst: CFLAGS += -DGEOM=2 -grid=multigrid
extrapolation2.c: extrapolation.c
	ln -s extrapolation.c extrapolation2.c

extrapolation2-adapt.tst: CFLAGS += -DGEOM=2 -grid=quadtree
extrapolation2-adapt.c: extrapolation.c
	ln -s extrapolation.c extrapolation2-adapt.c

extrapolation2-ls.tst: CFLAGS += -DLS=1 -DGEOM=2 -grid=multigrid
extrapolation2-ls.c: extrapolation.c
	ln -s extrapolation.c extrapolation2-ls.c

extrapolation2-ls-adapt.tst: CFLAGS += -DLS=1 -DGEOM=2 -grid=quadtree
extrapolation2-ls-adapt.c: extrapolation.c
	ln -s extrapolation.c extrapolation2-ls-adapt.c

extrapolation3.tst: CFLAGS += -DGEOM=3 -grid=multigrid
extrapolation3.c: extrapolation.c
	ln -s extrapolation.c extrapolation3.c

extrapolation3-adapt.tst: CFLAGS += -DGEOM=3 -grid=quadtree
extrapolation3-adapt.c: extrapolation.c
	ln -s extrapolation.c extrapolation3-adapt.c

extrapolation3-ls.tst: CFLAGS += -DLS=1 -DGEOM=3 -grid=multigrid
extrapolation3-ls.c: extrapolation.c
	ln -s extrapolation.c extrapolation3-ls.c

extrapolation3-ls-adapt.tst: CFLAGS += -DLS=1 -DGEOM=3 -grid=quadtree
extrapolation3-ls-adapt.c: extrapolation.c
	ln -s extrapolation.c extrapolation3-ls-adapt.c

## 3D extrapolation

extrapolation3D.tst: CFLAGS += -grid=multigrid3D

extrapolation3D-adapt.tst: CFLAGS += -grid=octree
extrapolation3D-adapt.c: extrapolation3D.c
	ln -s extrapolation3D.c extrapolation3D-adapt.c

extrapolation3D-ls.tst: CFLAGS += -DLS=1 -grid=multigrid3D
extrapolation3D-ls.c: extrapolation3D.c
	ln -s extrapolation3D.c extrapolation3D-ls.c

extrapolation3D-ls-adapt.tst: CFLAGS += -DLS=1 -grid=octree
extrapolation3D-ls-adapt.c: extrapolation3D.c
	ln -s extrapolation3D.c extrapolation3D-ls-adapt.c

extrapolation3D1.tst: CFLAGS += -DGEOM=1 -grid=multigrid3D
extrapolation3D1.c: extrapolation3D.c
	ln -s extrapolation3D.c extrapolation3D1.c

extrapolation3D1-adapt.tst: CFLAGS += -DGEOM=1 -grid=octree
extrapolation3D1-adapt.c: extrapolation3D.c
	ln -s extrapolation3D.c extrapolation3D1-adapt.c

extrapolation3D1-ls.tst: CFLAGS += -DLS=1 -DGEOM=1 -grid=multigrid3D
extrapolation3D1-ls.c: extrapolation3D.c
	ln -s extrapolation3D.c extrapolation3D1-ls.c

extrapolation3D1-ls-adapt.tst: CFLAGS += -DLS=1 -DGEOM=1 -grid=octree
extrapolation3D1-ls-adapt.c: extrapolation3D.c
	ln -s extrapolation3D.c extrapolation3D1-ls-adapt.c

embed-extrapolation-tests: extrapolation.tst extrapolation-adapt.tst extrapolation1.tst extrapolation1-adapt.tst extrapolation2.tst extrapolation2-adapt.tst extrapolation3.tst extrapolation3-adapt.tst extrapolation-ls.tst extrapolation-ls-adapt.tst extrapolation1-ls.tst extrapolation1-ls-adapt.tst extrapolation2-ls.tst extrapolation2-ls-adapt.tst extrapolation3-ls.tst extrapolation3-ls-adapt.tst extrapolation3D.tst extrapolation3D-adapt.tst extrapolation3D1.tst extrapolation3D1-adapt.tst extrapolation3D-ls.tst extrapolation3D-ls-adapt.tst extrapolation3D1-ls.tst extrapolation3D1-ls-adapt.tst
