all:	embed-ns-fixed-2D-tests \
	embed-ns-fixed-3D-tests \
	embed-ns-multi-fixed-2D-tests \
	embed-ns-moving-2D-tests \
	embed-ns-moving-3D-tests

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3 -DDISPLAY=1

# Check

embed-ns-check-tests:   hydrostatic2.tst hydrostatic3.tst \
			cylinder-steady.tst cylinder-steady-sl.tst cylinder-steady-o2.tst cylinder-steady-adapt.tst cylinder-steady-adapt-sl.tst cylinder-steady-adapt-o2.tst \
			cylinder-steady-moving.tst cylinder-steady-moving-o2.tst cylinder-steady-moving-sl.tst cylinder-steady-moving-adapt.tst cylinder-steady-moving-adapt-o2.tst cylinder-steady-moving-adapt-sl.tst \
			cylinder-rotating.tst cylinder-accelerating.tst \
			cylinder-oscillating.tst cylinder-oscillating1.tst \
			sphere-rotating.tst sphere-accelerating.tst \
			naca0015-pitching.tst naca2414-starting.tst naca2414-starting3D.tst

# Minimal tests

embed-ns-minimal-tests: hydrostatic2.tst hydrostatic3.tst \
			cylinder-steady.tst cylinder-steady-adapt.tst \
			cylinder-steady-moving.tst cylinder-steady-adapt.tst \
			cylinder-rotating.tst sphere-rotating.tst \
			cylinder-accelerating.tst sphere-accelerating.tst \
			cylinder-oscillating.tst cylinder-oscillating1.tst sphere-oscillating.tst \
			starting.tst starting-moving.tst \
			cylinder-oscillation-vertically.tst

# The following tests have special dependencies/compilation requirements

## 2D

hydrostatic2.tst: CC = mpicc -D_MPI=8

starting-2p.tst: CFLAGS += -DDOUBLEPROJECTION=1
starting-2p.c: starting.c
	ln -s starting.c starting-2p.c

starting-bmr.tst: CFLAGS += -DBMR=1
starting-bmr.c: starting.c
	ln -s starting.c starting-bmr.c

starting-sc1.tst: CFLAGS += -DSMALLCELL=1
starting-sc1.c: starting.c
	ln -s starting.c starting-sc1.c

starting-sc2.tst: CFLAGS += -DSMALLCELL=2
starting-sc2.c: starting.c
	ln -s starting.c starting-sc2.c

starting-sc3.tst: CFLAGS += -DSMALLCELL=3
starting-sc3.c: starting.c
	ln -s starting.c starting-sc3.c

starting1.tst: CC = mpicc -D_MPI=4
starting1.tst: CFLAGS += -DRE=1
starting1.c: starting.c
	ln -s starting.c starting1.c

starting1-2p.tst: CC = mpicc -D_MPI=4
starting1-2p.tst: CFLAGS += -DRE=1 -DDOUBLEPROJECTION=1
starting1-2p.c: starting.c
	ln -s starting.c starting1-2p.c

cylinder-steady.tst: CFLAGS += -grid=multigrid

cylinder-steady-sl.tst: CFLAGS += -grid=multigrid -DSLOPELIMITER=1
cylinder-steady-sl.c: cylinder-steady.c
	ln -s cylinder-steady.c cylinder-steady-sl.c

cylinder-steady-o2.tst: CFLAGS += -grid=multigrid -DORDER2=1
cylinder-steady-o2.c: cylinder-steady.c
	ln -s cylinder-steady.c cylinder-steady-o2.c

cylinder-steady-adapt.tst: CFLAGS += -grid=quadtree
cylinder-steady-adapt.c: cylinder-steady.c
	ln -s cylinder-steady.c cylinder-steady-adapt.c

cylinder-steady-adapt-sl.tst: CFLAGS += -grid=quadtree -DSLOPELIMITER=1
cylinder-steady-adapt-sl.c: cylinder-steady.c
	ln -s cylinder-steady.c cylinder-steady-adapt-sl.c

cylinder-steady-adapt-o2.tst: CFLAGS += -grid=quadtree -DORDER2=1
cylinder-steady-adapt-o2.c: cylinder-steady.c
	ln -s cylinder-steady.c cylinder-steady-adapt-o2.c

cylinder-unbounded1.tst: CFLAGS += -DRE=1
cylinder-unbounded1.c: cylinder-unbounded.c
	ln -s cylinder-unbounded.c cylinder-unbounded1.c

cylinder-unbounded2.tst: CFLAGS += -DRE=2
cylinder-unbounded2.c: cylinder-unbounded.c
	ln -s cylinder-unbounded.c cylinder-unbounded2.c

cylinder-strouhal1.tst: CFLAGS += -DSHAPE=1
cylinder-strouhal1.c: cylinder-strouhal.c
	ln -s cylinder-strouhal.c cylinder-strouhal1.c

cylinder-strouhal2.tst: CFLAGS += -DSHAPE=2
cylinder-strouhal2.c: cylinder-strouhal.c
	ln -s cylinder-strouhal.c cylinder-strouhal2.c

embed-ns-fixed-2D-tests: hydrostatic2.tst starting.tst starting-2p.tst starting-bmr.tst starting-sc1.tst starting-sc2.tst starting-sc3.tst starting1.tst starting1-2p.tst cylinder-steady.tst cylinder-steady-sl.tst cylinder-steady-o2.tst cylinder-steady-adapt.tst cylinder-steady-adapt-sl.tst cylinder-steady-adapt-o2.tst cylinder-unbounded.tst cylinder-unbounded1.tst cylinder-unbounded2.tst cylinder-strouhal.tst cylinder-strouhal1.tst cylinder-strouhal2.tst

## 3D

hydrostatic3.tst: CC = mpicc -D_MPI=8

spheroid-unbounded.tst: CC = mpicc -D_MPI=4
spheroid-unbounded.tst: CFLAGS += -DGAMMA=6 -DPHI=45 -DRE=0 -DDTMAX=1 -DDLENGTH=32 -DLMAX=10 -DCMAX=10

embed-ns-fixed-3D-tests: hydrostatic3.tst spheroid-unbounded.tst

## 2D

cylinder-confined1.tst: CFLAGS += -DRE=1
cylinder-confined1.c: cylinder-confined.c
	ln -s cylinder-confined.c cylinder-confined1.c

embed-ns-multi-fixed-2D-tests: cylinder-confined.tst cylinder-confined1.tst

# Moving embedded boundary tests

## 2D

cylinder-steady-moving.tst: CFLAGS += -grid=multigrid

cylinder-steady-moving-o2.tst: CFLAGS += -grid=multigrid -DORDER2=1
cylinder-steady-moving-o2.c: cylinder-steady-moving.c
	ln -s cylinder-steady-moving.c cylinder-steady-moving-o2.c

cylinder-steady-moving-sl.tst: CFLAGS += -grid=multigrid -DSLOPELIMITER=1
cylinder-steady-moving-sl.c: cylinder-steady-moving.c
	ln -s cylinder-steady-moving.c cylinder-steady-moving-sl.c

cylinder-steady-moving-adapt.tst: CFLAGS += -grid=quadtree
cylinder-steady-moving-adapt.c: cylinder-steady-moving.c
	ln -s cylinder-steady-moving.c cylinder-steady-moving-adapt.c

cylinder-steady-moving-adapt-o2.tst: CFLAGS += -grid=quadtree -DORDER2=1
cylinder-steady-moving-adapt-o2.c: cylinder-steady-moving.c
	ln -s cylinder-steady-moving.c cylinder-steady-moving-adapt-o2.c

cylinder-steady-moving-adapt-sl.tst: CFLAGS += -grid=quadtree -DSLOPELIMITER=1
cylinder-steady-moving-adapt-sl.c: cylinder-steady-moving.c
	ln -s cylinder-steady-moving.c cylinder-steady-moving-adapt-sl.c

starting-moving-bmr.tst: CFLAGS += -DBMR=1
starting-moving-bmr.c: starting-moving.c
	ln -s starting-moving.c starting-moving-bmr.c

starting1-moving.tst: CC = mpicc -D_MPI=4
starting1-moving.tst: CFLAGS += -DRE=1
starting1-moving.c: starting-moving.c
	ln -s starting-moving.c starting1-moving.c

cylinder-oscillating1.tst: CFLAGS += -DRE=1
cylinder-oscillating1.c: cylinder-oscillating.c
	ln -s cylinder-oscillating.c cylinder-oscillating1.c

cylinder-oscillating-vertically-bmr.tst: CFLAGS += -DBMR=1
cylinder-oscillating-vertically-bmr.c: cylinder-oscillating-vertically.c
	ln -s cylinder-oscillating-vertically.c cylinder-oscillating-vertically-bmr.c

naca0015-pitching.tst: CC = mpicc -D_MPI=4

naca2414-starting.tst: CC = mpicc -D_MPI=4

embed-ns-moving-2D-tests: cylinder-steady-moving.tst cylinder-steady-moving-o2.tst cylinder-steady-moving-sl.tst cylinder-steady-moving-adapt.tst cylinder-steady-moving-adapt-o2.tst cylinder-steady-moving-adapt-sl.tst cylinder-rotating.tst cylinder-accelerating.tst starting-moving.tst starting-moving-bmr.tst starting1-moving.tst cylinder-oscillating.tst cylinder-oscillating1.tst cylinder-oscillating-vertically.tst cylinder-oscillating-vertically-bmr.tst naca0015-pitching.tst naca2414-starting.tst

## 3D

sphere-rotating.tst: CC = mpicc -D_MPI=4

sphere-accelerating.tst: CC = mpicc -D_MPI=4

sphere-oscillating.tst: CC = mpicc -D_MPI=4
sphere-oscillating.tst: CFLAGS += -DSNAP=1	

naca2414-starting3D.tst: CC = mpicc -D_MPI=4

embed-ns-moving-3D-tests: sphere-rotating.tst sphere-accelerating.tst sphere-oscillating.tst naca2414-starting3D.tst

include $(BASILISK)/Makefile.defs
