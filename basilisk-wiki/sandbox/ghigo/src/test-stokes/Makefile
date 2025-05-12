all:	embed-stokes-2D-fixed-tests \
	embed-stokes-2D-moving-tests \
	embed-stokes-3D-fixed-tests \
	embed-stokes-3D-moving-tests

include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -Wdimensions -DMTRACE=3 -DDISPLAY=1

# Check

#embed-stokes-check-tests: CFLAGS += -catch -DTRASH=1
embed-stokes-check-tests:	embed-stokes-2D-fixed-tests \
				porous3D.tst \
				spheres.tst \
				spheres1.tst \
				spheres2.tst \
				torque.tst \
				embed-stokes-2D-moving-tests

# Minimal tests

embed-stokes-minimal-tests:	uf.tst uf-noavg.tst \
				porous.tst porous1.tst porous3D.tst \
				poiseuille.tst poiseuille-adapt.tst \
				poiseuille45.tst poiseuille45-adapt.tst \
				couette.tst couette-adapt.tst \
				wannier.tst wannier-adapt.tst \
				cylinders.tst spheres.tst \
				cylinder-steady.tst cylinder-steady-adapt.tst \
				cylinder-steady-moving.tst cylinder-steady-moving-adapt.tst \
				torque.tst \
				sphere-towards-wall.tst

# The following tests have special dependencies/compilation requirements

## 2D

uf-noavg.tst: CFLAGS += -DNOAVG=1
uf-noavg.c: uf.c
	ln -s uf.c uf-noavg.c

couette.tst: CFLAGS += -grid=multigrid

couette-adapt.tst: CFLAGS += -grid=quadtree
couette-adapt.c: couette.c
	ln -s couette.c couette-adapt.c

couette1.tst: CFLAGS += -grid=multigrid -DCENTERED=1
couette1.c: couette.c
	ln -s couette.c couette1.c

wannier.tst: CFLAGS += -grid=multigrid

wannier-adapt.tst: CFLAGS += -grid=quadtree
wannier-adapt.c: wannier.c
	ln -s wannier.c wannier-adapt.c

wannier1.tst: CFLAGS += -grid=multigrid -DCENTERED=1
wannier1.c: wannier.c
	ln -s wannier.c wannier1.c

poiseuille45.tst: CFLAGS += -grid=multigrid

poiseuille45-adapt.tst: CFLAGS += -grid=quadtree
poiseuille45-adapt.c: poiseuille45.c
	ln -s poiseuille45.c poiseuille45-adapt.c

poiseuille451.tst: CFLAGS += -grid=multigrid -DCENTERED=1
poiseuille451.c: poiseuille45.c
	ln -s poiseuille45.c poiseuille451.c

poiseuille451-adapt.tst: CFLAGS += -grid=quadtree -DCENTERED=1
poiseuille451-adapt.c: poiseuille45.c
	ln -s poiseuille45.c poiseuille451-adapt.c

poiseuille.tst: CFLAGS += -grid=multigrid

poiseuille-adapt.tst: CFLAGS += -grid=quadtree
poiseuille-adapt.c: poiseuille.c
	ln -s poiseuille.c poiseuille-adapt.c

poiseuille1.tst: CFLAGS += -grid=multigrid -DCENTERED=1
poiseuille1.c: poiseuille.c
	ln -s poiseuille.c poiseuille1.c

porous.tst: CFLAGS += -DORDER2=1

porous1.tst: CFLAGS += -DORDER2=1

cylinders.tst: CFLAGS += -DORDER2=1

cylinder-steady.tst: CFLAGS += -grid=multigrid

cylinder-steady-o2.tst: CFLAGS += -grid=multigrid -DORDER2=1
cylinder-steady-o2.c: cylinder-steady.c
	ln -s cylinder-steady.c cylinder-steady-o2.c

cylinder-steady-adapt.tst: CFLAGS += -grid=quadtree
cylinder-steady-adapt.c: cylinder-steady.c
	ln -s cylinder-steady.c cylinder-steady-adapt.c

cylinder-steady-adapt-o2.tst: CFLAGS += -grid=quadtree -DORDER2=1
cylinder-steady-adapt-o2.c: cylinder-steady.c
	ln -s cylinder-steady.c cylinder-steady-adapt-o2.c

embed-stokes-2D-fixed-tests: uf.tst uf-noavg.tst couette.tst couette-adapt.tst wannier.tst wannier-adapt.tst poiseuille45.tst poiseuille45-adapt.tst porous.tst porous1.tst cylinders.tst cylinder-steady.tst cylinder-steady-o2.tst cylinder-steady-adapt.tst cylinder-steady-adapt-o2.tst

## 3D

porous3D.tst: CC = mpicc -D_MPI=4
porous3D.tst: CFLAGS += -DORDER2=1

spheres.tst: CC = mpicc -D_MPI=8
spheres.tst: CFLAGS += -DLVL=5

spheres1.tst: CC = mpicc -D_MPI=8
spheres1.tst: CFLAGS += -DLVL=6
spheres1.c: spheres.c
	ln -s spheres.c spheres1.c

spheres2.tst: CC = mpicc -D_MPI=8
spheres2.tst: CFLAGS += -DLVL=7
spheres2.c: spheres.c
	ln -s spheres.c spheres2.c

torque.tst: CC = mpicc -D_MPI=4

sphere-towards-wall.tst: CC = mpicc -D_MPI=4
sphere-towards-wall.tst: CFLAGS += -DTRACE=2 -DDUMP=1 -DGAP=300 -DDTMAX=10 -DDLENGTH=128 -DLMAX=12

embed-stokes-3D-fixed-tests: porous3D.tst spheres.tst spheres1.tst spheres2.tst torque.tst sphere-towards-wall.tst

## 2D

cylinder-steady-moving.tst: CFLAGS += -grid=multigrid

cylinder-steady-moving-o2.tst: CFLAGS += -grid=multigrid -DORDER2=1
cylinder-steady-moving-o2.c: cylinder-steady-moving.c
	ln -s cylinder-steady-moving.c cylinder-steady-moving-o2.c

cylinder-steady-moving-adapt.tst: CFLAGS += -grid=quadtree
cylinder-steady-moving-adapt.c: cylinder-steady-moving.c
	ln -s cylinder-steady-moving.c cylinder-steady-moving-adapt.c

cylinder-steady-moving-adapt-o2.tst: CFLAGS += -grid=quadtree -DORDER2=1
cylinder-steady-moving-adapt-o2.c: cylinder-steady-moving.c
	ln -s cylinder-steady-moving.c cylinder-steady-moving-adapt-o2.c

embed-stokes-2D-moving-tests: cylinder-steady-moving.tst cylinder-steady-moving-o2.tst cylinder-steady-moving-adapt.tst cylinder-steady-moving-adapt-o2.tst
