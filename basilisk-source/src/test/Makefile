all: check 3D-tests mpi-tests load-balancing curvature-tests	\
	axi-tests cadna-tests embed-tests multilayer-tests	\
	multilayer-stratified gotm-tests dimensions-tests gpu-tests

failed:
	$(MAKE) -k `find . -name fail | sed -e 's|\./||g' -e 's|/fail|.tst|g'`

# the default CFLAGS are set in ../config
CFLAGS += -O2 -DMTRACE=3

# list only non-default tests
check: halo.vtst interpolate.vtst events.vtst faces.vtst coarsen.vtst \
	lid-restore.tst poiseuille.vtst tag.ctst

# these tests have special dependencies/compilation requirements

basilisk.tst: basilisk.gnu
basilisk.tst: CFLAGS += -DDEBUG=1

bleck.tst: bleck.ctst

boundaries3.tst: boundaries3.ctst

boundary_face.tst: boundary_face.ctst boundary_face_cartesian.tst

boundary_face_cartesian.c: boundary_face.c boundary_face_cartesian.ref
	cp boundary_face.c boundary_face_cartesian.c
boundary_face_cartesian.ref: boundary_face.ref
	cp boundary_face.ref boundary_face_cartesian.ref
boundary_face_cartesian.s: CFLAGS += -grid=cartesian
boundary_face_cartesian.tst: CFLAGS += -grid=cartesian

bubble-spherical.c: bubble.c
	ln -sf bubble.c bubble-spherical.c 
bubble-spherical.s: CFLAGS += -grid=multigrid1D -DSPHERICAL=1
bubble-spherical.tst: CFLAGS += -grid=multigrid1D -DSPHERICAL=1

collapse-inviscid.c: collapse.c
	ln -sf collapse.c collapse-inviscid.c 
collapse-inviscid.s: CFLAGS += -DINVISCID=1 
collapse-inviscid.tst: CFLAGS += -DINVISCID=1

collapse-spherical.c: collapse.c
	ln -sf collapse.c collapse-spherical.c 
collapse-spherical.s: CFLAGS += -grid=multigrid1D -DSPHERICAL=1
collapse-spherical.tst: CFLAGS += -grid=multigrid1D -DSPHERICAL=1

collapse.tst: collapse-spherical.tst collapse-inviscid.tst

dry-explicit.c: dry.c
	ln -sf dry.c dry-explicit.c 
dry-explicit.s: CFLAGS += -DEXPLICIT=1 
dry-explicit.tst: CFLAGS += -DEXPLICIT=1

dry.tst: dry-explicit.tst

drybump-explicit.c: drybump.c
	ln -sf drybump.c drybump-explicit.c 
drybump-explicit.s: CFLAGS += -DEXPLICIT=1
drybump-explicit.tst: CFLAGS += -DEXPLICIT=1

drybump.tst: drybump-explicit.tst

drybump2D-implicit.c: drybump2D.c
	ln -sf drybump2D.c drybump2D-implicit.c 
drybump2D-implicit.s: CFLAGS += -DIMPLICIT=1
drybump2D-implicit.tst: CFLAGS += -DIMPLICIT=1

drybump2D.tst: drybump2D-implicit.tst

lake-tr.tst: lake-tr-explicit.tst lake-tr-ml.tst

lake-tr-explicit.c: lake-tr.c
	ln -sf lake-tr.c lake-tr-explicit.c
lake-tr-explicit.s: CFLAGS += -DEXPLICIT=1
lake-tr-explicit.tst: CFLAGS += -DEXPLICIT=1

lake-tr-ml.c: lake-tr.c
	ln -sf lake-tr.c lake-tr-ml.c
lake-tr-ml.s: CFLAGS += -DML=1
lake-tr-ml.tst: CFLAGS += -DML=1

multiriverinflow.tst: multiriverinflow.ctst
multiriverinflow.tst: CFLAGS += -fopenmp
multiriverinflow.ctst: CFLAGS += -fopenmp

parabola-explicit.c: parabola.c
	ln -sf parabola.c parabola-explicit.c 
parabola-explicit.s: CFLAGS += -DEXPLICIT=1
parabola-explicit.tst: CFLAGS += -DEXPLICIT=1

parabola.tst: parabola-explicit.tst parabola-ml.tst

seawallsv.c: seawall.c
	ln -sf seawall.c seawallsv.c 
seawallsv.s: CFLAGS += -DSAINT_VENANT=1
seawallsv.tst: CFLAGS += -DSAINT_VENANT=1

seawall.tst: seawallsv.tst seawall-ml.tst

conicalsv.c: conical.c
	ln -sf conical.c conicalsv.c 
conicalsv.s: CFLAGS += -DSAINT_VENANT=1
conicalsv.tst: CFLAGS += -DSAINT_VENANT=1

conical-implicit.c: conical.c
	ln -sf conical.c conical-implicit.c 
conical-implicit.s: CFLAGS += -DSAINT_VENANT=1 -DIMPLICIT=1
conical-implicit.tst: CFLAGS += -DSAINT_VENANT=1 -DIMPLICIT=1

conical.tst: conicalsv.tst conical-implicit.tst conical-ml.tst

curvature.tst: curvature.3D.tst

explicit.c: implicit.c
	ln -sf implicit.c explicit.c 
explicit.s: CFLAGS += -DEXPLICIT=1
explicit.tst: CFLAGS += -DEXPLICIT=1

explicit-ml.c: implicit.c
	ln -sf implicit.c explicit-ml.c 
explicit-ml.s: CFLAGS += -DML=1 -DEXPLICIT=1
explicit-ml.tst: CFLAGS += -DML=1 -DEXPLICIT=1

implicit-ml.c: implicit.c
	ln -sf implicit.c implicit-ml.c 
implicit-ml.s: CFLAGS += -DML=1
implicit-ml.tst: CFLAGS += -DML=1

implicit.tst: explicit.tst explicit-ml.tst implicit-ml.tst

bore.tst: bore1.tst

explosion3D.tst: CFLAGS=-grid=multigrid3D

explosion3D.c: explosion.c
	ln -sf explosion.c explosion3D.c

explosion.tst: explosion.ctst explosion.3D.tst explosion3D.tst

laplacian.tst: laplacian.ctst

lidmac.c: lid.c
	ln -sf lid.c lidmac.c
lidmac.s: CFLAGS += -DMAC=1
lidmac.s.d: CFLAGS += -DMAC=1
lidmac.tst: CFLAGS += -DMAC=1

lid.tst: lidmac.tst

lid-restore.c: lid.c
	ln -sf lid.c lid-restore.c
lid-restore.dump: lid/dump
	ln -sf lid/dump lid-restore.dump
lid-restore.tst: lid.tst lid-restore.dump

oscillation-compressible.c: oscillation.c
	ln -sf oscillation.c oscillation-compressible.c 
oscillation-compressible.s: CFLAGS += -DCOMPRESSIBLE=1
oscillation-compressible.tst: CFLAGS += -DCOMPRESSIBLE=1

oscillation-momentum.c: oscillation.c
	ln -sf oscillation.c oscillation-momentum.c 
oscillation-momentum.s: CFLAGS += -DMOMENTUM=1
oscillation-momentum.tst: CFLAGS += -DMOMENTUM=1

oscillation.tst: oscillation-momentum.tst oscillation-compressible.tst

ponds-implicit.c: ponds.c
	ln -sf ponds.c ponds-implicit.c
ponds-implicit.s: CFLAGS += -DIMPLICIT=1
ponds-implicit.tst: CFLAGS += -DIMPLICIT=1

ponds.tst: ponds-implicit.tst

poiseuille-fenep.c: poiseuille-oldroydb.c
	ln -sf poiseuille-oldroydb.c poiseuille-fenep.c
poiseuille-fenep.s: CFLAGS += -DFENE_P=1
poiseuille-fenep.tst: CFLAGS += -DFENE_P=1

poiseuille-oldroydb.tst: poiseuille-fenep.tst

poisson.tst: poisson.ctst

rising-axi.c: rising.c
	ln -sf rising.c rising-axi.c
rising-axi.s: CFLAGS += -DAXIS=1
rising-axi.tst: CFLAGS += -DAXIS=1

rising-axi-momentum.c: rising.c
	ln -sf rising.c rising-axi-momentum.c
rising-axi-momentum.s: CFLAGS += -DAXIS=1 -DMOMENTUM=1
rising-axi-momentum.tst: CFLAGS += -DAXIS=1 -DMOMENTUM=1

rising-axi-clsvof.c: rising.c
	ln -sf rising.c rising-axi-clsvof.c
rising-axi-clsvof.s: CFLAGS += -DAXIS=1 -DLEVELSET=1 -DCLSVOF=1
rising-axi-clsvof.tst: CFLAGS += -DAXIS=1 -DLEVELSET=1 -DCLSVOF=1

rising2.c: rising.c
	ln -sf rising.c rising2.c
rising2.s: CFLAGS += -DCASE2=1
rising2.tst: CFLAGS += -DCASE2=1

rising2-levelset.c: rising.c
	ln -sf rising.c rising2-levelset.c
rising2-levelset.s: CFLAGS += -DCASE2=1 -DLEVELSET=1
rising2-levelset.tst: CFLAGS += -DCASE2=1 -DLEVELSET=1

rising2-clsvof.c: rising.c
	ln -sf rising.c rising2-clsvof.c
rising2-clsvof.s: CFLAGS += -DCASE2=1 -DLEVELSET=1 -DCLSVOF=1
rising2-clsvof.tst: CFLAGS += -DCASE2=1 -DLEVELSET=1 -DCLSVOF=1

rising-reduced.c: rising.c
	ln -sf rising.c rising-reduced.c
rising-reduced.s: CFLAGS += -DREDUCED=1
rising-reduced.tst: CFLAGS += -DREDUCED=1

rising-levelset.c: rising.c
	ln -sf rising.c rising-levelset.c
rising-levelset.s: CFLAGS += -DLEVELSET=1
rising-levelset.tst: CFLAGS += -DLEVELSET=1

rising-clsvof.c: rising.c
	ln -sf rising.c rising-clsvof.c
rising-clsvof.s: CFLAGS += -DLEVELSET=1 -DCLSVOF=1
rising-clsvof.tst: CFLAGS += -DLEVELSET=1 -DCLSVOF=1

reversed.tst: reversed.ctst

rotate.tst: rotate.ctst

rt-reduced.c: rt.c
	ln -sf rt.c rt-reduced.c
rt-reduced.s: CFLAGS += -DREDUCED=1
rt-reduced.tst: CFLAGS += -DREDUCED=1

taylor-green-all-mach.c: taylor-green.c
	ln -sf taylor-green.c taylor-green-all-mach.c 
taylor-green-all-mach.s: CFLAGS += -DALL_MACH=1
taylor-green-all-mach.tst: CFLAGS += -DALL_MACH=1

taylor-green.tst: taylor-green-all-mach.tst

view.tst: CC = mpicc -D_MPI=4
view.tst: view.3D.tst

view.3D.tst: CC = mpicc -D_MPI=4

coarsen.vtst: CFLAGS = -DMTRACE=3 -progress

# Axisymmetric tests

axi-tests: axiadvection.tst axi.tst poiseuille-axi.tst \
	rising-axi.tst rising-axi-momentum.tst

# 3D tests

3D-tests: circle.3D.tst curvature.3D.tst hf.3D.tst periodic.3D.tst \
	poisson.3D.tst refineu.3D.tst solenoidal.3D.tst

# Curvature tests

curvature-tests: hf.tst hf3D.tst hf.ctst hf.3D.tst hf1.tst hf-mask.tst \
	curvature.tst curvature.3D.tst

hf-mask.tst: CFLAGS += -DTRASH=1

# Multilayer tests

multilayer-tests: dispersion-gn.tst dispersion.tst \
	beach-ml.tst seawall-ml.tst conical-ml.tst \
	bar-ml.tst parabola-ml.tst lonlat-ml.tst \
	wind-driven-nh.tst wind-driven.tst large.tst \
	gaussian-stvt.tst gaussian-hydro.tst gaussian.tst \
	stokes-ns.tst stokes.tst ponds-ml.tst \
	galilean_invariance.tst \
	stommel-ml.tst explicit-ml.tst implicit-ml.tst lake-tr-ml.tst \
	geo.tst nonlinear-ml.tst tsunami.tst

multilayer-stratified: horn.tst lock.tst overflow.tst kh.tst

overflow.tst: overflow-isopycnal.tst

overflow-isopycnal.c: overflow.c
	ln -sf overflow.c overflow-isopycnal.c 
overflow-isopycnal.s: CFLAGS += -DISOPYCNAL=1
overflow-isopycnal.tst: CFLAGS += -DISOPYCNAL=1

kh.tst: kh-ns.tst kh-hydro.tst

kh-ns.tst: CC = mpicc -D_MPI=8

kh-hydro.c: kh.c
	ln -sf kh.c kh-hydro.c 
kh-hydro.s: CFLAGS += -DHYDRO=1
kh-hydro.tst: CFLAGS += -DHYDRO=1

dispersion-gn.c: dispersion.c
	ln -sf dispersion.c dispersion-gn.c 
dispersion-gn.tst: dispersion-gn.s

dispersion.tst: CFLAGS += -DML=1

bar-ml.c: bar.c
	ln -sf bar.c bar-ml.c 
bar-ml.s: CFLAGS += -DML=1
bar-ml.tst: CFLAGS += -DML=1

wind-driven-nh.c: wind-driven.c
	ln -sf wind-driven.c wind-driven-nh.c 
wind-driven-nh.s: CFLAGS += -DML=1 -DNH=1
wind-driven-nh.tst: CFLAGS += -DML=1 -DNH=1

wind-driven-stvt.c: wind-driven.c
	ln -sf wind-driven.c wind-driven-stvt.c 

wind-driven.s: CFLAGS += -DML=1
wind-driven.tst: CFLAGS += -DML=1

large.tst: large-ns.tst

gaussian-hydro.c: gaussian.c
	ln -sf gaussian.c gaussian-hydro.c 
gaussian-hydro.s: CFLAGS += -DML=1 -DHYDRO=1
gaussian-hydro.tst: CFLAGS += -DML=1 -DHYDRO=1

gaussian-nometric.c: gaussian.c
	ln -sf gaussian.c gaussian-nometric.c 
gaussian-nometric.s: CFLAGS += -DML=1 -DNOMETRIC=1
gaussian-nometric.tst: CFLAGS += -DML=1 -DNOMETRIC=1

gaussian-stvt.c: gaussian.c
	ln -sf gaussian.c gaussian-stvt.c 
gaussian-stvt.tst: gaussian-stvt.s

gaussian.s: CFLAGS += -DML=1
gaussian.tst: CFLAGS += -DML=1

parabola-ml.c: parabola.c
	ln -sf parabola.c parabola-ml.c 
parabola-ml.s: CFLAGS += -DML=1
parabola-ml.tst: CFLAGS += -DML=1

seawall-ml.c: seawall.c
	ln -sf seawall.c seawall-ml.c 
seawall-ml.s: CFLAGS += -DML=1
seawall-ml.tst: CFLAGS += -DML=1

beach-ml.c: beach.c
	ln -sf beach.c beach-ml.c 
beach-ml.s: CFLAGS += -DML=1
beach-ml.tst: CFLAGS += -DML=1

conical-ml.c: conical.c
	ln -sf conical.c conical-ml.c 
conical-ml.s: CFLAGS += -DML=1
conical-ml.tst: CFLAGS += -DML=1

lonlat-ml.c: lonlat.c
	ln -sf lonlat.c lonlat-ml.c 
lonlat-ml.s: CFLAGS += -DML=1
lonlat-ml.tst: CFLAGS += -DML=1

ponds-ml.c: ponds.c
	ln -sf ponds.c ponds-ml.c 
ponds-ml.tst: ponds.s

ponds-ml.s: CFLAGS += -DML=1
ponds-ml.tst: CFLAGS += -DML=1

nonlinear.tst: nonlinear-ml.tst

nonlinear-ml.c: nonlinear.c
	ln -sf nonlinear.c nonlinear-ml.c 
nonlinear-ml.s: CFLAGS += -DML=1
nonlinear-ml.tst: CFLAGS += -DML=1

reductions.tst: CC=mpicc -D_MPI=4

shockwave.tst: CFLAGS += -grid=bitree
shockwave.ctst: CFLAGS += -grid=multigrid1D
shockwave.tst: shockwave.ctst

shrinking-axi.c: shrinking.c
	ln -sf shrinking.c shrinking-axi.c 
shrinking-axi.s: CFLAGS += -DAXIS=1 -grid=quadtree
shrinking-axi.tst: CFLAGS += -DAXIS=1 -grid=quadtree

shrinking.s: CFLAGS += -grid=multigrid1D
shrinking.tst: CFLAGS += -grid=multigrid1D
shrinking.tst: shrinking-axi.tst

## MPI tests

# bump2Dp.tst crashes with -DDEBUGCOND=1 but works with -DDEBUGCOND=0

mpi-tests: indexing.tst indexing.3D.tst \
	mpi-restriction.tst mpi-restriction.3D.tst \
	mpi-reduce.tst openmp-reduce.tst \
	mpi-refine.tst mpi-refine1.tst mpi-refine.3D.tst \
	mpi-laplacian.tst mpi-laplacian.3D.tst \
	mpi-circle.tst mpi-circle1.tst mpi-flux.tst \
	mpi-interpu.tst mpi-coarsen.tst mpi-coarsen1.tst \
	hf1.tst pdump.tst restore.tst \
	pdump-multigrid.tst restore-multigrid.tst \
	restore-tree.tst \
	poiseuille-periodic.tst \
	gfsi.tst gfs.tst \
	load-balancing \
	mpi-grid.tst mpi-periodic-3D.tst \
	source.tst tag.tst tag1.tst view.tst view.3D.tst \
	boundary_vertex.tst boundary_vertex3D.tst \
	foreach_bnd1.tst vertices-bc.tst

indexing.tst:		CC = mpicc -D_MPI=3
indexing.3D.tst:	CC = mpicc -D_MPI=3
mpi-restriction.tst:	CC = mpicc -D_MPI=3
mpi-restriction.3D.tst:	CC = mpicc -D_MPI=3
mpi-reduce.tst:		CC = mpicc -D_MPI=3

openmp-reduce.c: mpi-reduce.c
	ln -sf mpi-reduce.c openmp-reduce.c
openmp-reduce.s: CFLAGS += -fopenmp
openmp-reduce.tst: CFLAGS += -fopenmp

mpi-refine.tst:		CC = mpicc -D_MPI=4
mpi-refine1.tst:	CC = mpicc -D_MPI=11
mpi-refine.3D.tst:	CC = mpicc -D_MPI=4
mpi-laplacian.tst:	CC = mpicc -D_MPI=3
mpi-laplacian.3D.tst:	CC = mpicc -D_MPI=3
mpi-circle.tst:		CC = mpicc -D_MPI=5
foreach_bnd1.tst:	CC = mpicc -D_MPI=4
vertices-bc.tst:	CC = mpicc -D_MPI=4

mpi-circle1.c: mpi-circle.c
	ln -sf mpi-circle.c mpi-circle1.c

mpi-circle1.tst:  CC = mpicc -D_MPI=6
mpi-flux.tst:     CC = mpicc -D_MPI=6
mpi-interpu.tst:  CC = mpicc -D_MPI=5
mpi-coarsen.tst:  CC = mpicc -D_MPI=2
mpi-coarsen1.tst: CC = mpicc -D_MPI=5
bump2Dp.tst:      CC = mpicc -D_MPI=55
vortex.s:         CFLAGS = -DJACOBI=1
vortex.tst:	  CC = mpicc -D_MPI=7 -DJACOBI=1
axiadvection.s:   CFLAGS = -DJACOBI=1
axiadvection.tst: CC = mpicc -D_MPI=7 -DJACOBI=1
hf1.tst:	  CC = mpicc -D_MPI=7
poiseuille-periodic.tst: CC = mpicc -D_MPI=4
source.tst: CC = mpicc -D_MPI=7

tag.tst: CC = mpicc -D_MPI=7
tag1.tst: CC = mpicc -D_MPI=7

boundary_vertex.tst: CC = mpicc -D_MPI=7
boundary_vertex3D.tst: CC = mpicc -D_MPI=7

# parallel dump()/restore()

pdump.tst: CC = mpicc -D_MPI=7
restore.dump: pdump/restore.dump
	ln -sf pdump/restore.dump
restore.tst: pdump.tst restore.dump

pdump-multigrid.c: pdump.c
	ln -sf pdump.c pdump-multigrid.c
pdump-multigrid.tst: pdump.c
pdump-multigrid.tst: CFLAGS = -grid=multigrid
pdump-multigrid.tst: CC = mpicc -D_MPI=16

restore-multigrid.dump: pdump-multigrid/restore.dump
	ln -sf pdump-multigrid/restore.dump restore-multigrid.dump
restore-multigrid.tst: pdump-multigrid.tst restore-multigrid.dump
restore-multigrid.tst: CFLAGS = -grid=multigrid
restore-multigrid.tst: CC = mpicc -D_MPI=16

restore-tree.dump: pdump/restore.dump
	ln -sf pdump/restore.dump restore-tree.dump
restore-tree.tst: pdump.tst restore-tree.dump
restore-tree.tst: CFLAGS = -DDEBUGCOND=false
restore-tree.tst: CC = mpicc -D_MPI=23

restore-multigrid-serial.dump: pdump-multigrid/restore.dump
	ln -sf pdump-multigrid/restore.dump restore-multigrid-serial.dump
restore-multigrid-serial.tst: pdump-multigrid.tst restore-multigrid-serial.dump

restore-rectangle.dump: dump-rectangle/restore-rectangle.dump
	ln -sf dump-rectangle/restore-rectangle.dump
restore-rectangle.tst: dump-rectangle.tst restore-rectangle.dump

bump2Dp-restore.c: bump2Dp.c
	ln -sf bump2Dp.c bump2Dp-restore.c
bump2Dp-restore.dump: bump2Dp/dump
	ln -sf bump2Dp/dump bump2Dp-restore.dump
bump2Dp-restore.tst: bump2Dp.tst bump2Dp-restore.dump
bump2Dp-restore.tst: CC = mpicc -D_MPI=15

# parallel output_gfs()

gfs.tst: CC = mpicc -D_MPI=7
gfsi.gfs: gfs/gfsi.gfs
	ln -sf gfs/gfsi.gfs
gfsi.tst: gfs.tst gfsi.gfs

# load-balancing

load-balancing: balance5.tst balance6.tst balance7.tst \
		bump2Dp.tst bump2Dp-restore.tst vortex.tst axiadvection.tst

balance5.tst: CC = mpicc -D_MPI=9
balance6.c: balance5.c
	ln -sf balance5.c balance6.c
balance6.tst: CC = mpicc -D_MPI=17
balance7.tst: CC = mpicc -D_MPI=17

# MPI-parallel multigrid

mpi-grid.tst: CC = mpicc -D_MPI=4
mpi-periodic-3D.c: periodic.c
	ln -sf periodic.c mpi-periodic-3D.c
mpi-periodic-3D.tst: CC = mpicc -D_MPI=8
mpi-periodic-3D.tst: CFLAGS += -grid=multigrid3D

# CADNA tests
# These are ignored for the moment since nobody seems to be
# interested and this is complicated to maintain and expensive to run

cadna-tests: # poiseuille.CADNA.tst sag.CADNA.tst poisson.CADNA.tst

# Periodic boundary tests

periodic1.tst: CFLAGS += -grid=bitree -DTRASH=1
periodic2.c: periodic1.c
	ln -sf periodic1.c periodic2.c
periodic2.tst: CFLAGS += -DTRASH=1
periodic3.c: periodic1.c
	ln -sf periodic1.c periodic3.c
periodic3.tst: CFLAGS += -grid=octree -DTRASH=1

# Embedded boundary tests

dirichlet.tst: CFLAGS += -DDIRICHLET=1
dirichlet.c: neumann.c
	ln -sf neumann.c dirichlet.c

dirichlet-axi.c: neumann-axi.c
	ln -sf neumann-axi.c dirichlet-axi.c
dirichlet-axi.s: CFLAGS += -DDIRICHLET=1
dirichlet-axi.tst: CFLAGS += -DDIRICHLET=1

dirichlet3D.tst: CFLAGS += -DDIRICHLET=1
dirichlet3D.c: neumann3D.c
	ln -sf neumann3D.c dirichlet3D.c

hydrostatic2.tst: CC = mpicc -D_MPI=8

embed-tests: neumann.tst dirichlet.tst poiseuille45.tst \
	couette.tst wannier.tst hydrostatic2.tst cylinders.tst \
	porous.tst porous1.tst uf.tst starting.tst \
	neumann3D.tst dirichlet3D.tst hydrostatic3.tst spheres.tst \
        neumann-axi.tst dirichlet-axi.tst

# General Ocean Turbulence Model (GOTM) tests

gotm-tests: couette-gotm.tst entrainment.tst ows_papa.tst

ows_papa-kpp.c: ows_papa.c
	ln -sf ows_papa.c ows_papa-kpp.c 
ows_papa-kpp.s: CFLAGS += -DKPP=1
ows_papa-kpp.tst: CFLAGS += -DKPP=1

ows_papa.tst: ows_papa-kpp.tst

# Stencils

stencils.tst: CFLAGS += -DPRINTBOUNDARY=1

# CLSVOF

capwave-clsvof.c: capwave.c
	ln -sf capwave.c capwave-clsvof.c 
capwave-clsvof.s: CFLAGS += -DCLSVOF=1
capwave-clsvof.tst: CFLAGS += -DCLSVOF=1

# Dimensions

dimensions-tests: bump2D.s capwave.s lid.s rising-axi.s		\
	rising-reduced.s vortex.s sessile.s large.s lock.s	\
	wind-driven-stvt.s shock.s starting.s layered.s 	\
	stokes.s stokes-ns.s lonlat.s overflow.s bleck.s

# GPUs
# fixme: the dependencies *.gpu.s.d are not created automatically
gpu-tests: gpu.gpu.tst gpu-cpu.gpu.tst gpu-reduce.gpu.tst layers.gpu.tst nonlinear.gpu.tst \
	advection.gpu.tst bump2D-gpu.tst reversed.gpu.tst axi.gpu.tst lid.gpu.tst \
	stommel-ml.gpu.tst bleck.gpu.tst rising-clsvof.gpu.tst
gpu-tests: CFLAGS += -DPRINTNSHADERS

gpu-cpu.c: gpu.c
	ln -sf gpu.c gpu-cpu.c
gpu-cpu.gpu.tst: CFLAGS += -cpu

rising-clsvof.gpu.c: rising.c
	ln -sf rising.c rising-clsvof.gpu.c
rising-clsvof.gpu.s: CFLAGS += -DLEVELSET=1 -DCLSVOF=1
rising-clsvof.gpu.tst: CFLAGS += -DLEVELSET=1 -DCLSVOF=1

include ../Makefile.defs
