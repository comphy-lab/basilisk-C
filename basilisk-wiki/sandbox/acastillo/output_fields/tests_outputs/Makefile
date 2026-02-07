include $(BASILISK)/Makefile.defs
CFLAGS+= -O2 -march=native -mtune=native -disable-dimensions

VTU_TESTS = test_output_vtu.tst \
            test_output_vtu_box.tst \
            test_output_vtu_facets.tst

XMF_TESTS = test_output_xmf.tst \
            test_output_xmf2.tst

HDF_TESTS = test_output_vtkhdf.tst \
            test_output_vtkhdf_box.tst

# Grid-specific flags
vtu-tree3D:  GFLAGS = -grid=octree
vtu-mgrid:   GFLAGS = -grid=multigrid
vtu-mgrid3D: GFLAGS = -grid=multigrid3D

xmf-tree3D:  GFLAGS = -grid=octree
xmf-mgrid:   GFLAGS = -grid=multigrid
xmf-mgrid3D: GFLAGS = -grid=multigrid3D

hdf-tree3D:  GFLAGS = -grid=octree
hdf-mgrid:   GFLAGS = -grid=multigrid
hdf-mgrid3D: GFLAGS = -grid=multigrid3D

# Notes: 
# - box in hdf-tree-mpi is empty, 3D is fine 
# - box in hdf-mgrid-mpi is empty, 3D is fine 
# - hdf-mgrid-mpi and hdf-mgrid3D-mpi are fine
# - hdf-mgrid-mpi and hdf-tree3D-mpi are fine

# Runner for all grid types
vtu-tree vtu-tree3D vtu-mgrid vtu-mgrid3D:
	CFLAGS='$(GFLAGS)' $(MAKE) $(VTU_TESTS)

xmf-tree xmf-tree3D xmf-mgrid xmf-mgrid3D:
	CFLAGS='$(GFLAGS)' $(MAKE) $(XMF_TESTS)

hdf-tree hdf-tree3D hdf-mgrid hdf-mgrid3D:
	CFLAGS='$(GFLAGS)' $(MAKE) $(HDF_TESTS)

# MPI wrapper for any of the above
%-mpi:
	CC='mpicc -D_MPI=8' $(MAKE) $*
