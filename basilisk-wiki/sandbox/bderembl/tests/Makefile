# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -I$(DOCUMENT_ROOT)/sandbox

include $(BASILISK)/Makefile.defs

vertices.tst: CC = mpicc -D_MPI=4

fraction_memory.tst: CFLAGS += -grid=multigrid3D -DMTRACE=1

wavelet_vertex.tst: CFLAGS += -grid=multigrid1D

wavelet_vertex_2d.tst: CFLAGS += -grid=multigrid

boussinesq.tst: CFLAGS += -grid=multigrid
