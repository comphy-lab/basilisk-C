CFLAGS += -I$(DOCUMENT_ROOT)/sandbox
CFLAGS += -O2 -disable-dimensions # skip dimensions check for now

# > local mpi run
# ml_breaking.tst: CC = mpicc -D_MPI=16

# > gpu (remap_C)
# dont forget the env variables to select nvidia gpu !
# __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia make ml_breaking.gpu.tst

# > run on HPC
_ml_breaking.c: CFLAGS += -D_MPI=1

include $(BASILISK)/Makefile.defs
