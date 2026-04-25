CFLAGS += -I$(DOCUMENT_ROOT)/sandbox #-O2
CFLAGS += -disable-dimensions # skip dimensions check for now

EXEC=ml_breaking.gpu.tst

# > local mpi run
${EXEC}.tst: CC = mpicc -D_MPI=16

# > gpu (remap_C)
# dont forget the env variables to select nvidia gpu !
# __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia make ml_breaking.gpu.tst
${EXEC}.gpu.c: ${EXEC}.c
	ln -sf ${EXEC}.c ${EXEC}.gpu.c

# > run on HPC
_${EXEC}.c: CFLAGS += -D_MPI=1



include $(BASILISK)/Makefile.defs
