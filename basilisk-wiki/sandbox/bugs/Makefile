CFLAGS += -O2 -Wdimensions -DMTRACE=3 -I$(DOCUMENT_ROOT)/sandbox

include $(BASILISK)/Makefile.defs

wavelet_bug.tst: CC=mpicc -D_MPI=6

automatic.tst: CC=mpicc -D_MPI=4

restore.c: restart-restore.c
	ln -sf restart-restore.c restore.c
restart-restore.dump: restore.tst
	ln -sf restore/dump restart-restore.dump
restart-restore.tst: restart-restore.dump

foreach_neighbor_test.tst: CC=mpicc -D_MPI=8

stack.tst: CFLAGS += -D_MPI=1 -grid=octree

reductions.tst: CC=mpicc -D_MPI=4

multired.tst: CC=mpicc -D_MPI=4

bc_face_vector.tst: CC=mpicc -D_MPI=4 -DPRINTBOUNDARY=1

bc_vector.tst: CFLAGS += -DPRINTBOUNDARY=1

my_adapt_bug.tst: CC=mpicc -D_MPI=4
