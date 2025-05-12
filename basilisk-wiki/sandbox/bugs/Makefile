CFLAGS += -O2 -Wdimensions

include $(BASILISK)/Makefile.defs

wavelet_bug.tst: CC=mpicc -D_MPI=6

automatic.tst: CC=mpicc -D_MPI=4

restore.c: restart-restore.c
	ln -s restart-restore.c restore.c
restart-restore.dump: restore.tst
	rm -f restart-restore.dump; ln -s restore/dump restart-restore.dump
restart-restore.tst: restart-restore.dump

foreach_neighbor_test.tst: CC=mpicc -D_MPI=8

stack.tst: CFLAGS += -D_MPI=1 -grid=octree

reductions.tst: CC=mpicc -D_MPI=4
