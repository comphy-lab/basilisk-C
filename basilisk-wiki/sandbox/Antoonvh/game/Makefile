CFLAGS += -O2

include $(BASILISK)/Makefile.defs

libcur.a: cur.c
	cc -c cur.c
	ar cr libcur.a cur.o

snake2.s: libcur.a

itmov.tst: CC=mpicc -D_MPI=3

testd2.tst: CC=mpicc -D_MPI=6

MPI_and_embed.tst: CC=mpicc -D_MPI=4

kizner2.tst: CFLAGS += -fopenmp

tube.tst: CFLAGS += -fopenmp

ash.tst: CFLAGS += -fopenmp -grid=octree

test_tag.tst: CC=mpicc -D_MPI=7
