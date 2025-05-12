CFLAGS += -O2 -Wdimensions

include $(BASILISK)/Makefile.defs

test_prestore.tst: CC=mpicc -D_MPI=3

tpp.tst: CC=mpicc -D_MPI=6

itmov.tst: CC=mpicc -D_MPI=3

tp2.tst: CC=mpicc -D_MPI=5

test_manyparticles.tst: CC=mpicc -D_MPI=6

testd2.tst: CC=mpicc -D_MPI=6

MPI_and_embed.tst: CC=mpicc -D_MPI=4

kizner2.tst: CFLAGS += -fopenmp

tube.tst: CFLAGS += -fopenmp

test_tag.tst: CC=mpicc -D_MPI=7

cor.tst: CFLAGS += -fopenmp

ring.tst: CFLAGS += -fopenmp

compactring.tst: CFLAGS += -fopenmp

helical.tst: CFLAGS += -fopenmp

impinging.tst: CFLAGS += -fopenmp

test_Hilbert.tst: CC=mpicc -D_MPI=3

visl2.tst: CFLAGS += -fopenmp

moist_bubble.tst: CFLAGS += -fopenmp

inject_particles.tst: CC=mpicc -D_MPI=4

dcomp.tst: CC=mpicc -D_MPI=5

test6-mpi.tst: CC=mpicc -D_MPI=3

vary_threads: CFLAGS += fopenmp

omega.tst: slave-omgpsi.o

slave-omgpsi.o: slave-omgpsi.c Makefile
	qcc -O2 -fno-common -c slave-omgpsi.c
	objcopy -G slave_free -G slave_interpolate -G slave_solve_psi -G slave_level -G slave_init slave-omgpsi.o
	cp -f slave-omgpsi.o omega/

hilbert.tst: CFLAGS += -disable-dimensions

multiscale-gpu.tst: CFLAGS += -grid=gpu/multigrid