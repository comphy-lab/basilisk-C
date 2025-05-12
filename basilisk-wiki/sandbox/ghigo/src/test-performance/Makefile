all: 	starting.tst \
	starting-bmr.tst \
	starting-uniform.tst \
	starting-moving.tst \
	starting-moving-bmr.tst \
	starting-moving-uniform.tst \
	sphere-unbounded.tst	 \
	sphere-unbounded-2.tst	 \
	sphere-unbounded-4.tst	 \
	sphere-unbounded-8.tst 	\
	sphere-unbounded-16.tst \
	sphere-unbounded-bview.tst

include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=2 -DTRACE=2 -DDISPLAY=1

# The following tests have special dependencies/compilation requirements

starting-bmr.tst: CFLAGS += -DBMR=1
starting-bmr.c: starting.c
	ln -s starting.c starting-bmr.c

starting-uniform.tst: CFLAGS += -DUNIFORM=1
starting-uniform.c: starting.c
	ln -s starting.c starting-uniform.c

starting-moving-bmr.tst: CFLAGS += -DBMR=1
starting-moving-bmr.c: starting-moving.c
	ln -s starting-moving.c starting-moving-bmr.c

starting-moving-uniform.tst: CFLAGS += -DUNIFORM=1
starting-moving-uniform.c: starting-moving.c
	ln -s starting-moving.c starting-moving-uniform.c

sphere-unbounded-2.tst: CC = mpicc -D_MPI=2
sphere-unbounded-2.c: sphere-unbounded.c
	ln -s sphere-unbounded.c sphere-unbounded-2.c

sphere-unbounded-4.tst: CC = mpicc -D_MPI=4
sphere-unbounded-4.c: sphere-unbounded.c
	ln -s sphere-unbounded.c sphere-unbounded-4.c

sphere-unbounded-8.tst: CC = mpicc -D_MPI=8
sphere-unbounded-8.c: sphere-unbounded.c
	ln -s sphere-unbounded.c sphere-unbounded-8.c

sphere-unbounded-16.tst: CC = mpicc -D_MPI=16
sphere-unbounded-16.c: sphere-unbounded.c
	ln -s sphere-unbounded.c sphere-unbounded-16.c

sphere-unbounded-bview.tst: CC = mpicc -D_MPI=16
sphere-unbounded-bview.tst: CFLAGS += -DBVIEW=1
sphere-unbounded-bview.c: sphere-unbounded.c
	ln -s sphere-unbounded.c sphere-unbounded-bview.c

