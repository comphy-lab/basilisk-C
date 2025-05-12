CFLAGS += -O2
SOURCE := isotropic.c

all:
	make wclean
	make 2d mpi2d

2d: isotropic.c
	rm -f _isotropic.c
	make clean
	$(QCC) $(CFLAGS) -D_MPI=0 -o $@ $^ -lm

mpi2d: isotropic.c
	rm -f _isotropic.c
	make clean
	$(QCC) -source $(CFLAGS) -D_MPI=1 $<
	mpicc $(CFLAGS) -o $@ _isotropic.c -lm

wclean:
	rm -f *.exec _*.c 2d mpi2d rm -rf .qcc*
	make clean
	
include $(BASILISK)/Makefile.defs