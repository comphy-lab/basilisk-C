
include $(BASILISK)/Makefile.defs
QCC = $(BASILISK)/qcc

CC99 = gcc
CFLAGS += -Wall -std=c99 -O3
OPENGLIBS=-lglutils -lfb_osmesa -lGLU -lOSMesa -lgsl -lgslcblas

two_helices2_mpi:
		export OMP_NUM_THREADS=1;
		export KMP_NUM_THREADS=1;
		qcc -D_MPI=16 $(QCFLAGS) -source two_helices2.c ;
		mpiicc $(CFLAGS) _two_helices2.c -I$(BASILISK) -L$(BASILISK)/gl $(OPENGLIBS) -o two_helices2/two_helices2 -lm ;
		cd two_helices2 ; mpirun -np 16 -genv I_MPI_PIN_PROCESSOR_LIST=allcores ./two_helices2 > out 2> log ; cd ..
		#cd two_helices2 ; mpirun -np 16 -genv I_MPI_PIN_PROCESSOR_LIST=allcores xterm -e gdb ./two_helices2 -ex run > out 2> log ; cd ..
