#CC=gcc # the C compiler
CC=mpicc -std=c99 # the C compiler
CFLAGS += -O2 -lm # optimization level, link math
LIB_DIR=$(CURDIR) # current working directory. CURDIR is a built-in macro in make
COMPILER_LIBS=/usr/lib/gcc/x86_64-linux-gnu/9
INCLUDE_DIR=/usr/lib/openmpi/include # directory to check for header files
export BASILISK := /home/oystein/progs/basilisk/src
export PATH:=/home/oystein/progs/basilisk/src:${PATH}

KDT_LIBS=$(BASILISK)/kdt
PPR_LIBS=$(BASILISK)/ppr



# for omp: export OMP_NUM_THREADS=4

.PHONY: default clean

default: instability ;

instability: instability.c
	CC99='gcc -std=c99' qcc instability.c -o instability -lm -L. -ldl -lstdc++ -Wall -O2 -L$(KDT_LIBS) -L$(PPR_LIBS) -lkdt -lppr -lgfortran 

clean:
	rm *.o ; \
	rm instability; 



