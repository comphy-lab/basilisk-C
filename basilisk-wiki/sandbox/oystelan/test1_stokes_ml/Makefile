#makefile for severely breaking stokes wave in multilayer test


PPR_LIBS=$(BASILISK)/ppr

.PHONY: default clean

default: stokes;

stokes: stokes.c
	CC99='gcc -std=c99' qcc -fopenmp stokes.c -o stokes -O2 -fcommon -Wall -L. -L$(BASILISK) -lm -L$(PPR_LIBS) -lppr -lgfortran

clean:
	rm *.o ; \
	rm stokes; 



