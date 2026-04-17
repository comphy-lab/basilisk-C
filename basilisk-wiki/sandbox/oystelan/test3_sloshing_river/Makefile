# makefile for sloshing river testcase
PPR_LIBS=$(BASILISK)/ppr

.PHONY: default clean

default: river ;

river: river.c
	CC99='gcc -std=c99' qcc -fopenmp river.c -o river -lm -L. -ldl -lstdc++ -Wall -O2 -L$(PPR_LIBS) -lppr -lgfortran 

clean:
	rm *.o ; \
	rm river; 



