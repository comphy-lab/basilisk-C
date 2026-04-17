#makefile for instability.c testcase
PPR_LIBS=$(BASILISK)/ppr


.PHONY: default clean

default: instability ;
instability: instability.c
	CC99='gcc -std=c99' qcc instability.c -o instability -lm -L. -ldl -lstdc++ -Wall -O2 -L$(PPR_LIBS) -lppr -lgfortran 

clean:
	rm *.o ; \
	rm instability; 



