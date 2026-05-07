# makefile for sloshing river testcase

.PHONY: default clean

default: river ;

river: river.c
	CC99='gcc -std=c99' qcc -fopenmp river.c -o river -lm -L. -ldl -lstdc++ -Wall -O2

clean:
	rm *.o ; \
	rm river; 



