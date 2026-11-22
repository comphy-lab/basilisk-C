.PHONY: all clean

BASE_CFLAGS = -I$(BASILISK) -I$(BASILISK)/grid -I.. -O2 -g -Wall

all: libbuda.a libbuda64.a

libbuda.a: cuda.o
	ar cr libbuda.a cuda.o

libbuda64.a: cuda64.o
	ar cr libbuda64.a cuda64.o

cuda.o: cuda.c
	$(CC) $(BASE_CFLAGS) -DSINGLE_PRECISION=1 -c -o $@ $<

cuda64.o: cuda.c
	$(CC) $(BASE_CFLAGS) -DSINGLE_PRECISION=0 -DDOUBLE_PRECISION=1 -c -o $@ $<

clean:
	rm -f cuda.o cuda64.o libbuda.a libbuda64.a
