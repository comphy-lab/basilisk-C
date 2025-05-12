CFLAGS = -Wall -O2
LIB = -lm
HEADERPATH=-I../ 

all: stl2raw kelvin_parallel  
stl2raw: stl2raw.c 
	qcc $(CFLAGS) -o $@ $< $(LIB)

kelvin_parallel: kelvin_parallel.c   
	CC99='mpicc -std=c99' qcc $S $(CFLAGS) $(HEADERPATH) -D_MPI=1 -o $@ $< $(LIB)

clean:
	rm -rf vti*
	rm -rf vtp*
	rm -rf htg*
	rm *.html*
	rm -f stl2raw \
		_kelvin_parallel.c \
		kelvin_parallel \
		kelvin*.raw \
		perfs \

	m -f *.pvd
