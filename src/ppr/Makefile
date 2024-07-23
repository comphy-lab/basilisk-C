CFLAGS = -O3

libppr.a: ppr_1d.o ppr.o
	ar cr libppr.a $^

%.o: %.f90 *.f90
	gfortran -cpp $(CFLAGS) -c $<

clean:
	rm -f *.o *.a
