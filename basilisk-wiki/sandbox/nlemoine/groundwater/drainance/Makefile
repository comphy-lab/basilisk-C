CFLAGS += -O2
include $(BASILISK)/Makefile.defs

test_fortran.s: libexpokit.a

libexpokit.a: expokit.o
	ar rcs libexpokit.a expokit.o

expokit.o: expokit.f
	gfortran -cpp $(CFLAGS) -o $@ -c $<