CFLAGS += -O2
include $(BASILISK)/Makefile.defs
envelope.s: libaem.a

SRCF77 = $(wildcard *.f)
OBJF77 = $(SRCF77:.f=.f.o)

libaem.a: $(OBJF77)
	ar rcs $@ $(OBJF77)
	rm -f $(OBJF77)

%.f.o: %.f
	gfortran -o $@ -c $<