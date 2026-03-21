
CFLAGS += -O2 -I$(DOCUMENT_ROOT)/sandbox 
CFLAGS += -disable-dimensions # skip dimensions check for now

# Jiarong's way of compiling
# field_PM/field_PM: field_PM.c
# 	mkdir -p field_PM
# 	qcc $(CFLAGS) -autolink -DPARAVIEW=0 -O2 -g -Wall -o field_PM/field_PM field_PM.c -I$(BASILISK) -L${BASILISK}/gl -lfb_tiny -lgfortran -L${BASILISK}/ppr -lppr -lm
#
# field_PMmpi/field_PM: field_PM.c
# 	mkdir -p field_PMmpi
# 	CC99='mpicc -std=c99' qcc $(CFLAGS)-autolink -DPARAVIEW=0 -O2 -g -Wall -o field_PM/field_PM field_PM.c -I$(BASILISK) -L${BASILISK}/gl -lfb_tiny -lgfortran -L${BASILISK}/ppr -lppr -lm



all: field_PM/F_kxky field_PM.tst 

# initial fields
field_PM/F_kxky: F_kxky kx ky
	mkdir -p field_PM
	mkdir -p field_PM/surface
	mkdir -p field_PM/field
	cp F_kxky kx ky ./field_PM
	
# running in // using mpi
field_PM.tst: CC = mpicc -D_MPI=16

# generate source code with
# make _field_PM.c

_field_PM.c: CFLAGS += -D_MPI=1


include $(BASILISK)/Makefile.defs


