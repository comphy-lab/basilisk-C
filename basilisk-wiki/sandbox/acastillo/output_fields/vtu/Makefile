include $(BASILISK)/Makefile.defs
CFLAGS+= -O3

HDFLIBS=-I$(HDF5_INCDIR) -L$(HDF5_LIBDIR) -lhdf5 -lhdf5_hl
OPENGLIBS=-L$(BASILISK)/gl -lglutils -lfb_tiny -lgsl -lgslcblas
LIBS+=-I$(BASILISK) $(OPENGLIBS) 