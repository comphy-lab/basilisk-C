include $(BASILISK)/Makefile.defs
CFLAGS+= -O0 -g -disable-dimensions
HDFLIBS=-I$(HDF5_INCDIR) -L$(HDF5_LIBDIR) -lhdf5 -lhdf5_hl
OPENGLIBS=-L$(BASILISK)/gl -lglutils -lfb_tiny 
LIBS+=-I$(BASILISK) $(OPENGLIBS) $(HDFLIBS) -lgsl -lgslcblas