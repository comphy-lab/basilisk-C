include $(BASILISK)/Makefile.defs
OPTFLAGS=-fno-strict-aliasing -O3 -march=core-avx2 -axCORE-AVX2,CORE-AVX512 
CJSONLIBS=-I$(CJSON_INCDIR) -L$(CJSON_LIBDIR) -lcjson
CFLAGS+= $(OPTFLAGS) $(CJSONLIBS) -disable-dimensions
HDFLIBS=-I$(HDF5_INCDIR) -L$(HDF5_LIBDIR) -lhdf5 -lhdf5_hl
OPENGLIBS=-L$(BASILISK)/gl -lglutils -lfb_tiny 
LIBS+=-I$(BASILISK) $(OPENGLIBS) $(HDFLIBS) -lgsl -lgslcblas 
