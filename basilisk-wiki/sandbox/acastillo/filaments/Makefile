include $(BASILISK)/Makefile.defs
CFLAGS+= -O3 -disable-dimensions #-mavx2 -axCORE-AVX2,CORE-AVX512

HDFLIBS=-I$(HDF5_INCDIR) -L$(HDF5_LIBDIR) -lhdf5 -lhdf5_hl
#OPENGLIBS=-L$(BASILISK)/gl -lglutils -lfb_tiny #-L$(LOCAL_LIBS)/lib -lgsl -lgslcblas
OPENGLIBS=-L$(BASILISK)/gl -lglutils -lfb_tiny -lgsl -lgslcblas
#OPENGLIBS=-L$(BASILISK)/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lgsl -lgslcblas
LIBS+=-I$(BASILISK) $(OPENGLIBS) #$(HDFLIBS)

test_filament0:
	mkdir $@/
	$(BASILISK)/qcc -source -disable-dimensions -DMAXLEVEL=$(MAX) -D_MPI=1 $@.c
	mpicc -D_MPI=1 $(CFLAGS) _$@.c $(LIBS) -o $@/$@ -lm
	make clean

test_filament1:
	mkdir $@/
	$(BASILISK)/qcc -source -disable-dimensions -DMAXLEVEL=$(MAX) -D_MPI=1 $@.c
	mpicc -D_MPI=1 $(CFLAGS) _$@.c $(LIBS) -o $@/$@ -lm
	make clean

test_vortex0:
	mkdir $@/
	$(BASILISK)/qcc -source -disable-dimensions -DMAXLEVEL=$(MAX) -D_MPI=1 $@.c
	mpicc -D_MPI=1 $(CFLAGS) _$@.c $(LIBS) -o $@/$@ -lm
	make clean
