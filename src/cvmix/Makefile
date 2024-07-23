# where is CVMIX?
ifeq ($(CVMIX),)
	CVMIX = $(HOME)/local/src/CVMix-src
endif

# how to compile?
ifeq ($(F90),)
	F90 := gfortran
endif
ifeq ($(FCFLAGS),)
	FCFLAGS := -Wall -O2
endif

# This should be completed as new objects are added to CVMix
HEADERS = kinds_and_types.h 	\
          background.h 		\
          convection.h 		\
          ddiff.h 		\
          kpp.h 		\
          math.h 		\
          put_get.h 		\
          shear.h 		\
          tidal.h 		\
          utils.h

OBJECTS = $(subst .h,.o,$(HEADERS))

libcvmixc.a: $(CVMIX)/lib/libcvmix.a $(HEADERS) common.o $(OBJECTS)
	cp -f $(CVMIX)/lib/libcvmix.a libcvmixc.a
	ar -rs libcvmixc.a common.o $(OBJECTS)
	objcopy --globalize-symbol=__cvmix_kpp_MOD_compute_phi_inv libcvmixc.a

%.h: $(CVMIX)/src/shared/cvmix_%.F90 comments.awk cvmix.awk
	( tr '[:upper:]' '[:lower:]' < $< | awk -f comments.awk \
	| awk -f cvmix.awk -vbase=$* -v sourcefile="$<" > $@ ) || \
	( rm -f $*.F90 && exit 1 )

%.o: %.F90
	$(F90) $(FCFLAGS) -I$(CVMIX)/include -c $<

clean:
	rm -f *.o $(HEADERS)
