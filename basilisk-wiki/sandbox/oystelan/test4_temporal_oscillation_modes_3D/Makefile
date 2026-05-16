# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -fopenmp -Wdimensions -DMTRACE=3 -I$(DOCUMENT_ROOT)/sandbox

include $(BASILISK)/Makefile.defs

instability3D-new.c: instability3D.c
	ln -sf instability3D.c instability3D-new.c
instability3D-new.s: CFLAGS += -Dimplicitness=0.505
instability3D-new.tst: CFLAGS += -Dimplicitness=0.505

instability3D.tst: instability3D-new.tst
