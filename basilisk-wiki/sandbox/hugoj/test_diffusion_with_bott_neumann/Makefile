CFLAGS += -I$(DOCUMENT_ROOT)/sandbox -O2
CFLAGS += -disable-dimensions # skip dimensions check for now

all: prep diff_test.tst diff_test/plots


prep: hydro.h implicit.h diffusion.h
	ln -sf ${BASILISK}/layered/hydro.h .
	ln -sf ${BASILISK}/layered/implicit.h .
	ln -sf ${DOCUMENT_ROOT}/sandbox/hugoj/lib/diffusionRheo.h diffusion.h


include $(BASILISK)/Makefile.defs


