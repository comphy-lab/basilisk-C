CFLAGS += -I$(DOCUMENT_ROOT)/sandbox -O2
CFLAGS += -disable-dimensions # skip dimensions check for now

all: diff_test_H.tst diff_test.tst diff_test/plots

diff_test_H.c: diff_test.c
	ln -sf diff_test.c diff_test_H.c
diff_test_H.tst: CFLAGS+='-DHEATING'

# prep: hydro.h implicit.h diffusion.h
# 	ln -sf ${BASILISK}/layered/hydro.h .
# 	ln -sf ${BASILISK}/layered/implicit.h .
# 	ln -sf ${DOCUMENT_ROOT}/sandbox/hugoj/lib/diffusionRheo.h diffusion.h


include $(BASILISK)/Makefile.defs


