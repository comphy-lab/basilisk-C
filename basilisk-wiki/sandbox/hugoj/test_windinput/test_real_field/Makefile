CFLAGS += -I$(DOCUMENT_ROOT)/sandbox -O2
CFLAGS += -disable-dimensions # skip dimensions check for now

PM_seastate.tst: CFLAGS+="-DFORCING=1 -DDUMPS=1"

PM_seastate.cuda.tst: CFLAGS+="-DDUMPS=1"

include $(BASILISK)/Makefile.defs


