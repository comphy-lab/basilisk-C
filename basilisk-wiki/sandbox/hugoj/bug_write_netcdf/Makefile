CFLAGS += -I$(DOCUMENT_ROOT)/sandbox
CFLAGS += -disable-dimensions # skip dimensions check for now

all: main.tst  mainGPU.tst #main1D.tst

mainGPU.c: main.c
	ln -sf main.c mainGPU.c

mainGPU.tst: CFLAGS += -grid=gpu/multigrid

include $(BASILISK)/Makefile.defs

