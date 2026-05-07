CFLAGS += -O2 -I$(DOCUMENT_ROOT)/sandbox 
CFLAGS += -disable-dimensions # skip dimensions check for now

all: breaking.tst

breaking.tst: CC = mpicc -D_MPI=1


include $(BASILISK)/Makefile.defs
