CFLAGS += -I$(DOCUMENT_ROOT)/sandbox -O2
CFLAGS += -disable-dimensions # skip dimensions check for now
CFLAGS += -DDIAG=0
CFLAGS += -DTRACE=2

all: ml_convection/namelist.toml ml_convection.tst


ml_convection/namelist.toml: namelist.toml
	mkdir -p ml_convection
	ln -sf ../namelist.toml ml_convection/namelist.toml

ml_convection.gpu/namelist.toml: namelist.toml
	mkdir -p ml_convection.gpu
	ln -sf ../namelist.toml ml_convection.gpu/namelist.toml


ml_convection.tst: CC = mpicc -D_MPI=16

_ml_convection.c: CFLAGS += -D_MPI=1

include $(BASILISK)/Makefile.defs


