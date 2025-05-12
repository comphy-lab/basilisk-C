LIBS = -lgsl
include $(BASILISK)/Makefile.defs

QCC = $(BASILISK)/qcc
CFLAGS += -DMTRACE=3 -fopenmp -O3
BVIEW=-L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa

PYINC=/usr/include/python3.8
PYLIB=/usr/lib/python3.8/config-3.8-x86_64-linux-gnu -L/usr/lib

read_from_python:
	$(QCC) $(CFLAGS) -o read_from_python/read_from_python read_from_python.c $(BVIEW) -I$(PYINC) -L$(PYLIB) -lpython3.8  -lm
