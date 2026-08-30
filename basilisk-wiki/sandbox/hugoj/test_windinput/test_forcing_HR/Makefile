CFLAGS += -I$(DOCUMENT_ROOT)/sandbox -O2
CFLAGS += -disable-dimensions # skip dimensions check for now

all: linear_wave_wind_input.tst linear_wave_wind_input/plots

linear_wave_wind_input.tst: CFLAGS+="-DFORCING=1"


no_forcing.c: linear_wave_wind_input.c
	ln -sf linear_wave_wind_input.c no_forcing.c

no_forcing.tst: CFLAGS+="-fopenmp"

include $(BASILISK)/Makefile.defs


