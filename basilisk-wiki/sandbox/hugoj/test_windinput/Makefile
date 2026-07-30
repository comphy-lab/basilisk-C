CFLAGS += -I$(DOCUMENT_ROOT)/sandbox -O2
CFLAGS += -disable-dimensions # skip dimensions check for now

all: linear_wave_wind_input.tst linear_wave_wind_input/plots

linear_wave_wind_input.tst: CFLAGS+="-DFORCING=1"


no_forcing.c: linear_wave_wind_input.c
	ln -sf linear_wave_wind_input.c no_forcing.c

stokes.c: linear_wave_wind_input.c
	ln -sf linear_wave_wind_input.c stokes.c
stokes.tst: CFLAGS+="-DSTOKES=1 "

no_forcing_instrumented.tst: CFLAGS+="-DDIAG=1 -fopenmp"

no_forcing_instrumented.c: prep_instrumented


prep_instrumented:
	(	cd modified_layered;\
	ln -sf ${BASILISK}/layered/implicit.h .; \
	ln -sf ${BASILISK}/layered/diffusion.h .;\
	ln -sf ${BASILISK}/layered/remap.h .;)


include $(BASILISK)/Makefile.defs


