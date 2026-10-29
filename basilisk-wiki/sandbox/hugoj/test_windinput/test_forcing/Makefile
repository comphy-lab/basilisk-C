CFLAGS += -I$(DOCUMENT_ROOT)/sandbox -O2
CFLAGS += -disable-dimensions # skip dimensions check for now
CFLAGS += -DOUTNC

all: no_forcing.tst linear_wave_wind_input.tst exact_forcing.tst linear_wave_wind_input/plots

# 1D
linear_wave_wind_input.tst: CFLAGS+= -grid=multigrid1D -DHOLD_FORCING=1

no_forcing.c: linear_wave_wind_input.c
	ln -sf linear_wave_wind_input.c no_forcing.c
no_forcing.tst: CFLAGS += -grid=multigrid1D

exact_forcing.c: linear_wave_wind_input.c
	ln -sf linear_wave_wind_input.c exact_forcing.c
exact_forcing.tst: CFLAGS += -DEXACT_FORCING=1 -grid=multigrid1D

stokes.c: linear_wave_wind_input.c
	ln -sf linear_wave_wind_input.c stokes.c
stokes.tst: CFLAGS += -DSTOKES -grid=multigrid1D

hold_stokes.c: linear_wave_wind_input.c
	ln -sf linear_wave_wind_input.c hold_stokes.c
hold_stokes.tst: CFLAGS += -DSTOKES -DHOLD_FORCING=1 -grid=multigrid1D

# 2D
no_forcing2D.c: linear_wave_wind_input.c
	ln -sf linear_wave_wind_input.c no_forcing2D.c

linear_wave_wind_input2D.c: linear_wave_wind_input.c
	ln -sf linear_wave_wind_input.c linear_wave_wind_input2D.c
linear_wave_wind_input2D.tst: CFLAGS+= -DHOLD_FORCING=1
linear_wave_wind_input2D.cuda.tst: CFLAGS+= -DHOLD_FORCING=1
linear_wave_wind_input2D.gpu.tst: CFLAGS+= -DHOLD_FORCING=1

exact_forcing2D.c: linear_wave_wind_input.c
	ln -sf linear_wave_wind_input.c exact_forcing2D.c
exact_forcing2D.tst: CFLAGS += -DEXACT_FORCING=1 
exact_forcing2D.cuda.tst: CFLAGS += -DEXACT_FORCING=1 

stokes2D.c: linear_wave_wind_input.c
	ln -sf linear_wave_wind_input.c stokes2D.c
stokes2D.tst: CFLAGS += -DSTOKES 

hold_stokes2D.c: linear_wave_wind_input.c
	ln -sf linear_wave_wind_input.c hold_stokes2D.c
hold_stokes2D.tst: CFLAGS += -DSTOKES -DHOLD_FORCING=1 







include $(BASILISK)/Makefile.defs


