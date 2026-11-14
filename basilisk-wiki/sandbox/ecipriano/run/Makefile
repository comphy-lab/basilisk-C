include $(BASILISK)/Makefile.defs

CFLAGS += -O2 \
	-DMOVIE=1 \
	-DUSE_GSL=1 \
	-DUSE_CANTERA=1 \
	-DUSE_OPENSMOKE=0 \
	-DTEST_RESTART=0 \
	-DMTRACE=0 \
	-DTRACE=0 \
	-DBASILISK_SERVER=1 \
	-disable-dimensions \
	-I../src \
	$(BASILISK_INCLUDES) \
	$(BASILISK_LIBS)

ifdef CANTERA_LIBRARY_PATH
CFLAGS += -I$(CANTERA_INCLUDE_PATH) -L$(CANTERA_LIBRARY_PATH) -Wl,-rpath,$(CANTERA_LIBRARY_PATH)
endif

microgravity.yaml:
	ln -s normalgravity.yaml microgravity.yaml

microgravity.tst: microgravity.yaml

microgravity.tst: CFLAGS += -DJUMP=0 -DTEMPERATURE=565 -DTEMPERATURE_DROPLET=363 \
	-DPRESSURE=28.226 -DSIGMA=0.01 -DKINFOLDER=microgravity.yaml -DMAXLEVEL=6 \
	-DTTOL=1e-2 -DUTOL=1e-2 -DDIAMETER=5e-6 -DENVIRONMENT=4 -DMAX_DD02=0.05

normalgravity.tst: CFLAGS += -DJUMP=0 -DMOVIE_EVERY=0.002 -DTEMPERATURE=773 \
	-DPRESSURE=10 -DKINFOLDER=normalgravity.yaml -DGRAVITY=-9.81

burningdroplet.tst: CFLAGS += -DKINFOLDER=burningdroplet.yaml -DMOVIE_EVERY=0.0001

# Low-mach with pressurization
pressurization-embed.c: pressurization.c
	ln -s pressurization.c pressurization-embed.c
pressurization-embed.s: CFLAGS += -DSOLID=1
pressurization-embed.tst: CFLAGS += -DSOLID=1

pressurization-heatflux.c: pressurization.c
	ln -s pressurization.c pressurization-heatflux.c
pressurization-heatflux.s: CFLAGS += -DHEATFLUX=1
pressurization-heatflux.tst: CFLAGS += -DHEATFLUX=1

pressurization-heatflux-embed.c: pressurization.c
	ln -s pressurization.c pressurization-heatflux-embed.c
pressurization-heatflux-embed.s: CFLAGS += -DHEATFLUX=1 -DSOLID=1
pressurization-heatflux-embed.tst: CFLAGS += -DHEATFLUX=1 -DSOLID=1

pressurization.tst: pressurization-embed.tst pressurization-heatflux.tst \
	pressurization-heatflux-embed.tst piston.tst cavity.tst

