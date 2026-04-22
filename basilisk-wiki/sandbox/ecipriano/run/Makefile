include $(BASILISK)/Makefile.defs

CFLAGS += -O2 -disable-dimensions -DBASILISK_SERVER=1 -DUSE_GSL=1 -DUSE_CANTERA=1 -I../src

expansion-T375.c: expansion.c
	ln -s expansion.c expansion-T375.c
expansion-T375.s: CFLAGS += -DTEMPERATURE=375
expansion-T375.tst: CFLAGS += -DTEMPERATURE=375

expansion-T400.c: expansion.c
	ln -s expansion.c expansion-T400.c
expansion-T400.s: CFLAGS += -DTEMPERATURE=400
expansion-T400.tst: CFLAGS += -DTEMPERATURE=400

expansion.tst: expansion-T375.tst expansion-T400.tst

microgravity.yaml:
	ln -s normalgravity.yaml microgravity.yaml

microgravity.tst: microgravity.yaml

microgravity.tst: CFLAGS += -DJUMP=0 -DTEMPERATURE=565 -DTEMPERATURE_DROPLET=363 \
	-DPRESSURE=28.226 -DSIGMA=0.01 -DKINFOLDER=microgravity.yaml -DMAXLEVEL=6 \
	-DTTOL=1e-2 -DUTOL=1e-2 -DDIAMETER=5e-6 -DENVIRONMENT=4 -DMAX_DD02=0.05

normalgravity.tst: CFLAGS += -DJUMP=0 -DMOVIE=1 -DMOVIE_EVERY=0.002 -DTEMPERATURE=773 \
	-DPRESSURE=10 -DKINFOLDER=normalgravity.yaml -DGRAVITY=-9.81

burningdroplet.tst: CFLAGS += -DKINFOLDER=burningdroplet.yaml -DMOVIE_EVERY=0.0001

include Makefile.dynamic
