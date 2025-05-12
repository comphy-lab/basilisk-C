include $(BASILISK)/Makefile.defs

#BASILISK_INCLUDES ?= -I../src
#BASILISK_LIBS ?= 

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -disable-dimensions -DBASILISK_SERVER=1 -I../src

fixedbubblevelocity-jump.c: fixedbubblevelocity.c
	ln -s fixedbubblevelocity.c fixedbubblevelocity-jump.c
fixedbubblevelocity-jump.s: CFLAGS += -DJUMP=1
fixedbubblevelocity-jump.tst: CFLAGS += -DJUMP=1

fixedbubblevelocity.tst: fixedbubblevelocity-jump.tst


expansion-T375.c: expansion.c
	ln -s expansion.c expansion-T375.c
expansion-T375.s: CFLAGS += -DTEMPERATURE=375
expansion-T375.tst: CFLAGS += -DTEMPERATURE=375

expansion-T400.c: expansion.c
	ln -s expansion.c expansion-T400.c
expansion-T400.s: CFLAGS += -DTEMPERATURE=400
expansion-T400.tst: CFLAGS += -DTEMPERATURE=400

expansion.tst: expansion-T375.tst expansion-T400.tst


normalgravity.tst: CFLAGS += -DJUMP=1 -DMOVIE=1 -DDIAMETER=1e-3 -DPRESSURE=10 -DRADIATION_INTERFACE=0
