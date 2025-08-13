all: 	network-flow.tst \
	network-pulsed.tst \
	wind-turbine-flow.tst

include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3

network-pulsed.stl: network-flow.stl
	ln -s -f network-flow.stl network-pulsed.stl

network-pulsed.tst: network-pulsed.stl
