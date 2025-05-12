all:	embed-stokes-torques-tests

include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -Wdimensions -DMTRACE=3 -DDISPLAY=-1 -I$(DOCUMENT_ROOT)/sandbox

#embed-stokes-torques-tests: couette_torque.tst
embed-stokes-torques-tests: couette_torque_space_shift.tst
