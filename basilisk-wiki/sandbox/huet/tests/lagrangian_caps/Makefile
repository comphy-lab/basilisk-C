include $(BASILISK)/Makefile.defs
CFLAGS += -O2 -I../../src -Wdimensions

nh_shear.tst: CFLAGS += -DTEND=0.001

bending_shear.tst: CFLAGS += -DTEND=0.000025

# The simulation below is far too big to run on the website so it will just run for a few time steps. The results shown on the website are obtained using 96 processors for 3 days.
constricted_channel.tst: CC=mpicc -D_MPI=12 -DTEND=0.001

caps_interception.tst: CC=mpicc -D_MPI=12 -DTEND=0.00025

caps_interception_inertia.tst: CC=mpicc -D_MPI=12 -DTEND=0.001
