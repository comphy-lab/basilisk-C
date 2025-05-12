all:

include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3 -DDISPLAY=1

## 2D

cylinder-unbounded.tst: CC = mpicc -D_MPI=4

cylinder-near-wall-viscoplastic.tst: CC = mpicc -D_MPI=2
cylinder-near-wall-viscoplastic.tst: CFLAGS += -DTRACE=2 -DBMR=1 -DGAP=400 -DDLENGTH=256 -DLMAX=13

embed-vp-2D-tests:

## 3D

sphere-near-wall-viscoplastic.tst: CC = mpicc -D_MPI=4
sphere-near-wall-viscoplastic.tst: CFLAGS += -DTRACE=2 -DBMR=1 -DGAP=400 -DDLENGTH=256 -DLMAX=13

embed-vp-3D-tests:
