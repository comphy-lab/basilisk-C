include ../Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -fopenmp

all: check

madsen-sv.c: madsen.c
	ln -s madsen.c madsen-sv.c 
madsen-sv.tst: madsen.s
madsen-sv.tst: CFLAGS += -DSAINT_VENANT=1

madsen.tst: madsen-sv.tst

tohoku-gn.c: tohoku.c
	ln -s tohoku.c tohoku-gn.c 
tohoku-gn.tst: tohoku.s
tohoku-gn.tst: CFLAGS += -DGN=1

tohoku-hydro.c: tohoku.c
	ln -s tohoku.c tohoku-hydro.c 
tohoku-hydro.tst: tohoku.s
tohoku-hydro.tst: CFLAGS += -DHYDRO=1

shoal-ml.c: shoal.c
	ln -s shoal.c shoal-ml.c 
shoal-ml.tst: shoal.s
shoal-ml.tst: CFLAGS += -DML=1

# Coupled solvers

master.s: slave.o

slave.o: $(QCC) slave.c ../slave.h Makefile
	$(QCC) $(CFLAGS) -fno-common -D_OBJECT -c slave.c
	objcopy -G slave_step -G slave_stop -G slave_interpolate slave.o
