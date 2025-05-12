# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3

include $(BASILISK)/Makefile.defs

Poisson-9Point-UniformGrid-O4.c: Poisson-9Point-UniformGrid-O2.c
	ln -s Poisson-9Point-UniformGrid-O2.c Poisson-9Point-UniformGrid-O4.c
Poisson-9Point-UniformGrid-O4.tst: CFLAGS += -DHIGHER=1

Poisson-9Point-UniformGrid-O2.tst: Poisson-9Point-UniformGrid-O4.tst

