include $(BASILISK)/Makefile.defs

#BASILISK_INCLUDES ?= -I../src
#BASILISK_LIBS ?= 

#CFLAGS += -O2 -disable-dimensions $(BASILISK_INCLUDES) $(BASILISK_LIBS)
CFLAGS += -O2 -disable-dimensions -I../src

poisson-gfm.tst: CFLAGS += -DGFM=1
capwave.tst: CFLAGS += -DGFM=1
spurious.tst: CFLAGS += -DGFM=1
oscillation.tst: CFLAGS += -DGFM=1

pinning-axi.c: pinning.c
	ln -s pinning.c pinning-axi.c
pinning-axi.s: CFLAGS += -DAXI=1
pinning-axi.tst: CFLAGS += -DAXI=1

pinning.tst: pinning-axi.tst

# New targets
GFM = poisson-gfm.tst capwave.tst spurious.tst oscillation.tst
quick: aslam.tst aslamvof.tst fsolve-gsl.tst interfaceregression.tst \
  mapregion.tst radiation.tst spark.tst fickcorrected.tst \
  boilingtemperature.tst $(GFM)
gfm: $(GFM)
slow: pinning.tst
opensmoke-all: opensmoke.tst opensmoke-properties.tst batchreactor.tst
all: quick slow opensmoke-all