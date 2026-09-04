CFLAGS += -I$(DOCUMENT_ROOT)/sandbox -O2
CFLAGS += -disable-dimensions # skip dimensions check for now

# WAVE DECAY
PM_decay.c: PM_seastate.c
	ln -sf PM_seastate.c PM_decay.c

PM_decay.cuda.tst: CFLAGS += -DDUMPS=1

# CONSTANT FORCING
PM_seastate.tst: CFLAGS += -DFORCING=1 -DDUMPS=1

PM_seastate.cuda.tst: CFLAGS += -DDUMPS=1 -DFORCING=1 

# CONSTANT ENERGY
PM_seastate_hold.c: PM_seastate.c
	ln -sf PM_seastate.c PM_seastate_hold.c

PM_seastate_hold.tst: CFLAGS += -DFORCING=1 -DCST_NRJ=1 -DDUMPS=1

PM_seastate_hold.cuda.tst: CFLAGS += -DFORCING=1 -DCST_NRJ=1 -DDUMPS=1

include $(BASILISK)/Makefile.defs


