CFLAGS += -I$(DOCUMENT_ROOT)/sandbox -O2
CFLAGS += -disable-dimensions # skip dimensions check for now

# WAVE DECAY
PM_decay.c: PM_seastate.c
	ln -sf PM_seastate.c PM_decay.c

PM_decay.cuda.tst: 

# CONSTANT FORCING
PM_seastate.tst: CFLAGS += -DCST_FORCING=1 

PM_seastate.cuda.tst: CFLAGS += -DCST_FORCING=1 

# CONSTANT ENERGY
PM_seastate_hold.c: PM_seastate.c
	ln -sf PM_seastate.c PM_seastate_hold.c

PM_seastate_hold.tst: CFLAGS += -DHOLD_FORCING=1

PM_seastate_hold.cuda.tst: CFLAGS += -DHOLD_FORCING=1 

include $(BASILISK)/Makefile.defs


