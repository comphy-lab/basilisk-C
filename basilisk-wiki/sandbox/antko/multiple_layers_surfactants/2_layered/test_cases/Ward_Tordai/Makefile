CFLAGS += -O2 -events
include $(BASILISK)/Makefile.defs

Ward1.c: Ward.c
	ln -s Ward.c Ward1.c 
Ward1.tst: CFLAGS += -DSIGMA=1

Ward5.c: Ward.c
	ln -s Ward.c Ward5.c 
Ward5.tst: CFLAGS += -DSIGMA=5
  
Ward10.c: Ward.c
	ln -s Ward.c Ward10.c 
Ward10.tst: CFLAGS += -DSIGMA=10
  
Ward.tst: Ward1.tst Ward5.tst Ward10.tst

