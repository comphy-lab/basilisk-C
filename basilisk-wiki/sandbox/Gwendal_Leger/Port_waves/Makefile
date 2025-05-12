MPI = 1
ifeq ($(MPI), 0)
	CFLAGS += -O2 -fopenmp -events
else
#	CC = mpicc -D_MPI=$(MPI)
	CFLAGS += -O2 -DDISPLAY=1 -fopenmp
endif


include $(BASILISK)/Makefile.defs




MAIN = Port_waves
MAIN2 = post_treatment


# Removes EVERYTHING deemed "useless" excepts scripts
fullclean :
	rm *~ $(MAIN)/* $(MAIN).tst $(MAIN).s $(MAIN).s.d Makefile.deps Makefile.tests $(MAIN).c.html $(MAIN).c.tags $(MAIN).c.tags.d \
	$(MAIN2)/* $(MAIN2).tst $(MAIN2).s $(MAIN2).s.d $(MAIN2).c.html $(MAIN2).c.tags $(MAIN2).c.tags.d


clean_post_treatment :
	rm *~ $(MAIN2)/* $(MAIN2).tst $(MAIN2).s $(MAIN2).s.d $(MAIN2).c.html $(MAIN2).c.tags $(MAIN2).c.tags.d

