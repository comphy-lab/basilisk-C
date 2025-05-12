# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3 -fopenmp

include $(BASILISK)/Makefile.defs

b-flood.tgz:
	cd .. && tar czvf b-flood/b-flood.tgz \
		b-flood/Readme b-flood/Download \
		b-flood/*.h b-flood/Testcases/*.c b-flood/Testcases/*.ref \
		b-flood/*.c b-flood/*.ref \
		b-flood/ToceFluvial/hydro.asc \
		b-flood/ToceFluvial/*.ref \
		b-flood/ToceFluvial/Topo
