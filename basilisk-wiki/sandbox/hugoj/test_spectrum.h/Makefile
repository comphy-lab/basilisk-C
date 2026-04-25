CFLAGS += -I$(DOCUMENT_ROOT)/sandbox
CFLAGS += -disable-dimensions 

all: direction.tst test_parceval


test_direction: direction.tst
	python3 direction.py 

test_parceval: parceval.tst
	python3 parceval.py 

include $(BASILISK)/Makefile.defs
