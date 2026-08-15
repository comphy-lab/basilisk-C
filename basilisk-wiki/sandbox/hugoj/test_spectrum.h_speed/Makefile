CFLAGS += -I$(DOCUMENT_ROOT)/sandbox
CFLAGS += -disable-dimensions 

all: speed_wave_ini.tst

speed_wave_ini.tst: CFLAGS+="-DTRACE=2 -DCHECK_PARSEVAL"

include $(BASILISK)/Makefile.defs
