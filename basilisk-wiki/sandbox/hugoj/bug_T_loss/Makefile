CFLAGS += -I$(DOCUMENT_ROOT)/sandbox
CFLAGS += -O2 -disable-dimensions # skip dimensions check for now

all: T_cst_withwave.cuda.tst T_cst_nowave.cuda.tst

T_cst_withwave.c: PMwave.c
	ln -sf PMwave.c T_cst_withwave.c

T_cst_withwave.tst: CFLAGS += -DWITH_PM

T_cst_withwave.cuda.tst: CFLAGS += -DWITH_PM

T_cst_nowave.c: PMwave.c
	ln -sf PMwave.c T_cst_nowave.c


include $(BASILISK)/Makefile.defs
