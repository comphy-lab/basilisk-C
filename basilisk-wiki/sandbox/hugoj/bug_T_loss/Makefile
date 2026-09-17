CFLAGS += -I$(DOCUMENT_ROOT)/sandbox
CFLAGS += -O2 -disable-dimensions # skip dimensions check for now

all: T_cst_withwave.cuda.tst T_cst_nowave.cuda.tst

T_cst_withwave.c: ml_breaking_strat.c
	ln -sf ml_breaking_strat.c T_cst_withwave.c

T_cst_withwave.tst: CFLAGS += -DWITH_PM

T_cst_withwave.cuda.tst: CFLAGS += -DWITH_PM

T_cst_nowave.c: ml_breaking_strat.c
	ln -sf ml_breaking_strat.c T_cst_nowave.c


include $(BASILISK)/Makefile.defs
