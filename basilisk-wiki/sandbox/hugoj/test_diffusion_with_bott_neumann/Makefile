CFLAGS += -I$(DOCUMENT_ROOT)/sandbox -O2
CFLAGS += -disable-dimensions # skip dimensions check for now

all: T_cst_notperiodic.tst T_cst_notperiodic.cuda.tst T_cst.cuda.tst T_cst.tst with_diff.tst diff_test_H.tst neumann0.tst diff_test.tst diff_test/plots

diff_test_H.c: diff_test.c
	ln -sf diff_test.c diff_test_H.c
diff_test_H.tst: CFLAGS+= -DHEATING

neumann0.c: diff_test.c
	ln -sf diff_test.c neumann0.c
neumann0.tst: CFLAGS+= -DNEUMANN0

with_diff.c: diff_test.c
	ln -sf diff_test.c with_diff.c
with_diff.tst: CFLAGS+= -DNEUMANN0 -DHDIFF

with_diff.cuda.tst: CFLAGS+= -DNEUMANN0 -DHDIFF

T_cst.c: diff_test.c
	ln -sf diff_test.c T_cst.c
T_cst.tst: CFLAGS += -DNEUMANN0 -DT_CST
T_cst.cuda.tst: CFLAGS += -DNEUMANN0 -DT_CST


T_cst_notperiodic.c: diff_test.c
	ln -sf diff_test.c T_cst_notperiodic.c
T_cst_notperiodic.tst: CFLAGS += -DNEUMANN0 -DT_CST -DNOT_PERIODIC
T_cst_notperiodic.cuda.tst: CFLAGS += -DNEUMANN0 -DT_CST -DNOT_PERIODIC


include $(BASILISK)/Makefile.defs


