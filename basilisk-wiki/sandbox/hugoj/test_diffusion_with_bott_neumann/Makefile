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

T_cst_delta.c: diff_test.c
	ln -sf diff_test.c T_cst_delta.c
T_cst_delta.tst: CFLAGS += -DNEUMANN0 -DT_CST -DDIFF_DELTA
T_cst_delta.cuda.tst: CFLAGS += -DNEUMANN0 -DT_CST -DDIFF_DELTA

T_cst_notperiodic_delta.c: diff_test.c
	ln -sf diff_test.c T_cst_notperiodic_delta.c
T_cst_notperiodic_delta.tst: CFLAGS += -DNEUMANN0 -DT_CST -DNOT_PERIODIC -DDIFF_DELTA
T_cst_notperiodic_delta.cuda.tst: CFLAGS += -DNEUMANN0 -DT_CST -DNOT_PERIODIC -DDIFF_DELTA

# isolated vertical-diffusion test (see simple_iso.c)
simple_iso.c: diff_test.c
simple_iso_delta.c: simple_iso.c
	ln -sf simple_iso.c simple_iso_delta.c
simple_iso_delta.tst: CFLAGS += -DDIFF_DELTA
simple_iso_delta.cuda.c: simple_iso.c
	ln -sf simple_iso.c simple_iso_delta.cuda.c
simple_iso_delta.cuda.tst: CFLAGS += -DDIFF_DELTA
simple_iso_single.c: simple_iso.c
	ln -sf simple_iso.c simple_iso_single.c
simple_iso_single.tst: CFLAGS += -DSINGLE_PRECISION

# S1-style isolation of the stock vertical_diffusion() (simple.c, report1.md)
simple_rep.c: simple.c
	ln -sf simple.c simple_rep.c
simple_rep.tst: simple_rep.c
simple_rep.tst: CFLAGS += -DNOREP
simple_rep_single.c: simple.c
	ln -sf simple.c simple_rep_single.c
simple_rep_single.tst: CFLAGS += -DNOREP -DSINGLE_PRECISION
simple_rep.cuda.tst: CFLAGS += -DNOREP

# minimal coupled T_cst analogues (simple.c, N=1)
simple_adv.c: simple.c
	ln -sf simple.c simple_adv.c
simple_adv.tst: simple_adv.c
simple_adv.tst: CFLAGS += -DTRACER -DNO_REMAP
simple_remap.c: simple.c
	ln -sf simple.c simple_remap.c
simple_remap.tst: simple_remap.c
simple_remap.tst: CFLAGS += -DTRACER
simple_remap.cuda.tst: CFLAGS += -DTRACER
simple_hr.c: simple.c
	ln -sf simple.c simple_hr.c
simple_hr.tst: simple_hr.c
simple_hr.tst: CFLAGS += -DTRACER -DH_OVER_NL
simple_hr.cuda.tst: simple_hr.c
simple_hr.cuda.tst: CFLAGS += -DTRACER -DH_OVER_NL
# full T_cst-like parameters: zb=-100, T=20, G, CFL, h=100/nl
simple_tcl.c: simple.c
	ln -sf simple.c simple_tcl.c
simple_tcl.tst: simple_tcl.c
simple_tcl.tst: CFLAGS += -DTRACER -DH_OVER_NL -DT_CST_LIKE
simple_tcl.cuda.tst: simple_tcl.c
simple_tcl.cuda.tst: CFLAGS += -DTRACER -DH_OVER_NL -DT_CST_LIKE

# drift-attribution diagnostics (phase C)
diag.c: diff_test.c
	ln -sf diff_test.c diag.c
diag.tst: CFLAGS += -DNEUMANN0 -DT_CST -DDIAG_STEPS
diag_delta.c: diag.c
	ln -sf diff_test.c diag_delta.c
diag_delta.tst: CFLAGS += -DNEUMANN0 -DT_CST -DDIAG_STEPS -DDIFF_DELTA
diag.cuda.tst: CFLAGS += -DNEUMANN0 -DT_CST -DDIAG_STEPS
diag_delta.cuda.c: diag.c
	ln -sf diff_test.c diag_delta.cuda.c
diag_delta.cuda.tst: CFLAGS += -DNEUMANN0 -DT_CST -DDIAG_STEPS -DDIFF_DELTA

# minimal coupled run: hydro + dr only (no nh, no remap) -- culprit isolation
diag_min.c: diff_test.c
	ln -sf diff_test.c diag_min.c
diag_min.tst: CFLAGS += -DNEUMANN0 -DT_CST -DDIAG_STEPS -DMINIMAL
diag_min.cuda.tst: CFLAGS += -DNEUMANN0 -DT_CST -DDIAG_STEPS -DMINIMAL

# nh only (no remap)
diag_nh.c: diff_test.c
	ln -sf diff_test.c diag_nh.c
diag_nh.tst: CFLAGS += -DNEUMANN0 -DT_CST -DDIAG_STEPS -DNO_REMAP
diag_nh.cuda.tst: CFLAGS += -DNEUMANN0 -DT_CST -DDIAG_STEPS -DNO_REMAP

include $(BASILISK)/Makefile.defs


