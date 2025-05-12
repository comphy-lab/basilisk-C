# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2

include $(BASILISK)/Makefile.defs

poiseuille-eforcing.c: poiseuille.c
	ln -s poiseuille.c poiseuille-eforcing.c
poiseuille-eforcing.tst: CFLAGS += -DPOISEUILLE=1 -DDLM_alpha_coupling=1 -DSTENCIL_INTERIOR=1 -DDLMFD=1

poiseuille.tst: poiseuille-eforcing.tst

poiseuille-temporal-dlmfd.c: poiseuille-temporal.c
	ln -s poiseuille-temporal.c poiseuille-temporal-dlmfd.c
poiseuille-temporal-dlmfd.tst: CFLAGS += -DPOISEUILLE=1 -DDLMFD=1 -DSTENCIL_INTERIOR=1 

poiseuille-temporal.tst: poiseuille-temporal-dlmfd.tst

poiseuille-temporal-eforcing.c: poiseuille-temporal.c
	ln -s poiseuille-temporal.c poiseuille-temporal-eforcing.c
poiseuille-temporal-eforcing.tst: CFLAGS += -DPOISEUILLE=1 -DDLM_alpha_coupling=1 -DSTENCIL_INTERIOR=1 -DDLMFD=1

poiseuille-temporal.tst: poiseuille-temporal-eforcing.tst
