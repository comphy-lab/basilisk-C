CFLAGS += -O2 -disable-dimensions
include $(BASILISK)/Makefile.defs

capwave.tst: capclsvof.tst caphf.tst caphf2d.tst
#        @echo "Running capwave.tst after dependencies"

capclsvof.c: capwave.c
	ln -s capwave.c capclsvof.c
capclsvof.s: CFLAGS += -DCLSVOF=1
capclsvof.tst: CFLAGS += -DCLSVOF=1

caphf.c: capwave.c
	ln -s capwave.c caphf.c
caphf.s: CFLAGS += -DINTHF=1
caphf.tst: CFLAGS += -DINTHF=1

caphf2d.c: capwave.c
	ln -s capwave.c caphf2d.c
caphf2d.s: CFLAGS += -DCLSVOF=1 -DHF2D=1
caphf2d.tst: CFLAGS += -DCLSVOF=1 -DHF2D=1

rising.tst: rising-clsvof.tst rising-hf.tst rising-hf2d.tst
#	@echo "Running rising.tst after dependencies"

rising-clsvof.c: rising.c
	ln -s rising.c rising-clsvof.c
rising-clsvof.s: CFLAGS += -DCLSVOF=1
rising-clsvof.tst: CFLAGS += -DCLSVOF=1

rising-hf.c: rising.c
	ln -s rising.c rising-hf.c
rising-hf.s: CFLAGS += -DINTHF=1
rising-hf.tst: CFLAGS += -DINTHF=1

rising-hf2d.c: rising.c
	ln -s rising.c rising-hf2d.c
rising-hf2d.s: CFLAGS += -DHF2D=1
rising-hf2d.tst: CFLAGS += -DHF2D=1

marangoni.tst: marangoni-clsvof.tst marangoni-hf2d.tst marangoni-PF.tst

marangoni-PF.c: marangoni.c
	ln -s marangoni.c marangoni-PF.c
marangoni-PF.s: CFLAGS += -DPALAS=1
marangoni-PF.tst: CFLAGS += -DPALAS=1

marangoni-clsvof.c: marangoni.c
	ln -s marangoni.c marangoni-clsvof.c
marangoni-clsvof.s: CFLAGS += -DCLSVOF=1
marangoni-clsvof.tst: CFLAGS += -DCLSVOF=1

marangoni-hf2d.c: marangoni.c
	ln -s marangoni.c marangoni-hf2d.c
marangoni-hf2d.s: CFLAGS += -DHF2D=1
marangoni-hf2d.tst: CFLAGS += -DHF2D=1

marangoniwall.tst: marangoniwall-clsvof.tst marangoniwall-hf2d.tst marangoniwall2.tst marangoniwall-clsvof2.tst marangoniwall-hf2d2.tst

marangoniwall-clsvof.c: marangoniwall.c
	ln -s marangoniwall.c marangoniwall-clsvof.c
marangoniwall-clsvof.s: CFLAGS += -DCLSVOF=1
marangoniwall-clsvof.tst: CFLAGS += -DCLSVOF=1

marangoniwall-hf2d.c: marangoniwall.c
	ln -s marangoniwall.c marangoniwall-hf2d.c
marangoniwall-hf2d.s: CFLAGS += -DHF2D=1
marangoniwall-hf2d.tst: CFLAGS += -DHF2D=1

marangoniwall2.c: marangoniwall.c
	ln -s marangoniwall.c marangoniwall2.c
marangoniwall2.s: CFLAGS += -DCASE2=1
marangoniwall2.tst: CFLAGS += -DCASE2=1

marangoniwall-clsvof2.c: marangoniwall.c
	ln -s marangoniwall.c marangoniwall-clsvof2.c
marangoniwall-clsvof2.s: CFLAGS += -DCLSVOF=1 -DCASE2=1
marangoniwall-clsvof2.tst: CFLAGS += -DCLSVOF=1 -DCASE2=1

marangoniwall-hf2d2.c: marangoniwall.c
	ln -s marangoniwall.c marangoniwall-hf2d2.c
marangoniwall-hf2d2.s: CFLAGS += -DHF2D=1 -DCASE2=1
marangoniwall-hf2d2.tst: CFLAGS += -DHF2D=1 -DCASE2=1

marangoniad2bub.tst: marangoniad2bub-clsvof.tst

marangoniad2bub-clsvof.c: marangoniad2bub.c
	ln -s marangoniad2bub.c marangoniad2bub-clsvof.c
marangoniad2bub-clsvof.s: CFLAGS += -DCLSVOF=1
marangoniad2bub-clsvof.tst: CFLAGS += -DCLSVOF=1

spuriousMov.tst: spuriousMov-hf2d.tst spuriousMov-clsvof.tst

spuriousMov-hf2d.c: spuriousMov.c
	ln -sf spuriousMov.c spuriousMov-hf2d.c
spuriousMov-hf2d.s: CFLAGS += -DHF2D=1
spuriousMov-hf2d.tst: CFLAGS += -DHF2D=1

spuriousMov-clsvof.c: spuriousMov.c
	ln -sf spuriousMov.c spuriousMov-clsvof.c
spuriousMov-clsvof.s: CFLAGS += -DCLSVOF=1
spuriousMov-clsvof.tst: CFLAGS += -DCLSVOF=1
