CFLAGS += -I$(DOCUMENT_ROOT)/sandbox 
CFLAGS += -disable-dimensions

all: monochromatic.tst


monochromatic.c: gen_ini.c
	ln -sf gen_ini.c monochromatic.c

monochromatic.tst: CFLAGS += -DWAVE=0

stokes.c: gen_ini.c
	ln -sf gen_ini.c stokes.c

stokes.tst: CFLAGS += -DWAVE=1

synth_eta.c: gen_ini.c
	ln -sf gen_ini.c synth_eta.c

synth_eta.tst: CFLAGS += -DWAVE=2


gen_ini.tst: CC = mpicc -D_MPI=16



	# all: comparison_Fk.svg
#
# comparison_Fk.svg: F_k_py F_k_C
# 	python3 compare.py
#
#
# F_k_py: gen_ini.py
# 	python3 gen_ini.py
#
# F_k_C: main
# 	./main
#
# main: main.c spectrum.h
# 	qcc -autolink -O2 -I/home/jacqhugo/basilisk/wiki/sandbox -disable-dimensions -g -Wall -pipe -D_FORTIFY_SOURCE=2 -o main main.c -lm
#
# clean:
# 	rm F_* kx* ky* kmod* eta_C eta_py main *png *svg

include $(BASILISK)/Makefile.defs
