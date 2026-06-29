.PHONY: all clean

CFLAGS += -I$(BASILISK) -I.. -DSINGLE_PRECISION -O2

# By default we only make liberrors.a because it does not depend on
# external libraries, unlike libgpu.a

all: liberrors.a

# these are not Basilisk programs
EXCLUDE = errors.c  glad.c  opengl.c  reduction.c

# uncomment if errors.lex is modified
# errors.c: errors.lex
#	flex -P gpu_errors -o errors.c errors.lex

libgpu.a: glad.o opengl.o
	ar cr libgpu.a glad.o opengl.o

liberrors.a: errors.o
	ar cr liberrors.a errors.o

opengl.o: backend.h ../externals.h opengl.c reduction.c

check:

Benchmarks/plots: Benchmarks/advection Benchmarks/bump2D-gpu Benchmarks/advection.awk \
	Benchmarks/lid Benchmarks/reversed Benchmarks/turbulence
Benchmarks.md.html: Benchmarks/plots

include ../../Makefile.defs
