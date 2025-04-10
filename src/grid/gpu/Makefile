.PHONY: all clean

CFLAGS += -I.. -O2

all: libgpu.a

# uncomment if errors.lex is modified
# errors.c: errors.lex
#	flex -P gpu_errors -o errors.c errors.lex

libgpu.a: errors.o glad.o
	ar cr libgpu.a errors.o glad.o

check:

Benchmarks/plots: Benchmarks/advection Benchmarks/bump2D-gpu Benchmarks/advection.awk \
	Benchmarks/lid Benchmarks/reversed Benchmarks/turbulence
Benchmarks.md.html: Benchmarks/plots

include ../../Makefile.defs
