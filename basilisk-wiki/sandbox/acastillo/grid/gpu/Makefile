.PHONY: all clean

BASE_CFLAGS = -I$(BASILISK) -I$(BASILISK)/grid -I.. -O2 -g -Wall

all: libgpu.a libgpu64.a

libgpu.a: glad.o opengl.o
	ar cr libgpu.a glad.o opengl.o

libgpu64.a: glad.o opengl64.o
	ar cr libgpu64.a glad.o opengl64.o

glad.o: $(BASILISK)/grid/gpu/glad.c
	$(CC) $(BASE_CFLAGS) -c -o $@ $<

opengl.o: opengl.c reduction.c
	$(CC) $(BASE_CFLAGS) -DSINGLE_PRECISION=1 -c -o $@ $<

opengl64.o: opengl.c reduction.c
	$(CC) $(BASE_CFLAGS) -DSINGLE_PRECISION=0 -DDOUBLE_PRECISION=1 -c -o $@ $<

clean:
	rm -f glad.o opengl.o opengl64.o libgpu.a libgpu64.a
