CFLAGS += -I..

.PHONY: tinyrenderer all clean

all: tinyrenderer libglutils.a libfb_tiny.a

clean: tinyrenderer
tinyrenderer:
	$(MAKE) -C $@ $(filter-out tinyrenderer,$(MAKECMDGOALS))

libglutils.a: trackball.o utils.o polygonize.o \
		og_stroke_mono_roman.o parser.o TinyPngOut.o
	ar cr $@ $^

utils.o: utils.h

# Uncomment the recipe below if you need to re-generate parser.tab.c

# parser.tab.c: parser.y
#	bison parser.y

parser: parser.tab.c parser.h
	$(CC) $(CFLAGS) -DSTANDALONE=1 -o parser parser.tab.c -lm

parser.o: parser.tab.c parser.h
	$(CC) $(CFLAGS) -c -o parser.o parser.tab.c

clean:
	rm -f *.o *.a

libfb_tiny.a: fb_tiny.o tinyrenderer/tiny.o
	ar cr $@ $^

tinyrenderer/tiny.o:
	cd tinyrenderer && $(MAKE)

fb_tiny.o: tinygl.h tinyrenderer/tiny.h tinyrenderer/geometry.h 

## These libraries depend on OpenGL and are not built by default

libfb_osmesa.a: fb_osmesa.o
	ar cr $@ $^

libfb_glx.a: fb_glx.o OffscreenContextGLX.o fbo.o
	ar cr $@ $^

check:
