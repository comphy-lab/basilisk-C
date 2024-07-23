#include <GL/osmesa.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "framebuffer.h"

struct _framebuffer {
  OSMesaContext ctx;
  unsigned char * image;
  float * depth;
  unsigned width, height;
};

void framebuffer_destroy (framebuffer * p)
{
  free (p->image);
  free (p->depth);
  OSMesaDestroyContext (p->ctx);
  free (p);
}

framebuffer * framebuffer_new (unsigned width, unsigned height)
{
  framebuffer * p = malloc (sizeof(framebuffer));

  p->width = width;
  p->height = height;
  p->image = malloc (p->width*p->height*4*sizeof (char));
  p->depth = NULL;
  
  /* Create an RGBA-mode context for OSMesa */
#if 1 // OSMESA_MAJOR_VERSION * 100 + OSMESA_MINOR_VERSION >= 305
  /* specify Z, stencil, accum sizes */
  p->ctx = OSMesaCreateContextExt (OSMESA_RGBA, 32, 0, 0, NULL);
#else
  p->ctx = OSMesaCreateContext (OSMESA_RGBA, NULL);
#endif
  if (!p->ctx) {
    fprintf (stderr, "framebuffer_new(): OSMesaCreateContext failed!\n");
    exit (1);
  }

  if (!OSMesaMakeCurrent (p->ctx, p->image, GL_UNSIGNED_BYTE,
			  p->width, p->height)) {
    fprintf (stderr, "framebuffer_new(): OSMesaMakeCurrent failed!\n");
    exit (1);
  }

  return p;
}

unsigned char * framebuffer_image (framebuffer * p)
{
  return p->image;
}

float * framebuffer_depth (framebuffer * p)
{
  unsigned int * depth;
  GLint width, height, bytesPerValue;
  OSMesaGetDepthBuffer (p->ctx, &width, &height, &bytesPerValue,
			(void **)&depth);
  assert (p->width == width && p->height == height && bytesPerValue == 4);
#if GALLIUM
  // fix for bug in gallium/libosmesa
  // the depth buffer is flipped vertically
  // see https://gitlab.freedesktop.org/mesa/mesa/-/issues/885
  // Fixed since
  GLint i, j;
  for (j = 0; j < height/2; j++)
    for (i = 0; i < width; i++) {
      unsigned int tmp = depth[j*width + i];
      depth[j*width + i] = depth[(height - 1 - j)*width + i];
      depth[(height - 1 - j)*width + i] = tmp;
    }
#endif // GALLIUM
  if (!p->depth)
    p->depth = (float *) malloc (p->width*p->height*sizeof (float));
  float * d = p->depth;
  unsigned int * z = depth;
  for (unsigned long i = 0; i < p->width*p->height; i++, z++, d++)
    *d = *z;
  return p->depth;
}
