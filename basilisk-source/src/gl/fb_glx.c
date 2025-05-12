#include <stdlib.h>
#include "framebuffer.h"
#include "OffscreenContext.h"

struct _framebuffer {
  OffscreenContext * ctx;
  float * depth;
  unsigned width, height;
};

void framebuffer_destroy (framebuffer * p)
{
  free (p->depth);
  teardown_offscreen_context (p->ctx);
  free (p);
}

framebuffer * framebuffer_new (unsigned width, unsigned height)
{
  framebuffer * p = malloc (sizeof (struct _framebuffer));
  p->ctx = create_offscreen_context (width, height);
  p->depth = NULL;
  p->width = width, p->height = height;
  return p;
}

unsigned char * framebuffer_image (framebuffer * p)
{
  return get_framebuffer_pixels (p->ctx);
}

float * framebuffer_depth (framebuffer * p)
{
  unsigned int * depth = get_framebuffer_depth (p->ctx);
  if (!p->depth)
    p->depth = (float *) malloc (p->width*p->height*sizeof (float));
  float * d = p->depth;
  unsigned int * z = depth;
  for (unsigned long i = 0; i < p->width*p->height; i++, z++, d++)
    *d = *z;
  return p->depth;
}
