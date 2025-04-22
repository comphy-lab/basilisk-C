#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "geometry.h"
#include "tiny.h"
#define sq(x) ((x)*(x))
static mat4 Viewport;

void tiny_viewport (const int x, const int y, const int w, const int h)
{
  Viewport = (mat4){{w/2., 0, 0, x + w/2.}, {0, h/2., 0, y + h/2.}, {0,0,1,0}, {0,0,0,1}};
}

/**
## Framebuffer */

framebuffer * TinyFramebuffer = NULL;

void framebuffer_destroy (framebuffer * p) {
  free (p->image);
  free (p->zbuffer);
  free (p->depth);
  free (p);
  TinyFramebuffer = NULL;
}

framebuffer * framebuffer_new (unsigned width, unsigned height)
{
  framebuffer * f = (framebuffer *) malloc (sizeof (framebuffer));
  f->image = (unsigned char *) calloc (width*height*4, sizeof (unsigned char));
  f->depth = NULL;
  f->width = width, f->height = height;
  f->zbuffer = (real *) malloc (width*height*sizeof (real));
  real * p = f->zbuffer;
  for (unsigned int i = 0; i < width*height; i++, p++)
    *p = 1e30;
  TinyFramebuffer = f;
  tiny_viewport (0., 0., width, height);
  return f;
}

unsigned char * framebuffer_image (framebuffer * p) {
  return p->image;
}

float * framebuffer_depth (framebuffer * p)
{
  if (!p->depth)
    p->depth = (float *) malloc (p->width*p->height*sizeof (float));
  real * z = p->zbuffer;
  float * d = p->depth;
  for (unsigned long i = 0; i < p->width*p->height; i++, z++, d++)
    *d = *z;
  return p->depth;
}

static inline
void framebuffer_set_depth (framebuffer * f, int x, int y, const TinyColor * color, real frag_depth)
{
  unsigned char * c = f->image + 4*(x + y*f->width);
  for (int i = 0; i < 4; i++)
    c[i] = ((unsigned char *)color)[i];
  f->zbuffer[x + y*f->width] = frag_depth;
}

/**
## Primitives */

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))
#define swap(a,b) do { typeof (a) c; c = a; a = b; b = c; } while (0)

static int constant_color (const void * data, const vec3 bar, TinyColor * frag_color)
{
  const TinyColor * color = data;
  *frag_color = *color;
  return 0;
}

static inline
real orient2d (const vec2 a, const vec2 b, const vec2 c)
{
  return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
}

void tiny_triangle (const vec4 clip_verts[3],
		    const void * shader, const TinyShader fragment,
		    const int face,
		    framebuffer * image)
{
  vec4 pts[3] = {
    mat4_mul (Viewport, clip_verts[0]),
    mat4_mul (Viewport, clip_verts[1]),
    mat4_mul (Viewport, clip_verts[2])
  };  // triangle screen coordinates before persp. division
  vec2 v[3] = {
    vec4_proj2(vec4_div (pts[0], pts[0].t)),
    vec4_proj2(vec4_div (pts[1], pts[1].t)),
    vec4_proj2(vec4_div (pts[2], pts[2].t))
  };  // triangle screen coordinates after  persp. division

  // backface culling
  real area = orient2d (v[0], v[1], v[2]);
  if (!area || face*area < 0) return;
  
  int bboxmin[2] = {100000, 100000};
  int bboxmax[2] = {- 100000, - 100000};
  for (int i = 0; i < 3; i++) {
    if (bboxmin[0] > v[i].x) bboxmin[0] = v[i].x;
    if (bboxmax[0] < v[i].x) bboxmax[0] = v[i].x;
    if (bboxmin[1] > v[i].y) bboxmin[1] = v[i].y;
    if (bboxmax[1] < v[i].y) bboxmax[1] = v[i].y;
  } // fixme: make this a function and return, also add to tiny_line 

  if (bboxmin[0] < 0) bboxmin[0] = 0;
  if (bboxmin[1] < 0) bboxmin[1] = 0;  
  if (bboxmax[0] >= image->width) bboxmax[0] = image->width - 1;
  if (bboxmax[1] >= image->height) bboxmax[1] = image->height - 1;

#if 0 // OpenMP does accelerate somewhat but not much
#pragma omp parallel for
#endif
  for (int y = bboxmin[1]; y <= bboxmax[1]; y++)
    for (int x = bboxmin[0]; x <= bboxmax[0]; x++) {
      vec3 bc; // barycentric coordinates
      if ((bc.x = orient2d (v[1], v[2], (vec2){x, y})/area) >= 0 &&
	  (bc.y = orient2d (v[2], v[0], (vec2){x, y})/area) >= 0 &&
	  (bc.z = orient2d (v[0], v[1], (vec2){x, y})/area) >= 0) {
#if 1
	// check https://github.com/ssloy/tinyrenderer/wiki/Technical-difficulties-linear-interpolation-with-perspective-deformations	
	bc = (vec3){bc.x/pts[0].t, bc.y/pts[1].t, bc.z/pts[2].t};
	bc = vec3_div (bc, bc.x + bc.y + bc.z);
#endif
	real frag_depth = vec3_scalar ((vec3){clip_verts[0].z, clip_verts[1].z, clip_verts[2].z}, bc);
	if (frag_depth >= image->zbuffer[x + y*image->width])
	  continue;
	TinyColor color;
	if (fragment (shader, bc, &color)) continue; // fragment shader can discard current fragment
	framebuffer_set_depth (image, x, y, &color, frag_depth);
      }
    }
}

void tiny_line (const vec4 clip_verts0, const vec4 clip_verts1, const TinyColor * color, float thickness,
		framebuffer * image)
{
  vec4 pts[2]  = {
    mat4_mul (Viewport, clip_verts0),
    mat4_mul (Viewport, clip_verts1)
  };  // line screen coordinates before persp. division
  vec2 v[2] = {
    vec4_proj2 (vec4_div (pts[0], pts[0].t)),
    vec4_proj2 (vec4_div (pts[1], pts[1].t))
  };  // line screen coordinates after  persp. division

  if (thickness > 1) {
    vec2 t = vec2_normalized (vec2_sub (v[1], v[0]));
    mat4 i = mat4_invert (Viewport);
    thickness /= 2.;
    float ext = 0.;
    vec4 v1 = mat4_mul (i, (vec4){ pts[0].t*(v[0].x + (- t.x*ext + t.y)*thickness),
				   pts[0].t*(v[0].y + (- t.y*ext - t.x)*thickness),
				   pts[0].z, pts[0].t });
    vec4 v2 = mat4_mul (i, (vec4){ pts[1].t*(v[1].x + (+ t.x*ext + t.y)*thickness),
				   pts[1].t*(v[1].y + (+ t.y*ext - t.x)*thickness),
				   pts[1].z, pts[1].t });
    vec4 v3 = mat4_mul (i, (vec4){ pts[1].t*(v[1].x + (+ t.x*ext - t.y)*thickness),
				   pts[1].t*(v[1].y + (+ t.y*ext + t.x)*thickness),
				   pts[1].z, pts[1].t });
    vec4 v4 = mat4_mul (i, (vec4){ pts[0].t*(v[0].x + (- t.x*ext - t.y)*thickness),
				   pts[0].t*(v[0].y + (- t.y*ext + t.x)*thickness),
				   pts[0].z, pts[0].t });
    tiny_triangle ((vec4[3]){v1, v2, v3}, color, constant_color, 0, image);
    tiny_triangle ((vec4[3]){v3, v4, v1}, color, constant_color, 0, image);
    return;
  }
  
  int x0 = v[0].x, y0 = v[0].y, x1 = v[1].x, y1 = v[1].y;
  if (x1 == x0 && y1 == y0) return;
  real z0 = clip_verts0.z, z1 = clip_verts1.z;
  // from: http://members.chello.at/~easyfilter/bresenham.html
  int dx = abs (x1 - x0), sx = x0 < x1 ? 1 : -1;
  int dy = abs (y1 - y0), sy = y0 < y1 ? 1 : -1;
  real a = (z1 - z0)/(real)(dx > dy ? (x1 - x0) : (y1 - y0));
  int err = dx - dy; /* error value e_xy */
  int x = x0, y = y0;
  while (x != x1 || y != y1) {
    real frag_depth = z0 - 0.01 + a*(dx > dy ? (x - x0) : (y - y0));
    if (x >= 0 && y >= 0 && x < image->width && y < image->height &&
	frag_depth < image->zbuffer[x + y*image->width])
      framebuffer_set_depth (image, x, y, color, frag_depth);
    int e2 = 2*err;
    if (e2 >= - dy) { err -= dy; x += sx; } /* e_xy+e_x > 0 */
    if (e2 <= dx) { err += dx; y += sy; } /* e_xy+e_y < 0 */
  }
}

void tiny_point (const vec4 clip_verts0, const TinyColor * color, float radius,
		 framebuffer * image) {
  vec4 a  =  mat4_mul (Viewport, clip_verts0);
  // point screen coordinates before persp. division
  vec2 b = vec4_proj2 (vec4_div (a, a.t));
  // line screen coordinates after  persp. division
  for (real x = b.x - radius; x <= b.x + radius; x++) 
    for (real y = b.y - radius; y <= b.y + radius; y++) 
      if (sq(x - b.x) + sq(y - b.y) < sq(radius) &&
	  x >= 0 && y >= 0 && x < image->width &&  y < image->height &&
	  clip_verts0.z < image->zbuffer[(int)x + (int)y*image->width])
	framebuffer_set_depth (image, x, y, color, clip_verts0.z);
}
