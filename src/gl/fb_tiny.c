/**
# Tiny OpenGL implementation

This provides the minimal set of OpenGL functions necessary for
Basilisk View. It has no dependencies other than the standard C
library. 

The "hardware" for this implementation is the [tiny
renderer](tinyrenderer/README.md) originally written by Dmitry
V. Sokolov. 

Also this page:
https://fgiesen.wordpress.com/2013/02/06/the-barycentric-conspirac/
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "tinygl.h"
#include "tinyrenderer/geometry.h"
#include "tinyrenderer/tiny.h"

/**
## "Unused" OpenGL functions */

void glBindTexture (GLenum target, GLuint texture) {
  // fixme: add "undefined" message option
}

void glDisable (GLenum cap) {}
void glEnable (GLenum cap) {}
void glFinish (void) {}
void glGetDoublev (GLenum pname, GLdouble * params) {}
void glHint (GLenum target, GLenum mode) {}
void glLightModeli (GLenum pname, GLint param) {}
void glOrtho (GLdouble left, GLdouble right, GLdouble bottom, GLdouble top,
	      GLdouble nearVal, GLdouble farVal) {}
void glShadeModel (GLenum mode) {}
void glTexCoord2f (GLfloat s, GLfloat t) {}
void glTexParameteri (GLenum target, GLenum pname, GLint param) {}

#define TEXTURE_WIDTH 256
static float texture[3*TEXTURE_WIDTH];

void glTexImage1D (GLenum target, GLint level, GLint internalFormat,
		   GLsizei width, GLint border, GLenum format, GLenum type,
		   const void * data)
{
  assert (target == GL_TEXTURE_1D && level == 0 && internalFormat == GL_RGB &&
	  width == TEXTURE_WIDTH && border == 0 && format == GL_RGB && type == GL_FLOAT);
  memcpy (texture, data, 3*TEXTURE_WIDTH*sizeof (float));
}

static real modelview0[16] = {
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1
}, * modelview = modelview0;

static real projection0[16] = {
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1
}, * projection = projection0;;

/**
## Lights 

Very loosely based on the Mesa implementation in
mesa-17.2.4/src/mesa/tnl/t_vb_lighttmp.h:light_fast_rgba_single(). Note
that we allow only for "grayscale" materials (i.e. the material
properties are floats not vec3). This could be changed easily if
necessary. 

A single light (GL_LIGHT0) is allowed and specular reflection is not
implemented. */

typedef struct {
  float ambient, diffuse;
  int two_sides;
  vec4 position;
  vec3 _VP_inf_norm; // in modelview coordinates  
} light_t;

static
light_t Light0 = {
  .ambient = 0.2, .diffuse = 1., // OpenGL defaults
  .two_sides = 1, // fixme: there seems to be an orientation problem compared to OSMesa when this is unset
  .position = {0, 0, 50, 0} // not OpenGL default
};

static const int not_implemented = 0;

void glLightfv (GLenum light, GLenum pname, const GLfloat *params)
{
  assert (light == GL_LIGHT0); // only one light is implemented
  switch (pname) {
    
  case GL_POSITION:
    Light0.position = (vec4){ params[0], params[1], params[2], params[3] };
    Light0._VP_inf_norm = vec3_normalized (vec4_proj3 (mat4_mul (*((mat4 *)modelview), Light0.position)));
    break;

#if 0 // fixme: does not seem to match with OSMesa when changed from the default (0.2) above   
  case GL_AMBIENT:
    assert (params[0] == params[1] && params[1] == params[2]); // only white lights are implemented
    Light0.ambient = params[0];
    break;
#endif
    
  case GL_DIFFUSE:
    assert (params[0] == params[1] && params[1] == params[2]); // only white lights are implemented
    Light0.diffuse = params[0];
    break;
    
  default:
    assert (not_implemented); // only the attributes above are implemented
  }
}

#define clamp(a,b,c) ((a) < (b) ? (b) : (a) > (c) ? (c) : (a))

static inline
float normal_shade (const vec3 normal)
{
  float n_dot_VP = vec3_scalar (normal, Light0._VP_inf_norm); // normal . light direction
  float sum = Light0.ambient;
  if (n_dot_VP < 0) {
    if (Light0.two_sides)
      n_dot_VP = - n_dot_VP;
    else
      n_dot_VP = 0.;
  }
  sum += n_dot_VP*Light0.diffuse;
  return sum > 1. ? 1. : sum;
}

/**
## Matrix transformations */

static real * current = projection0;

GLenum glGetError (void) { return GL_NO_ERROR; }

void glGetFloatv (GLenum pname, GLfloat * params)
{
  switch (pname) {

  case GL_MODELVIEW_MATRIX:
    for (int i = 0; i < 16; i++)
      params[i] = modelview[i];
    Light0._VP_inf_norm = vec3_normalized (vec4_proj3 (mat4_mul (*((mat4 *)modelview), Light0.position)));
    break;
    
  case GL_PROJECTION_MATRIX:
    for (int i = 0; i < 16; i++)
      params[i] = projection[i];
    break;

  default:
    assert (not_implemented);
  }
}

void glMultMatrixf (const GLfloat * m)
{
  GLfloat r[16];
  int i, j, k;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++) {
      GLfloat a = 0.;
      for (k = 0; k < 4; k++)
	a += current[4*k + i]*m[4*j + k];
      r[4*j + i] = a;
    }
  for (i = 0; i < 16; i++)
    current[i] = r[i];
}

void glLoadIdentity (void)
{
  static GLdouble identity[16] = {
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
  };
  int i;
  for (i = 0; i < 16; i++)
    current[i] = identity[i];
}

void glLoadMatrixd (const GLdouble * m)
{
  int i;
  for (i = 0; i < 16; i++)
    current[i] = m[i];
}

void glScalef (GLfloat x, GLfloat y, GLfloat z)
{
  glMultMatrixf ((GLfloat []){
      x, 0, 0, 0,
      0, y, 0, 0,
      0, 0, z, 0,
      0, 0, 0, 1
  });
}

void glTranslatef (GLfloat x, GLfloat y, GLfloat z)
{
  glMultMatrixf ((GLfloat []){
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      x, y, z, 1
  });
}

void glRotatef (GLfloat angle, GLfloat x, GLfloat y, GLfloat z)
{
  angle *= M_PI/180.;
  float c = cos(angle), s = sin(angle), n = sqrt(x*x + y*y + z*z);
  if (n > 0.)
    x /= n, y /= n, z /= n;
  glMultMatrixf ((GLfloat []){
      x*x*(1. - c) + c,    y*x*(1. - c) + z*s,  x*z*(1. - c) - y*s,   0,
      x*y*(1. - c) - z*s,  y*y*(1. - c) + c,    y*z*(1. - c) + x*s,   0,
      x*z*(1. - c) + y*s,  y*z*(1. - c) - x*s,  z*z*(1. - c) + c,     0,
      0,                   0,                   0,                    1
  });  
}

void glMatrixMode (GLenum mode)
{
  switch (mode) {
  case GL_MODELVIEW:  current = modelview;  break;
  case GL_PROJECTION: current = projection; break;
  default: assert (0);
  }
}

typedef struct {
  void * stack;
  int len;
} Stack;

static Stack modelview_stack = {0}, projection_stack = {0};

#define SIZE (16*sizeof (real))

void glPopMatrix (void)
{
  Stack * s = (current == modelview) ? &modelview_stack : &projection_stack;
  assert (s->len >= SIZE);
  memcpy (current, s->stack + s->len - SIZE, SIZE);
  s->len -= SIZE;
  if (s->len == 0)
    free (s->stack), s->stack = NULL;
}

void glPushMatrix (void)
{
  Stack * s = (current == modelview) ? &modelview_stack : &projection_stack;
  s->stack = realloc (s->stack, s->len + SIZE);
  memcpy (s->stack + s->len, current, SIZE);
  s->len += SIZE;
}

void glGetIntegerv (GLenum pname, GLint * params)
{
  switch (pname) {

  case GL_VIEWPORT:
    assert (TinyFramebuffer);
    params[0] = 0;
    params[1] = 0;
    params[2] = TinyFramebuffer->width;
    params[3] = TinyFramebuffer->height;
    break;
    
  default: assert (not_implemented);
  }
}

/**
## Drawing primitives */

static int Mode = -1;
static int nvertex = 0, nnormal = 0, ncolor = 0, ntexture = 0;

static inline
void reset_vertices() {
  nvertex = nnormal = ncolor = ntexture = 0;
}

void glBegin (GLenum mode) {
  assert (Mode < 0); // glBegins cannot be imbricated
  Mode = mode;
  reset_vertices();
}

static TinyColor clear = {255,255,255,0};

void glClear (GLbitfield mask)
{
  assert (TinyFramebuffer);
  if (mask & GL_COLOR_BUFFER_BIT) {
    unsigned char * p = TinyFramebuffer->image;
    for (int i = 0; i < TinyFramebuffer->width*TinyFramebuffer->height; i++, p += 4)
      for (int j = 0; j < 4; j++)
	*(p + j) = ((unsigned char *)&clear)[j];
  }
  if (mask & GL_DEPTH_BUFFER_BIT) {
    real * p = TinyFramebuffer->zbuffer;
    for (int i = 0; i < TinyFramebuffer->width*TinyFramebuffer->height; i++, p++)
      *p = 1e30;
  }
}

void glClearColor (float red, float green, float blue, float alpha) {
  //  fprintf (stderr, "%g %g %g %g\n", red, green, blue, alpha);
  clear.r = red*255;
  clear.g = green*255;
  clear.b = blue*255;
  clear.a = alpha*255;
}

static float LineWidth = 1.;

void glLineWidth (GLfloat width) {
  LineWidth = width;
}

static float PointSize = 4.;

void glPointSize (GLfloat size) {
  PointSize = size;
}

// Maximum number of vertices in any single POLYGON or TRIANGLE_FAN
#define NVERTMAX 1024
static vec3 color[NVERTMAX];
static TinyColor FgColor = {0, 0, 0, 255 };

void glColor3f (GLfloat red, GLfloat green, GLfloat blue)
{
  assert (TinyFramebuffer);
  assert (ncolor < NVERTMAX);
  color[ncolor++] = (vec3){red, green, blue};
  FgColor.r = red*255;
  FgColor.g = green*255;
  FgColor.b = blue*255;
}

static
int Face = 0;

void glColorMaterial (GLenum face, GLenum mode) {
  switch (face) {
  case GL_FRONT_AND_BACK: Face =  0; break;
  case GL_FRONT:          Face =  1; break;
  case GL_BACK:           Face = -1; break;
  default: assert (0);
  }
}

/**
## Shaders */

static real texcoord[NVERTMAX];

void glTexCoord1d (GLdouble s)
{
  assert (TinyFramebuffer);
  assert (ntexture < NVERTMAX);
  texcoord[ntexture++] = s;
}

static vec4 vertex[NVERTMAX], startv = {0,0,0,1};
static vec3 normal[NVERTMAX] = {0};
static float constant_normal_shade = 1.;

void glNormal3d (GLdouble nx, GLdouble ny, GLdouble nz)
{
  assert (TinyFramebuffer);
  assert (nnormal < NVERTMAX);
  normal[nnormal++] = vec3_normalized ((vec3){nx, ny, nz});
  constant_normal_shade = normal_shade (normal[nnormal - 1]);
}

static inline
void shaded_color (const TinyColor * color, float shade, TinyColor * frag_color)
{
  unsigned char * f = &frag_color->r;
  const unsigned char * c = &color->r;
  for (int i = 0; i < 3; i++)
    f[i] = shade*c[i];
  f[3] = c[3]; // alpha channel
}

static
int constant_normal_shader (const void * data, const vec3 bar, TinyColor * frag_color)
{
  const TinyColor * color = data;
  shaded_color (color, constant_normal_shade, frag_color);
  return 0; // the pixel is not discarded
}

static
int constant_normal_color_shader (const void * data, const vec3 bar, TinyColor * frag_color)
{
  const mat3 * nc = data;
  vec3 bc = mat3_mul (mat3_transpose (nc[0]), bar); // per-vertex color interpolation
  TinyColor color = { bc.x*255, bc.y*255, bc.z*255, 255 };
  shaded_color (&color, constant_normal_shade, frag_color);
  return 0; // the pixel is not discarded
}

static
int constant_normal_texture_shader (const void * data, const vec3 bar, TinyColor * frag_color)
{
  const vec3 * t = data;
  real bt = vec3_scalar (*t, bar); // per-vertex texture interpolation
  int i = clamp (bt, 0, 1)*(TEXTURE_WIDTH - 1);
  TinyColor color = { 255*texture[3*i], 255*texture[3*i+1], 255*texture[3*i+2], 255 };
  shaded_color (&color, constant_normal_shade, frag_color);
  return 0; // the pixel is not discarded
}

static
int vertex_normal_shader (const void * data, const vec3 bar, TinyColor * frag_color)
{
  const mat3 * nm = data;
  vec3 bn = vec3_normalized (mat3_mul (mat3_transpose (nm[0]), bar)); // per-vertex normal interpolation
  shaded_color (&FgColor, normal_shade (bn), frag_color);
  return 0; // the pixel is not discarded
}

static
int vertex_normal_color_shader (const void * data, const vec3 bar, TinyColor * frag_color)
{
  const mat3 * nm = data;
  vec3 bn = vec3_normalized (mat3_mul (mat3_transpose (nm[0]), bar)); // per-vertex normal interpolation
  vec3 bc = mat3_mul (mat3_transpose (nm[1]), bar); // per-vertex color interpolation
  TinyColor color = { bc.x*255, bc.y*255, bc.z*255, 255 };
  shaded_color (&color, normal_shade (bn), frag_color);
  return 0; // the pixel is not discarded
}

static
int vertex_normal_texture_shader (const void * data, const vec3 bar, TinyColor * frag_color)
{
  const mat3 * nm = data;
  vec3 bn = vec3_normalized (mat3_mul (mat3_transpose (nm[0]), bar)); // per-vertex normal interpolation
  real bt = vec3_scalar (nm[1].x, bar); // per-vertex texture interpolation
  int i = clamp (bt, 0, 1)*(TEXTURE_WIDTH - 1);
  TinyColor color = { 255*texture[3*i], 255*texture[3*i+1], 255*texture[3*i+2], 255 };
  shaded_color (&color, normal_shade (bn), frag_color);
  return 0; // the pixel is not discarded
}

void glVertex3d (GLdouble x, GLdouble y, GLdouble z)
{
  assert (TinyFramebuffer);
  assert (Mode >= 0);
  assert (nvertex < NVERTMAX);
  vertex[nvertex++] = mat4_mul (mat4_transpose (*((mat4 *)projection)),
				mat4_mul (mat4_transpose (*((mat4 *)modelview)),
					  (vec4){x, y, z, 1}));

  switch (Mode) {

  case GL_LINES:
    if (nvertex == 2) {
      tiny_line (vertex[0], vertex[1], &FgColor, LineWidth, TinyFramebuffer);
      reset_vertices();
    }
    break;

  case GL_LINE_STRIP: case GL_LINE_LOOP:
    if (nvertex == 1)
      startv = vertex[0];
    else {
      tiny_line (vertex[0], vertex[1], &FgColor, LineWidth, TinyFramebuffer);
      vertex[0] = vertex[1];
      nvertex = 1;
    }
    break;

  case GL_QUADS:
    if (nvertex == 4) {
      assert (nnormal == 0); // only constant shading is implemented
      assert (ntexture == 0); // textures are not implemented on quads
      tiny_triangle ((vec4[3]){vertex[0], vertex[1], vertex[3]},
		     &FgColor, constant_normal_shader, Face, // fixme: swap NULL and constant_normal_shader
		     TinyFramebuffer);
      tiny_triangle ((vec4[3]){vertex[1], vertex[2], vertex[3]},
		     &FgColor, constant_normal_shader, Face,
		     TinyFramebuffer);
      reset_vertices();
    }
    break;

  case GL_POINTS:
    tiny_point (vertex[0], &FgColor, PointSize, TinyFramebuffer);
    reset_vertices();
    break;
    
  case GL_POLYGON:
  case GL_TRIANGLE_FAN:
    break;
    
  default:
    assert (not_implemented);
  }
}

void glVertex3f (GLfloat x, GLfloat y, GLfloat z) {
  glVertex3d (x, y, z);
}

void glEnd (void)
{
  assert (TinyFramebuffer && Mode >= 0); // fixme: replace assertions
  
  switch (Mode) {

  case GL_LINE_LOOP:
    tiny_line (vertex[0], startv, &FgColor, LineWidth, TinyFramebuffer);
    break;

  case GL_TRIANGLE_FAN:
  case GL_POLYGON:
    if (nnormal == 0) {
      if (ntexture == nvertex)
	for (int i = 1; i < nvertex - 1; i++) {
	  vec3 t = { texcoord[i], texcoord[i + 1], texcoord[0] };
	  tiny_triangle ((vec4[3]){vertex[i], vertex[i + 1], vertex[0]},
			 &t, constant_normal_texture_shader, Face, TinyFramebuffer);
	}
      else if (ncolor == 0)
	for (int i = 1; i < nvertex - 1; i++)
	  tiny_triangle ((vec4[3]){vertex[i], vertex[i + 1], vertex[0]},
			 &FgColor, constant_normal_shader, Face, TinyFramebuffer);
      else if (ncolor == nvertex)
	for (int i = 1; i < nvertex - 1; i++) {
	  mat3 mc = { color[i], color[i + 1], color[0] };
	  tiny_triangle ((vec4[3]){vertex[i], vertex[i + 1], vertex[0]},
			 &mc, constant_normal_color_shader, Face, TinyFramebuffer);
	}
      else
	fprintf (stderr, "%s:%d: warning: %d != %d\n", __FILE__, __LINE__, ncolor, nvertex);
    }
    else if (nnormal == nvertex) {
      if (ntexture == nvertex)
	for (int i = 1; i < nvertex - 1; i++) {
	  mat3 nm[2] = { { normal[i], normal[i + 1], normal[0] },
			 { { texcoord[i], texcoord[i + 1], texcoord[0] } } };
	  tiny_triangle ((vec4[3]){vertex[i], vertex[i + 1], vertex[0]},
			 nm, vertex_normal_texture_shader, Face, TinyFramebuffer);
	}
      else if (ncolor == 0)
	for (int i = 1; i < nvertex - 1; i++) {
	  mat3 nm = { normal[i], normal[i + 1], normal[0] };
	  tiny_triangle ((vec4[3]){vertex[i], vertex[i + 1], vertex[0]},
			 &nm, vertex_normal_shader, Face, TinyFramebuffer);
	}
      else if (ncolor == nvertex)
	for (int i = 1; i < nvertex - 1; i++) {
	  mat3 nm[2] = { { normal[i], normal[i + 1], normal[0] },
			 { color[i],  color[i + 1],  color[0] } };
	  tiny_triangle ((vec4[3]){vertex[i], vertex[i + 1], vertex[0]},
			 nm, vertex_normal_color_shader, Face, TinyFramebuffer);
	}
      else
	fprintf (stderr, "%s:%d: warning: %d, %d, %d\n", __FILE__, __LINE__, ntexture, ncolor, nvertex);
    }
    else
      fprintf (stderr, "%s:%d: warning: %d != %d\n", __FILE__, __LINE__, nnormal, nvertex);
    break;
    
  }
  Mode = -1;
  reset_vertices();
}
