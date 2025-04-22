#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "utils.h"
#include "TinyPngOut.h"

/**
## Various helper functions

A helper function to write a PPM file from a RGB buffer. Downsampling
by averaging is performed if *samples* is larger than one. */

void gl_write_image (FILE * fp, const GLubyte * buffer,
		     unsigned width, unsigned height, unsigned samples)
{
  const GLubyte *ptr = buffer;

  if (samples < 1)
    samples = 1;
  if (samples > 4)
    samples = 4;

  width /= samples, height /= samples;
  fprintf (fp, "P6 %d %d 255\n", width, height);
  int x, y, j, k;
  for (y = height - 1; y >= 0; y--)
    for (x = 0; x < width; x++) {
      int r = 0, g = 0, b = 0;
      for (j = 0; j < samples; j++)
	for (k = 0; k < samples; k++) {
	  int i = (((y*samples + j)*width + x)*samples + k)*4;
	  if (ptr)
	    r += ptr[i], g += ptr[i+1], b += ptr[i+2];
	}
      fputc (r/samples/samples, fp); /* write red */
      fputc (g/samples/samples, fp); /* write green */
      fputc (b/samples/samples, fp); /* write blue */
    }
}

/**
A helper function to write a PNG file from a RGB buffer. Downsampling
by averaging is performed if *samples* is larger than one. */

void gl_write_image_png (FILE * fp, const GLubyte * buffer,
			 unsigned width, unsigned height, unsigned samples)
{
  const GLubyte *ptr = buffer;

  if (samples < 1)
    samples = 1;
  if (samples > 4)
    samples = 4;

  width /= samples, height /= samples;
  struct TinyPngOut pngout;
  enum TinyPngOut_Status status = TinyPngOut_init (&pngout, width, height, fp);
  if (status != TINYPNGOUT_OK) {
    fprintf (stderr, "error: TinyPngOut init failed\n");
    return;
  }

  int x, y, j, k;
  for (y = height - 1; y >= 0; y--)
    for (x = 0; x < width; x++) {
      int r = 0, g = 0, b = 0;
      for (j = 0; j < samples; j++)
	for (k = 0; k < samples; k++) {
	  int i = (((y*samples + j)*width + x)*samples + k)*4;
	  if (ptr)
	    r += ptr[i], g += ptr[i+1], b += ptr[i+2];
	}
      uint8_t pixel[3] = { r/samples/samples, g/samples/samples, b/samples/samples };
      status = TinyPngOut_write (&pngout, pixel, 1);
      if (status != TINYPNGOUT_OK) {
	fprintf (stderr, "error: TinyPngOut write failed\n");
	return;
      }
    }
}

/**
This is the basic OpenGL setup. */

void init_gl() {
  GLfloat light0_pos[4]   = { 0.0, 0.0, 50.0, 0.0 };

  glDisable (GL_CULL_FACE);
  glEnable (GL_DEPTH_TEST);
  glEnable (GL_NORMALIZE);

  /* speedups */
  glEnable (GL_DITHER);
  glShadeModel (GL_SMOOTH);
  glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_FASTEST);
  glHint (GL_POLYGON_SMOOTH_HINT, GL_FASTEST);

  /* light */
  glLightModeli (GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  glLightfv (GL_LIGHT0, GL_POSITION, light0_pos);
  GLfloat diffuse[4] = { 0.8, 0.8, 0.8, 1 };
  glLightfv (GL_LIGHT0, GL_DIFFUSE, diffuse);
  glEnable (GL_LIGHT0);
  glEnable (GL_LIGHTING);

  glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable (GL_COLOR_MATERIAL);
}

void gl_draw_texture (GLuint id, int width, int height)
{
  glMatrixMode (GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho (0.0, width, 0.0, height, -1.0, 1.0);
  glMatrixMode (GL_MODELVIEW);
  glPushMatrix();

  glLoadIdentity();
  glDisable (GL_LIGHTING);

  glColor3f (1,1,1);
  glEnable (GL_TEXTURE_2D);
  glBindTexture (GL_TEXTURE_2D, id);

  glBegin(GL_QUADS);
  glTexCoord2f(0, 0); glVertex3f (0, 0, 0);
  glTexCoord2f(0, 1); glVertex3f (0, 100, 0);
  glTexCoord2f(1, 1); glVertex3f (100, 100, 0);
  glTexCoord2f(1, 0); glVertex3f (100, 0, 0);
  glEnd();

  glDisable(GL_TEXTURE_2D);
  glPopMatrix();

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  glMatrixMode(GL_MODELVIEW);
}

#define RC(r,c) m[(r)+(c)*4]
#define RCM(m,r,c) (m)[(r)+(c)*4]

void matrix_multiply (float * m, const float * n)
{
  float o[16];
  int i;
  for (i = 0; i < 16; i++) o[i] = m[i];
  RC(0,0)=RCM(o,0,0)*RCM(n,0,0)+RCM(o,0,1)*RCM(n,1,0)+
          RCM(o,0,2)*RCM(n,2,0)+RCM(o,0,3)*RCM(n,3,0);
  RC(0,1)=RCM(o,0,0)*RCM(n,0,1)+RCM(o,0,1)*RCM(n,1,1)+
          RCM(o,0,2)*RCM(n,2,1)+RCM(o,0,3)*RCM(n,3,1);
  RC(0,2)=RCM(o,0,0)*RCM(n,0,2)+RCM(o,0,1)*RCM(n,1,2)+
          RCM(o,0,2)*RCM(n,2,2)+RCM(o,0,3)*RCM(n,3,2);
  RC(0,3)=RCM(o,0,0)*RCM(n,0,3)+RCM(o,0,1)*RCM(n,1,3)+
          RCM(o,0,2)*RCM(n,2,3)+RCM(o,0,3)*RCM(n,3,3);
  RC(1,0)=RCM(o,1,0)*RCM(n,0,0)+RCM(o,1,1)*RCM(n,1,0)+
          RCM(o,1,2)*RCM(n,2,0)+RCM(o,1,3)*RCM(n,3,0);
  RC(1,1)=RCM(o,1,0)*RCM(n,0,1)+RCM(o,1,1)*RCM(n,1,1)+
          RCM(o,1,2)*RCM(n,2,1)+RCM(o,1,3)*RCM(n,3,1);
  RC(1,2)=RCM(o,1,0)*RCM(n,0,2)+RCM(o,1,1)*RCM(n,1,2)+
          RCM(o,1,2)*RCM(n,2,2)+RCM(o,1,3)*RCM(n,3,2);
  RC(1,3)=RCM(o,1,0)*RCM(n,0,3)+RCM(o,1,1)*RCM(n,1,3)+
          RCM(o,1,2)*RCM(n,2,3)+RCM(o,1,3)*RCM(n,3,3);
  RC(2,0)=RCM(o,2,0)*RCM(n,0,0)+RCM(o,2,1)*RCM(n,1,0)+
          RCM(o,2,2)*RCM(n,2,0)+RCM(o,2,3)*RCM(n,3,0);
  RC(2,1)=RCM(o,2,0)*RCM(n,0,1)+RCM(o,2,1)*RCM(n,1,1)+
          RCM(o,2,2)*RCM(n,2,1)+RCM(o,2,3)*RCM(n,3,1);
  RC(2,2)=RCM(o,2,0)*RCM(n,0,2)+RCM(o,2,1)*RCM(n,1,2)+
          RCM(o,2,2)*RCM(n,2,2)+RCM(o,2,3)*RCM(n,3,2);
  RC(2,3)=RCM(o,2,0)*RCM(n,0,3)+RCM(o,2,1)*RCM(n,1,3)+
          RCM(o,2,2)*RCM(n,2,3)+RCM(o,2,3)*RCM(n,3,3);
  RC(3,0)=RCM(o,3,0)*RCM(n,0,0)+RCM(o,3,1)*RCM(n,1,0)+
          RCM(o,3,2)*RCM(n,2,0)+RCM(o,3,3)*RCM(n,3,0);
  RC(3,1)=RCM(o,3,0)*RCM(n,0,1)+RCM(o,3,1)*RCM(n,1,1)+
          RCM(o,3,2)*RCM(n,2,1)+RCM(o,3,3)*RCM(n,3,1);
  RC(3,2)=RCM(o,3,0)*RCM(n,0,2)+RCM(o,3,1)*RCM(n,1,2)+
          RCM(o,3,2)*RCM(n,2,2)+RCM(o,3,3)*RCM(n,3,2);
  RC(3,3)=RCM(o,3,0)*RCM(n,0,3)+RCM(o,3,1)*RCM(n,1,3)+
          RCM(o,3,2)*RCM(n,2,3)+RCM(o,3,3)*RCM(n,3,3);
}

void vector_multiply (float * v, const float * m)
{
  float o[4];
  int i;
  for (i = 0; i < 4; i++) o[i] = v[i];  
  v[0]=RC(0,0)*o[0]+RC(0,1)*o[1]+RC(0,2)*o[2]+RC(0,3)*o[3];
  v[1]=RC(1,0)*o[0]+RC(1,1)*o[1]+RC(1,2)*o[2]+RC(1,3)*o[3];
  v[2]=RC(2,0)*o[0]+RC(2,1)*o[1]+RC(2,2)*o[2]+RC(2,3)*o[3];
  v[3]=RC(3,0)*o[0]+RC(3,1)*o[1]+RC(3,2)*o[2]+RC(3,3)*o[3];
}

void gl_check_error()
{
  switch (glGetError()) {
  case GL_NO_ERROR: return;
  case GL_INVALID_ENUM: fprintf (stderr, "OpenGL: invalid enum\n"); break;
  case GL_INVALID_VALUE: fprintf (stderr, "OpenGL: invalid value\n"); break;
  case GL_INVALID_OPERATION: fprintf (stderr, "OpenGL: invalid operation\n");
    break;
  case GL_INVALID_FRAMEBUFFER_OPERATION:
    fprintf (stderr, "OpenGL: invalid framebuffer operation\n"); break;
  case GL_OUT_OF_MEMORY:
    fprintf (stderr, "OpenGL: out of memory\n"); break;
  case GL_STACK_UNDERFLOW:
    fprintf (stderr, "OpenGL: stack underflow\n"); break;
  case GL_STACK_OVERFLOW:
    fprintf (stderr, "OpenGL: stack overflow\n"); break;
  }
  abort();
}

void gl_get_frustum (Frustum * f)
{
  GLint v[4];
  glGetIntegerv (GL_VIEWPORT, v);
  gl_check_error();
  f->width = v[2];
  glGetFloatv (GL_MODELVIEW_MATRIX, f->m);
  gl_check_error();
  glGetFloatv (GL_PROJECTION_MATRIX, f->p);
  gl_check_error();
  float p[16];
  int i;
  for (i = 0; i < 16; i++) p[i] = f->p[i];
  matrix_multiply (p, f->m);

  /* right */
  f->n[0][0] = p[3] - p[0];
  f->n[0][1] = p[7] - p[4];
  f->n[0][2] = p[11] - p[8];
  f->d[0]    = p[15] - p[12];
   
  /* left */
  f->n[1][0] = p[3] + p[0];
  f->n[1][1] = p[7] + p[4];
  f->n[1][2] = p[11] + p[8];
  f->d[1]    = p[15] + p[12];
  
  /* top */
  f->n[2][0] = p[3] - p[1];
  f->n[2][1] = p[7] - p[5];
  f->n[2][2] = p[11] - p[9];
  f->d[2]    = p[15] - p[13];

  /* bottom */
  f->n[3][0] = p[3] + p[1];
  f->n[3][1] = p[7] + p[5];
  f->n[3][2] = p[11] + p[9];
  f->d[3]    = p[15] + p[13];
  
  /* front */
  f->n[4][0] = p[3] + p[2];
  f->n[4][1] = p[7] + p[6];
  f->n[4][2] = p[11] + p[10];
  f->d[4]    = p[15] + p[14];
  
  /* back */
  f->n[5][0] = p[3] - p[2];
  f->n[5][1] = p[7] - p[6];
  f->n[5][2] = p[11] - p[10];
  f->d[5]    = p[15] - p[14];
  
  for (i = 0; i < 6; i++) {
    float n = sqrt(f->n[i][0]*f->n[i][0] +
		   f->n[i][1]*f->n[i][1] +
		   f->n[i][2]*f->n[i][2]);
    if (n > 0.) {
      f->n[i][0] /= n; f->n[i][1] /= n; f->n[i][2] /= n;
      f->d[i] /= n;
    }
  }
}

/*
  Returns 0 if the sphere is outside the view frustum, 1, if it is
  inside and -1 if it is partly inside. 
*/
int sphere_in_frustum (double x, double y, double z, double r, Frustum * f)
{
  int I1 = 0, i, I = 1;
  for (i = 0; i < 6; i++) {
    double d = f->n[i][0]*x + f->n[i][1]*y + f->n[i][2]*z + f->d[i];
    if (d < -r) {
      I = 0;
      break;
    }
    if (d < r)
      I = -1;
  }
  if (I == 1)
    return 1;
  if (I == -1)
    I1 = -1;
  return I1;
}

/*
  Returns the diameter (in pixels) of a sphere projected on the
  screen.
*/
float sphere_diameter (double x, double y, double z, double r, Frustum * f)
{
  float v[4];
  v[0] = x; v[1] = y; v[2] = z; v[3] = 1.;
  vector_multiply (v, f->m);
  v[0] = r;
  vector_multiply (v, f->p);
  float rp = v[3] == 0. ? 0 : v[0]*f->width/v[3];
  return rp;
}

/*
  Replacement for gluPerspective in GLU from:
  https://www.khronos.org/opengl/wiki/GluPerspective_code
*/
static
void glhFrustumf2(double *matrix, double left, double right,
		  double bottom, double top,
                  double znear, double zfar)
{
  double temp, temp2, temp3, temp4;
  temp = 2.0 * znear;
  temp2 = right - left;
  temp3 = top - bottom;
  temp4 = zfar - znear;
  matrix[0] = temp / temp2;
  matrix[1] = 0.0;
  matrix[2] = 0.0;
  matrix[3] = 0.0;
  matrix[4] = 0.0;
  matrix[5] = temp / temp3;
  matrix[6] = 0.0;
  matrix[7] = 0.0;
  matrix[8] = (right + left) / temp2;
  matrix[9] = (top + bottom) / temp3;
  matrix[10] = (-zfar - znear) / temp4;
  matrix[11] = -1.0;
  matrix[12] = 0.0;
  matrix[13] = 0.0;
  matrix[14] = (-temp * zfar) / temp4;
  matrix[15] = 0.0;
}

void gl_perspective (double fovy, double aspect, double znear, double zfar)
{
  double matrix[16];
  double ymax, xmax;
  ymax = znear * tanf(fovy*M_PI/360.0);
  // ymin = -ymax;
  // xmin = -ymax * aspectRatio;
  xmax = ymax * aspect;
  glhFrustumf2 (matrix, -xmax, xmax, -ymax, ymax, znear, zfar);
  glMatrixMode (GL_PROJECTION);
  glLoadMatrixd (matrix);
}

// derived from: Mesa-7.8.2/src/glu/sgi/libutil/project.c:234
int gl_project (float objx, float objy, float objz, 
		const float modelMatrix[16], 
		const float projMatrix[16],
		const int viewport[4],
		float *winx, float *winy, float *winz)
{
  float in[4] = { objx, objy, objz, 1. };

  vector_multiply (in, modelMatrix);
  vector_multiply (in, projMatrix);
  
  if (in[3] == 0.0) return 0;

  *winx = viewport[0] + viewport[2]*(in[0]/in[3] + 1.)/2.;
  *winy = viewport[1] + viewport[3]*(in[1]/in[3] + 1.)/2.;
  *winz = (in[2]/in[3] + 1.)/2.;
  return 1;
}
