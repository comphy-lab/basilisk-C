/**
# Basilisk View

This module defines functions which compute various graphical
representations ("drawings") of Basilisk fields, including
reconstructed Volume-Of-Fluid facets and colorscale representations of
cross-sections of scalar fields. These representations are rendered
using [OpenGL](https://en.wikipedia.org/wiki/OpenGL) and can be saved
in various formats (PPM, Gnuplot, OBJ, KML, PDF, SVG etc.).

## Installation

In contrast with other Basilisk modules, this module relies on
additional libraries which needs to be installed and linked with the
Basilisk program. See [INSTALL#visualisation]() for instructions.

## Usage

A simple example would look like:

~~~literatec
...
#include "view.h"
...
event image (t = 1) {
  clear();
  draw_vof ("f");
  box();
  save ("image.ppm");
}
~~~

The *clear()* function resets the image, *draw_vof()* and *box()* add
two graphical representations and *save()* saves the resulting image
in [PPM](https://en.wikipedia.org/wiki/Netpbm_format) format. See
[User functions](view.h#user-functions) for a detailed documentation.

The resulting program needs to be linked with the appropriate
libraries. This can be done automatically using something like:

~~~bash
qcc -autolink -Wall -O2 program.c -o program -lm
~~~

or manually using e.g.:

~~~bash
qcc -Wall -O2 program.c -o program -L$BASILISK/gl -lglutils -lfb_tiny -lm
~~~

# Implementation

We include the various helper functions defined either by the system
or by the Basilisk libraries in gl/. */

#include <gl/framebuffer.h>
#include <gl/trackball.h>
#include <gl/utils.h>
#pragma autolink -L$BASILISK/gl -lglutils $OPENGLIBS

#include "utils.h"
#include "input.h"

/**
## A cache of "compiled" expressions

The cache has a maximum size and least-used expressions are discarded
first. */

typedef struct {
  char * expr;
  scalar s;
} cexpr;

static scalar get_cexpr (cexpr * cache, const char * expr)
{
  cexpr * c = cache;
  while (c->expr) {
    if (!strcmp (c->expr, expr)) {
      // move this expression to the top of the cache.
      // the "top" is the last element.
      cexpr tmp = *c;
      while ((c + 1)->expr)
	*c = *(c + 1), c++;
      *c = tmp;
      return c->s;
    }
    c++;
  }
  return (scalar){-1};
}

static cexpr * add_cexpr (cexpr * cache, int maxlen,
			  const char * expr, scalar s)
{
  cexpr * c = cache;
  while (c->expr) c++;
  int len = c - cache;
  if (len < maxlen) {
    cache = realloc (cache, sizeof(cexpr)*(len + 2));
    c = &cache[len];
  }
  else {
    // discard first expression
    c = cache;
    free (c->expr);
    scalar s = c->s;
    delete ({s});
    // shift remaining expressions
    while ((c + 1)->expr)
      *c = *(c + 1), c++;
  }
  c->expr = strdup (expr);
  c->s = s;
  (c + 1)->expr = NULL;
  return cache;
}

static void free_cexpr (cexpr * cache)
{
  cexpr * c = cache;
  while (c->expr) {
    free (c->expr);
    scalar s = c->s;
    delete ({s});
    c++;
  }
  free (cache);
}

/**
## The *bview* class

Contains the definition of the current view. */

typedef void (* MapFunc) (coord *);

struct _bview {
  float tx, ty, sx, sy, sz;
  float quat[4];
  float fov;
  float tz, near, far;

  bool gfsview;   // rotate axis to match gfsview
  bool reversed;  // reverse normals
  
  float bg[3];
  float lc;
  float res;

  unsigned width, height, samples;

  framebuffer * fb;
  Frustum frustum; // the view frustum

  MapFunc map; // an optional mapping function
  
  int ni; // number of items drawn
  
  bool active;

  cexpr * cache; // a cache of compiled expressions
  int maxlen; // the maximum number of cached expressions
};

typedef struct _bview bview;

/**
The allocator method. */

bview * bview_new()
{
  bview * p = qcalloc (1, bview);

  p->tx = p->ty = 0;
  p->sx = p->sy = p->sz = 1.;
  p->quat[0] = p->quat[1] = p->quat[2] = 0; p->quat[3] = 1;
  p->fov = 24.;
  gl_trackball (p->quat, 0.0, 0.0, 0.0, 0.0);

#if dimension <= 2
  p->bg[0] = 1; p->bg[1] = 1; p->bg[2] = 1;
#else
  p->bg[0] = 0.3; p->bg[1] = 0.4; p->bg[2] = 0.6;
#endif
  p->res = 1.;
  p->lc = 0.004;

  p->samples = 4;
  p->width = 600*p->samples, p->height = 600*p->samples;

  /* OpenGL somehow generates floating-point exceptions... turn them off */
  disable_fpe (FE_DIVBYZERO|FE_INVALID);

  p->fb = framebuffer_new (p->width, p->height);

  init_gl();
  p->active = false;
  
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
  
  return p;
}

/**
The destructor method. */

void bview_destroy (bview * p)
{
  framebuffer_destroy (p->fb);
  if (p->cache)
    free_cexpr (p->cache);
  free (p);
}

/**
For the moment there is a single (static) current view. */

static bview * _view = NULL;

/**
The current view needs to be destroyed when we exit Basilisk. This is
done by adding this callback to the free_solver() lists of
destructors. */

static void destroy_view()
{
  assert (_view);
  bview_destroy (_view);
}

bview * get_view() {
  if (!_view) {
    _view = bview_new();
    free_solver_func_add (destroy_view);
  }
  return _view;
}

/**
The main drawing function. */

static void redraw() {
  bview * view = get_view();
    
  /* OpenGL somehow generates floating-point exceptions... turn them off */
  disable_fpe (FE_DIVBYZERO|FE_INVALID);

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();

  if (view->far <= view->near) { // "traditional" camera parameters
    double max = 2.;    
    gl_perspective (view->fov, view->width/(float)view->height, 1., 1. + 2.*max);

    glMatrixMode (GL_MODELVIEW);	    
    glLoadIdentity ();
    glTranslatef (view->tx, view->ty, - (1. + max));
  }
  else { // camera parameters compatible with interactive Basilisk View
    gl_perspective (view->fov, view->width/(float)view->height,
		    view->near, view->far);

    glMatrixMode (GL_MODELVIEW);	    
    glLoadIdentity ();
    glTranslatef (view->tx, view->ty, view->tz);
  }
    
  GLfloat m[4][4];
  gl_build_rotmatrix (m, view->quat);
  glMultMatrixf (&m[0][0]);
  
  if (view->gfsview) { // rotate to match gfsview parameters
    m[0][0] = 0., m[0][1] =  0., m[0][2] = -1.;
    m[1][0] = 0., m[1][1] = -1., m[1][2] =  0.;
    m[2][0] = 1., m[2][1] =  0., m[2][2] =  0.;
    glMultMatrixf (&m[0][0]);
  }
  
  glScalef (view->sx/L0, view->sy/L0, view->sz/L0);

  glClearColor (view->bg[0], view->bg[1], view->bg[2], 0.);
  glClear (GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

  gl_get_frustum (&view->frustum);
  
  view->active = true;
  view->ni = 0;
}

/**
This is called by graphics primitives before drawing. */

bview * draw() {
  bview * view = get_view();
  if (!view->active)
    redraw();
  else
    /* OpenGL somehow generates floating-point exceptions... turn them off
       See also: https://bugs.freedesktop.org/show_bug.cgi?id=108856 */
    disable_fpe (FE_DIVBYZERO|FE_INVALID);
  glMatrixMode (GL_PROJECTION);
  glTranslatef (0, 0, - 1e-4); // next object is drawn below the current one
  return view;
}

/**
## Helper function for parallel image composition

compose_image() returns an image buffer made by composition of the
framebuffer images on each of the MPI processes. */

typedef void * pointer; // fixme: trace is confused by pointers

#if !_MPI
trace
static pointer compose_image (bview * view) {
  return framebuffer_image((view)->fb);
}
#else // _MPI
#if dimension <= 2
typedef struct {
  GLubyte a[4];
} RGBA;

static void compose_image_op (void * pin, void * pout, int * len,
			      MPI_Datatype * dptr)
{
  RGBA * rin = pin, * out = pout;
  for (int i = 0; i < *len; i++,rin++,out++)
    if (out->a[3] == 0)
      *out = *rin;
}

trace
static pointer compose_image (bview * view)
{
  unsigned char * image = framebuffer_image (view->fb);
  assert (image);
  if (npe() > 1) {
    MPI_Op op;
    MPI_Op_create (compose_image_op, true, &op);    
    MPI_Datatype rgba;
    MPI_Type_contiguous (4, MPI_BYTE, &rgba);
    MPI_Type_commit (&rgba);
    int size = view->width*view->height;
    if (pid() == 0)
      MPI_Reduce (MPI_IN_PLACE, image, size, rgba, op, 0, MPI_COMM_WORLD);
    else
      MPI_Reduce (image, image, size, rgba, op, 0, MPI_COMM_WORLD);
    MPI_Op_free (&op);
    MPI_Type_free (&rgba);
  }
  return image;
}
#else /* 3D */
typedef struct {
  GLubyte a[4];
  float depth;
} RGBA;

static void compose_image_op (void * pin, void * pout, int * len,
			      MPI_Datatype * dptr)
{
  RGBA * rin = pin, * out = pout;
  for (int i = 0; i < *len; i++,rin++,out++)
    if (out->depth > rin->depth)
      *out = *rin;
}

trace
static pointer compose_image (bview * view)
{
  unsigned char * image = framebuffer_image (view->fb);
  assert (image);
  if (npe() > 1) {
    MPI_Op op;
    MPI_Op_create (compose_image_op, true, &op);
    MPI_Datatype rgba;
    MPI_Type_create_struct (2,
			    (int[]){4,1},
			    (MPI_Aint[]){0,4},
			    (MPI_Datatype[]){MPI_BYTE, MPI_FLOAT},
			    &rgba);
    MPI_Type_commit (&rgba);
    float * depth = framebuffer_depth (view->fb);
    int size = view->width*view->height;
    RGBA * buf = malloc (size*sizeof(RGBA));
    unsigned char * ptr = image;
    float * dptr = depth;
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < 4; j++)
	buf[i].a[j] = *ptr++;
      buf[i].depth = *dptr++;
    }
    if (pid() == 0) {
      MPI_Reduce (MPI_IN_PLACE, buf, size, rgba, op, 0, MPI_COMM_WORLD);
      unsigned char * ptr = image;
      for (int i = 0; i < size; i++)
	for (int j = 0; j < 4; j++)
	  *ptr++ = buf[i].a[j];
    }
    else
      MPI_Reduce (buf, buf, size, rgba, op, 0, MPI_COMM_WORLD);
    free (buf);
    MPI_Op_free (&op);
    MPI_Type_free (&rgba);
  }
  return image;
}
#endif /* 3D */
#endif /* _MPI */

#include "vertexbuffer.h"

/**
# User functions

Drawing user functions are defined in [draw.h](). */

#include "draw.h"

/**
## *load()*: read drawing commands from a file or buffer

The commands are calls of user functions. They can be read from a
file defined by *fp* or *file*, or from the memory buffer *buf*.

Besides the *load()*, *save()* and drawing functions defined in
[draw.h](), valid drawing commands also include:

* [*restore()*](output.h#dump-basilisk-snapshots)
* [*dump()*](output.h#dump-basilisk-snapshots)
* [*input_gfs()*](input.h#input_gfs-gerris-simulation-format)
*/

bool load (FILE * fp = NULL, char * file = NULL, Array * buf = NULL);

static void bview_draw (bview * view)
{
  if (!view->active)
    return;
  view->active = false;
  glFinish ();
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}

/**
## *save()*: saves drawing to a file in various formats

The file to write to is given either using its name *file* or the file
pointer *fp*. If neither is specified, the default is *stdout*.

If *file* is used, options for 'convert' or 'ffmpeg' can be given in
*opt*.

The format to use is given either explicitly using *format*, or, if a
*file* name is given, using the file name extension (i.e. *.ppm*,
*.bv*, etc.). If neither is specified, the default is "ppm".

For vector graphics, the base line width can be specified using
*lw*. The default is one.

The recognised file formats are:

* "ppm": [Portable PixMap](https://en.wikipedia.org/wiki/Netpbm_format) 
         format. A basic uncompressed image format.
* "png", "jpg": Compressed image formats. Will only work if the
                *convert* command from
                [ImageMagick](http://imagemagick.org) is installed on
                the system.
* "mp4", "gif", "ogv": Compressed animation formats. Will only work if
                       [ffmpeg](https://www.ffmpeg.org) is installed on 
                       the system.

The following formats are no longer supported, but this could be fixed
later:

* "bv": Basilisk View format. Saves all Basilisk function calls necessary to 
        reproduce the figure. Use together with 
        [load()](view.h#load-read-drawing-commands-from-a-file-or-buffer).
* "gnu": Gnuplot format. Saves a vector graphics (3D) representation 
         of the objects.
* "obj": [Wavefront 3D object format](https://en.wikipedia.org/wiki/Wavefront_.obj_file). Can be read by a number of 3D visualisation tools.
* "kml": [Keyhole Markup Language](https://en.wikipedia.org/wiki/Keyhole_Markup_Language). Can be used with Google Earth and other [GIS](https://en.wikipedia.org/wiki/Geographic_information_system).
* "ps", "eps", "tex", "pdf", "svg", "pgf": the various 
         [vector graphics](https://en.wikipedia.org/wiki/Vector_graphics) 
         formats supported by [gl2ps](http://www.geuz.org/gl2ps).

Note that MPI-parallel output is only implemented for the "ppm" format
at the moment. Other animation and image formats will be automatically
converted to PPM when using MPI. */

trace
bool save (char * file = NULL, char * format = "ppm", char * opt = NULL,
	   FILE * fp = NULL,
	   float lw = 0,
	   int sort = 0, int options = 0,
	   
	   bview * view = NULL)
{
  if (file) {
    char * s = strchr (file, '.'), * dot = s;
    while (s) {
      dot = s;
      s = strchr (s + 1, '.');
    }
    if (dot)
      format = dot + 1;
  }

  if (!view)
    view = get_view();

  if ((!strcmp (format, "png") && which ("convert")) ||
      !strcmp (format, "jpg") ||
      (file && is_animation (file))) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0) {
      FILE * fp = open_image (file, opt);
      if (!fp) {
	perror (file);
	return false;
      }      
      gl_write_image (fp, image, view->width, view->height, view->samples);
      close_image (file, fp);
    }
    return true;
  }  
  
  if (file && (fp = fopen (file, "w")) == NULL) {
    perror (file);
    return false;
  }
  if (!fp)
    fp = stdout;

  if (!strcmp (format, "ppm")) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0)
      gl_write_image (fp, image, view->width, view->height, view->samples);    
  }

  else if (!strcmp (format, "png")) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0)
      gl_write_image_png (fp, image, view->width, view->height, view->samples);
  }

  else if (!strcmp (format, "bv")) {
#if 1 // fixme: not implemented yet
    fprintf (ferr, "save(): error: the '%s' format is no longer supported\n",
	     format);
#else
    assert (history);
    fprintf (fp,
	     "view (fov = %g, quat = {%g,%g,%g,%g}, "
	     "tx = %g, ty = %g, "
	     "bg = {%g,%g,%g}, "
	     "width = %d, height = %d, samples = %d"
	     ");\n",
	     view->fov,
	     view->quat[0], view->quat[1], view->quat[2], view->quat[3],
	     view->tx, view->ty,
	     view->bg[0], view->bg[1], view->bg[2],
	     view->width/view->samples, view->height/view->samples,
	     view->samples);
    fwrite (history->p, 1, history->len, fp);
#endif
  }
  
  else if (!strcmp (format, "gnu") ||
	   !strcmp (format, "obj") ||
	   !strcmp (format, "kml") ||
	   !strcmp (format, "ps")  ||
	   !strcmp (format, "eps") ||
	   !strcmp (format, "tex") ||
	   !strcmp (format, "pdf") ||
	   !strcmp (format, "svg") ||
	   !strcmp (format, "pgf"))
    fprintf (ferr, "save(): error: the '%s' format is no longer supported\n",
	     format);

  else {
    fprintf (ferr, "save(): unknown format '%s'\n", format);
    if (file) {
      fclose (fp);
      remove (file);
    }
    return false;
  }

  fflush (fp);
  if (file)
    fclose (fp);

  return true;
}

/**
## Implementation of the *load()* function.

The functions below parse a text file and perform the corresponding
function calls. */

static char * remove_blanks (char * line)
{
  while (strchr (" \t", *line)) line++;
  char * s = line, * cur = line;
  bool instring = false;
  while (*s != '\0' && *s != '#') {
    if (*s == '"')
      instring = !instring;
    if (instring || !strchr (" \t", *s))
      *cur++ = *s;
    s++;
  }
  *cur = '\0';
  return line;
}

#include "parse.h"

bool process_line (char * line)
{
  if (line[0] == '\0')
    return true;
  char * s = mystrtok (remove_blanks (line), "(");
  if (!s)
    return true;

  if (!strcmp (s, "restore")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file) {
      bview * view = get_view();
      if (view->cache) {
	free_cexpr (view->cache);
	view->cache = calloc (1, sizeof (cexpr));
      }
      if (!restore (file = file, list = all))
	fprintf (ferr, "could not restore from '%s'\n", file);
      else {
	restriction (all);
	fields_stats();
	clear();
      }
    }
  }

  else if (!strcmp (s, "dump")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    dump (file = file);
  }
  
  else if (!strcmp (s, "input_gfs")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file) {
      input_gfs (file = file, list = all);
      restriction (all);
      fields_stats();
      clear();
    }
  }

  else if (!strcmp (s, "save")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file)
      save (file = file);
  }

  else if (!strcmp (s, "load")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file)
      load (file = file);
  }

  /**
  The [draw_get.h]() file is generated automatically by [params.awk]()
  and contains parsing commands for the functions defined in
  [draw.h](). */

  #include "draw_get.h"

  else if (!strcmp (s, "end_mirror"))
    end_mirror();
  
  else if (!strcmp (s, "end_translate"))
    end_translate();
  
  else if (!strcmp (s, "clear"))
    clear();

  else if (s[0] != '\n' && s[0] != '\0')
    fprintf (ferr, "load(): syntax error: '%s'\n", s);

  return true;
}

bool load (FILE * fp = NULL, char * file = NULL, Array * buf = NULL)
{
  if (file) {
    fp = fopen (file, "r");
    if (!fp) {
      perror (file);
      return false;
    }
  }

  if (fp) { // read lines from file
    char line[256];
    while (fgets (line, 256, fp) && process_line (line));
  }
  else if (buf) { // read lines from buffer
    int i = 0;
    char * s = (char *) buf->p;
    while (i < buf->len) {
      char * start = s;
      while (i < buf->len && *s != '\n')
	s++, i++;
      if (*s == '\n' && ++s > start) {
	char line[s - start + 1];
	strncpy (line, start, s - start);
	line[s - start] = '\0';
	process_line (line);
      }
    }
  }
  return true;
}
