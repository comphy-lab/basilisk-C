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
Basilisk program. See [gl/INSTALL]() for instructions.

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
libraries using e.g.:

~~~bash
qcc -Wall -O2 program.c -o program \
    -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
~~~

or

~~~bash
qcc -Wall -O2 program.c -o program \
    -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11
~~~

depending on which version of OpenGL should be used.

# Implementation

We include the various helper functions defined either by the system
or by the Basilisk libraries in gl/. */

#if defined(__APPLE__)
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif

#include <gl/framebuffer.h>
#include <gl/gl2ps/gl2ps.h>
#include <gl/trackball.h>
#include <gl/utils.h>
#pragma autolink -L$BASILISK/gl -lglutils $OPENGLIBS

#include "utils.h"
#include "input.h"

/**
## The *bview* class

Contains the definition of the current view. */

struct _bview {
  float tx, ty, sx, sy, sz;
  float quat[4];
  float fov;

  bool gfsview;   // rotate axis to match gfsview
  bool vector;    // are we doing vector graphics?
  bool reversed;  // reverse normals
  
  float bg[3];
  float lc;
  float res;

  unsigned width, height, samples;

  framebuffer * fb;
  Frustum frustum; // the view frustum

  void (* map) (coord *); // an optional mapping function
  
  int ni; // number of items drawn

  bool active;
};

typedef struct _bview bview;

/**
The allocator method. */

bview * bview_new() {
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
  p->lc = 0.001;

  p->vector = false;
  
  p->samples = 4;
  p->width = 600*p->samples, p->height = 600*p->samples;

  p->fb = framebuffer_new (p->width, p->height);
  
  init_gl();
  p->active = false;
  
  return p;
}

/**
The destructor method. */

void bview_destroy (bview * p)
{
  framebuffer_destroy (p->fb);
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
  double max = 2.;    
#if 0  
  GList * symmetries = get_symmetries (list);
  max = gfs_gl_domain_extent (domain, symmetries);
#endif
  
  gluPerspective (view->fov, view->width/(float)view->height, 1., 1. + 2.*max);
  glMatrixMode (GL_MODELVIEW);
	    
  glLoadIdentity ();
  glTranslatef (view->tx, view->ty, - (1. + max));
    
  GLfloat m[4][4];
  gl_build_rotmatrix (m, view->quat);
  glMultMatrixf (&m[0][0]);
  
  if (view->gfsview) { // rotate to match gfsview parameters
    m[0][0] = 0., m[0][1] =  0., m[0][2] =  -1.;
    m[1][0] = 0., m[1][1] = -1., m[1][2] =   0.;
    m[2][0] = 1., m[2][1] =  0., m[2][2] =   0.;
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
  RGBA * in = pin, * out = pout;
  for (int i = 0; i < *len; i++,in++,out++)
    if (out->a[3] == 0)
      *out = *in;
}

trace
static pointer compose_image (bview * view)
{
  unsigned char * image = framebuffer_image (view->fb);
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
  RGBA * in = pin, * out = pout;
  for (int i = 0; i < *len; i++,in++,out++)
    if (out->depth > in->depth)
      *out = *in;
}

trace
static pointer compose_image (bview * view)
{
  unsigned char * image = framebuffer_image (view->fb);
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
    fbdepth_t * depth = framebuffer_depth (view->fb);
    int size = view->width*view->height;
    RGBA * buf = malloc (size*sizeof(RGBA));
    unsigned char * ptr = image;
    fbdepth_t * dptr = depth;
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

/**
# User functions

Drawing user functions are defined in [draw.h](). */

#include "draw.h"

/**
## *load()*: read drawing commands from a file or buffer

The commands are calls of user functions. They can be read from a
file defined by *fp* or *file*, or from the memory buffer *buf*.

If *history* is not *NULL*, the commands are appended to this buffer. 

Besides the *load()*, *save()* and drawing functions defined in
[draw.h](), valid drawing commands also include:

* [*restore()*](output.h#dump-basilisk-snapshots)
* [*dump()*](output.h#dump-basilisk-snapshots)
* [*input_gfs()*](input.h#input_gfs-gerris-simulation-format)
* *display()*: forces redrawing.
* *show()*: writes on standard error the content of the current command history.
* *quit()*: stops parsing commands.
*/

struct _load {
  FILE * fp;   // read commands from this file
  char * file; // read commands from this file
  Array * buf; // read commands from this buffer
  Array * history; // append history to this one
};

bool load (struct _load p);

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

struct _save {
  char * file, * format, * opt;
  FILE * fp;
  float lw; /* base line width for vector drawings */
  int sort, options;

  Array * history; // command history
  bview * view;
};

static void bview_draw (bview * view)
{
  if (!view->active)
    return;
  view->active = false;
  glFinish ();
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}

static void redraw_feedback (struct _save * p)
{
  bview * view = p->view ? p->view : get_view();
  assert (p->history);
  if (p->history->len) {
    float res = view->res;
    view->res = 0.;
    view->vector = true; // vector graphics
    redraw();
    // we use a display list, just as a workaround for some buggy
    // feedback buffer implementations
    int list = glGenLists (1);
    glNewList (list, GL_COMPILE);
    load (buf = p->history);
    glEndList();
    glCallList (list);
    glFinish ();
    glDeleteLists (list, 1);
    enable_fpe (FE_DIVBYZERO|FE_INVALID);
    view->active = false;
    view->vector = false;
    view->res = res;
  }
}

#define MAXBUFFSIZE (1 << 28) // 1 GB of feedback buffer

trace
bool save (struct _save p)
{
  char ppm[] = "ppm";
  if (!p.format) {
    p.format = ppm;
    if (p.file) {
      char * s = strchr (p.file, '.'), * dot = s;
      while (s) {
	dot = s;
	s = strchr (s + 1, '.');
      }
      if (dot)
	p.format = dot + 1;
    }
  }

  bview * view = p.view ? p.view : get_view();
  
  if (!strcmp (p.format, "png") ||
      !strcmp (p.format, "jpg") ||
      (p.file && is_animation (p.file))) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0) {
      FILE * fp = open_image (p.file, p.opt);
      gl_write_image (fp, image, view->width, view->height, view->samples);
      close_image (p.file, fp);
    }
    return true;
  }  
  
  if (p.file && (p.fp = fopen (p.file, "w")) == NULL) {
    perror (p.file);
    return false;
  }
  if (!p.fp)
    p.fp = stdout;

  if (!strcmp (p.format, "ppm")) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0)
      gl_write_image (p.fp, image, view->width, view->height, view->samples);    
  }

  else if (!strcmp (p.format, "bv")) {
    assert (p.history);
    fprintf (p.fp,
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
    fwrite (p.history->p, 1, p.history->len, p.fp);
  }
  
  else if (!strcmp (p.format, "gnu") ||
	   !strcmp (p.format, "obj") ||
	   !strcmp (p.format, "kml")) {
    int format = (!strcmp (p.format, "gnu") ? FEEDBACK_GNU :
		  !strcmp (p.format, "obj") ? FEEDBACK_OBJ :
		  !strcmp (p.format, "kml") ? FEEDBACK_KML :
		  -1);
    unsigned buffsize = 1 << 24;
    bool done = false;
    while (!done && buffsize <= MAXBUFFSIZE) {
      float * f = gl_feedback_begin (buffsize);
      redraw_feedback (&p);
      done = gl_feedback_end (f, p.fp, format);
      buffsize *= 2;
    }
    if (!done)
      fprintf (ferr, "save(): error: exceeded maximum feedback buffer size\n");
  }
  
  else if (!strcmp (p.format, "ps") ||
	   !strcmp (p.format, "eps") ||
	   !strcmp (p.format, "tex") ||
	   !strcmp (p.format, "pdf") ||
	   !strcmp (p.format, "svg") ||
	   !strcmp (p.format, "pgf")) {
    GLint format = (!strcmp (p.format, "ps") ? GL2PS_PS :
		    !strcmp (p.format, "eps") ? GL2PS_EPS :
		    !strcmp (p.format, "tex") ? GL2PS_TEX :
		    !strcmp (p.format, "pdf") ? GL2PS_PDF :
		    !strcmp (p.format, "svg") ? GL2PS_SVG :
		    !strcmp (p.format, "pgf") ? GL2PS_PGF :
		    -1);
    GLint state = GL2PS_OVERFLOW;
    GLint sort = p.sort ? p.sort : GL2PS_SIMPLE_SORT;
    GLint options = p.options ? p.options : (GL2PS_SIMPLE_LINE_OFFSET |
					     GL2PS_SILENT |
					     GL2PS_BEST_ROOT |
					     GL2PS_OCCLUSION_CULL |
					     GL2PS_USE_CURRENT_VIEWPORT |
					     GL2PS_TIGHT_BOUNDING_BOX);
    unsigned buffsize = 1 << 24;
    while (state == GL2PS_OVERFLOW && buffsize <= MAXBUFFSIZE) {
      gl2psBeginPage ("", "bview",
		      NULL,
		      format, sort, options, 
		      GL_RGBA, 0, NULL, 
		      0, 0, 0,
		      buffsize, p.fp, "");
      redraw_feedback (&p);
      disable_fpe (FE_DIVBYZERO|FE_INVALID);
      state = gl2psEndPage();
      enable_fpe (FE_DIVBYZERO|FE_INVALID);
      buffsize *= 2;
    }
    if (state == GL2PS_OVERFLOW)
      fprintf (ferr, "save(): error: exceeded maximum feedback buffer size\n");
  }

  else {
    fprintf (ferr, "save(): unknown format '%s'\n", p.format);
    if (p.file) {
      fclose (p.fp);
      remove (p.file);
    }
    return false;
  }

  fflush (p.fp);
  if (p.file)
    fclose (p.fp);

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

static void fields_stats()
{
  fprintf (ferr, "# t = %g, fields = {", t);
  for (scalar s in all)
    fprintf (ferr, " %s", s.name);
  fputs (" }\n", ferr);
  fprintf (ferr, "# %12s: %12s %12s %12s %12s\n",
	   "name", "min", "avg", "stddev", "max");
  for (scalar s in all) {
    stats ss = statsf (s);
    fprintf (ferr, "# %12s: %12g %12g %12g %12g\n",
	     s.name, ss.min, ss.sum/ss.volume, ss.stddev, ss.max);
  }
}

static void draw_append (char * buf, Array * history, FILE * interactive)
{
  if (interactive) {
    if (history->len)
      load (buf = history);
    save (fp = interactive);
  }
  array_append (history, buf, strlen(buf)*sizeof(char));
}

/**
The [draw_get.h]() file is generated automatically by [params.awk]()
and contains parsing commands for the functions defined in
[draw.h](). */

#include "draw_get.h"

//#if dimension == 2
#include "front2d.h"	
//#endif

static bool process_line (char * line, Array * history, FILE * interactive)
{
  if (line[0] == '\0')
    return true;
  char * buf = strdup (line);
  char * s = strtok (remove_blanks (line), "(");
  if (!s) {
    free (buf);
    return true;
  }

  if (!strcmp (s, "restore")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file) {
      if (!restore (file = file, list = all))
	fprintf (ferr, "could not restore from '%s'\n", file);
      else {
	restriction (all);
	fields_stats();
	clear();
	// rebuild display list using history
	if (history->len && load (buf = history) && interactive)
	  save (fp = interactive);
      }
    }
  }

  else if (!strcmp (s, "dump")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    dump (file = file);
  }

  else if (!strcmp (s, "draw_front")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file) {	
    	if (!draw_front (file = file))
				fprintf (ferr, "could not restore from '%s'\n", file);
      else {
				fprintf (ferr, "Restored front from '%s'\n", file);
			}
			
		}
    draw_append (buf, history, interactive);
  }
  
  else if (!strcmp (s, "input_gfs")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file) {
      input_gfs (file = file, list = all);
      restriction (all);
      fields_stats();
      clear();
      // rebuild display list using history
      if (history->len && load (buf = history) && interactive)
	save (fp = interactive);
    }
  }

  else if (!strcmp (s, "save")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file)
      save (file = file, history = history);
  }

  else if (!strcmp (s, "load")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file && load (file = file, history = history) && interactive) {
      load (buf = history);
      save (fp = interactive);
    }
  }
        
  else if (!strcmp (s, "cells")) {
    struct _cells p = {{0}};
    _cells_get (&p);
    cells (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "vectors")) {
    struct _vectors p = {0};
    _vectors_get (&p);
    vectors (p);
    draw_append (buf, history, interactive);
  }
        
  else if (!strcmp (s, "draw_vof")) {
    struct _draw_vof p = {0};
    _draw_vof_get (&p);
    if (draw_vof (p))
      draw_append (buf, history, interactive);
  }
    
  else if (!strcmp (s, "isoline")) {
    struct _isoline p = {0};
    _isoline_get (&p);
    if (isoline (p))
      draw_append (buf, history, interactive);
  }
    
  else if (!strcmp (s, "squares")) {
    struct _squares p = {0};
    _squares_get (&p);
    squares (p);
    draw_append (buf, history, interactive);
  }
  
  else if (!strcmp (s, "begin_translate")) {
    struct _translate p = {0};
    _translate_get (&p);
    begin_translate (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "end_translate")) {
    end_translate();
    draw_append (buf, history, interactive);
  }
  
  else if (!strcmp (s, "begin_mirror")) {
    struct _mirror p = {{0}};
    _mirror_get (&p);
    begin_mirror (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "end_mirror")) {
    end_mirror();
    draw_append (buf, history, interactive);
  }
  
  else if (!strcmp (s, "squares")) {
    struct _squares p = {0};
    _squares_get (&p);
    squares (p);
    draw_append (buf, history, interactive);
  }
    
  else if (!strcmp (s, "isosurface")) {
    struct _isosurface p = {0};
    _isosurface_get (&p);
    isosurface (p);
    draw_append (buf, history, interactive);
  }
    
  else if (!strcmp (s, "draw_string")) {
    struct _draw_string p = {0};
    _draw_string_get (&p);
    draw_string (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "display")) {
    if (interactive && history->len && load (buf = history))
      save (fp = interactive);
  }

  else if (!strcmp (s, "clear")) {
    clear();
    if (interactive)
      save (fp = interactive);
    history->len = 0;
  }

  else if (!strcmp (s, "show")) {
    if (interactive && history->len)
      save (fp = ferr, format = "bv", history = history);
  }

  else if (!strcmp (s, "box")) {
    struct _box p = {0};
    _box_get (&p);
    box (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "view")) {
    struct _view_set p = {0};
    _view_set_get (&p);
    view (p);
    if (p.width || p.height || p.samples) {
      // rebuild display list using history
      if (history->len && load (buf = history) && p.samples && interactive)
	save (fp = interactive);
    }
  }

  else if (!strcmp (s, "quit")) {
    free (buf);
    return false; // quit
  }
  
  else if (s[0] != '\n')
    fprintf (ferr, "load(): syntax error: '%s'\n", s);

  free (buf);
  return true;
}

bool load (struct _load p) {
  if (p.file) {
    p.fp = fopen (p.file, "r");
    if (!p.fp) {
      perror (p.file);
      return false;
    }
  }

  Array * history = array_new();
  if (p.fp) { // read lines from file
    char line[256];
    while (fgets (line, 256, p.fp) && process_line (line, history, NULL));
  }
  else if (p.buf) { // read lines from buffer
    int i = 0;
    char * s = (char *) p.buf->p;
    while (i < p.buf->len) {
      char * start = s;
      while (i < p.buf->len && *s != '\n')
	s++, i++;
      if (*s == '\n' && ++s > start) {
	char line[s - start + 1];
	strncpy (line, start, s - start);
	line[s - start] = '\0';
	process_line (line, history, NULL);
      }
    }
  }
  if (p.history)
    array_append (p.history, history->p, history->len);
  array_free (history);

  return true;
}