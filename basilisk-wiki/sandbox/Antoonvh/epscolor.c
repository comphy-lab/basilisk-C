/**
# A `bview` hack for vector graphics with colorbar options

![scaleable vector graphics](epscolor/nja.svg){width=500px}

Concept: place `draw_commands` in a new `save()`-like function. 
 */
#include "view.h"

bool save_PS (struct _save p); //prototype

void draw_commands (void) { //Do not call `view()`
  squares ("s", map = cool_warm, min = -.3, max = 0.3);
  box();
}

scalar s[];

int main () {
  X0 = Y0 = -L0/2;
  init_grid (16);
  foreach()
    s[] = sq(x) + y;
  view (fov = 35, psi = 0.5);
  save_PS ("nja.svg");
}
/**
The implementation:
 */
bool save_PS (struct _save p) {
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
  if (p.file && (p.fp = fopen (p.file, "w")) == NULL) {
    perror (p.file);
    return false;
  }
  if (!p.fp)
    p.fp = stdout;
  if (!strcmp (p.format, "ps") ||
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
      
      float res = view->res;
      view->res = 0.;
      view->vector=true;
    
      draw_commands();
       
      glFinish ();
      enable_fpe (FE_DIVBYZERO|FE_INVALID);
      view->active = false;
      view->vector = false;
      view->res = res;
      //draw();
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
