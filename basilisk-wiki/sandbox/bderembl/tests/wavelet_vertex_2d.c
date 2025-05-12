/**
# Wavelet transforms and filtering

This simple example illustrates how to compute wavelet transforms and
perform filtering and mesh adaptation. We do this in one dimension,
using bi-trees. */

#include "grid/multigrid.h"
#include "bderembl/libs/inner-vertex.h"
#include "Antoonvh/my_vertex.h"

/**
   Testing several restriction functions
 */

static inline void restriction_coarsen_vert2 (Point point, scalar s) {
#if (dimension == 1)
  s[] = (fine(s,1,0,0) + 2*fine(s,0,0,0) + fine(s,-1,0,0))/4.;
#elif (dimension == 2)
  s[] = (4*fine(s,0,0,0) + 
         2*fine(s,1,0,0) + 2*fine(s,-1,0,0) +
	 2*fine(s,0,1,0) + 2*fine(s,0,-1,0) + 
         fine(s,1,1,0) + fine(s,-1,1,0) +
	 fine(s,1,-1,0) + fine(s,-1,-1,0))/16.;
#endif
}


/**
   Rewriting the multigrid wavelet and inverse wavelet functions with vertex loops
 */


void wavelet_vertex (scalar s, scalar w)
{
  restriction ({s});
  for (int l = grid->maxdepth - 1; l >= 0; l--) {
    foreach_coarse_vertex_level (l) {
      foreach_child()
        w[] = s[];
      s.prolongation (point, s);
      foreach_child() {
        double sp = s[];
        s[] = w[];
        /* difference between fine value and its prolongation */
        w[] -= sp;
      }
    }
    boundary_level ({w}, l + 1);
  }
  /* root cell */
  foreach_vertex_level(0) 
    w[] = s[];
  boundary_level ({w}, 0);
}

void inverse_wavelet_vertex (scalar s, scalar w)
{
  foreach_vertex_level(0) 
    s[] = w[];
  boundary_level ({s}, 0);
  for (int l = 0; l <= grid->maxdepth - 1; l++) {
    foreach_coarse_vertex_level (l) {
      s.prolongation (point, s);
      foreach_child()
        s[] += w[];
    }
    boundary_level ({s}, l + 1);
  }
}





/**
A simple function to output fields on each level. */

void write_level (scalar * list, const char * tag, FILE * fp)
{
  for (int l = 0; l <= depth(); l++)
    foreach_vertex_level (l, serial) {
      fprintf (fp, "%s%d %g %g ", tag, l, x, y);
      for (scalar s in list)
	fprintf (fp, "%g ", s[]);
      fputc ('\n', fp);
    }
}

void write_vertex (scalar * list, const char * tag, FILE * fp)
{
    foreach_vertex (serial) {
      fprintf (fp, "%g %g ", x, y);
      for (scalar s in list)
	fprintf (fp, "%g ", s[]);
      fputc ('\n', fp);
    }
}


int main()
{

  /**
  We consider a periodic 1D domain. */
  
  init_grid (128);
  periodic (right);
  size (2);
  
  vertex scalar s[], w[], s2[];

    w[top]    = 0;
    w[bottom] = 0;
    w[right]  = 0;
    w[left]   = 0;

    s[top]    = 0;
    s[bottom] = 0;
    s[right]  = 0;
    s[left]   = 0;

    s2[top]    = 0;
    s2[bottom] = 0;
    s2[right]  = 0;
    s2[left]   = 0;


  /**
  We can optionally change the prolongation operator. We need to make
  sure that all fields use the same prolongation. */

#if 1
  for (scalar i in {s,w,s2}) {
    i.prolongation = refine_vert;
//    i.restriction = restriction_coarsen_vert2;
    i.restriction = restriction_vert;
  }
#endif

  /**
  We initialise the field with a function containing low-frequency and
  localised high-frequency components.  */
  
  foreach_vertex()
    s[] = ((sin(2.*pi*x) + 0.4*sin(15*pi*x)*max(sin(2.*pi*x), 0))*\
           (sin(2.*pi*y) + 0.4*sin(15*pi*y)*max(sin(2.*pi*y), 0)));

  /**
  The *w* field contains the wavelet transform of *s*. */

  wavelet_vertex (s, w);  

  /**
  We check that the inverse wavelet transform recovers the initial
  signal (almost) exactly. */

  inverse_wavelet_vertex (s2, w);
  foreach_vertex()
    assert (fabs (s2[] - s[]) < 1e-12);



  wavelet_vertex (s, w);
  for (int l = 0; l < 5; l++) {
    foreach_vertex_level (l)
      w[] = 0.;
    boundary_level ({w}, l);
  }
  inverse_wavelet_vertex (s2, w);
  write_vertex ({s,w}, "b", stderr);
  write_vertex ({s2,w}, "b", stdout);

  /**
  ~~~pythonplot Original signal
  import numpy as np
  import matplotlib.pyplot as plt
  log = np.loadtxt('log')
  s = log[:,2].reshape(129,129)
  plt.figure()
  plt.imshow(s,origin='lower')
  plt.savefig('./original.svg')
  ~~~

  ~~~pythonplot Filtered signal
  import numpy as np
  import matplotlib.pyplot as plt
  out = np.loadtxt('out')
  s2 = out[:,2].reshape(129,129)
  plt.figure()
  plt.imshow(s2,origin='lower')
  plt.savefig('./filter.svg')
  ~~~

  */

}

