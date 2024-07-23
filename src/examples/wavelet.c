/**
# Wavelet transforms and filtering

This simple example illustrates how to compute wavelet transforms and
perform filtering and mesh adaptation. We do this in one dimension,
using bi-trees. */

#include "grid/bitree.h"

/**
A simple function to output fields on each level. */

void write_level (scalar * list, const char * tag, FILE * fp)
{
  for (int l = 0; l <= depth(); l++)
    foreach_level (l, serial) {
      fprintf (fp, "%s%d %g ", tag, l, x);
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
  
  scalar s[], w[], s2[];

  /**
  We can optionally change the prolongation operator. We need to make
  sure that all fields use the same prolongation. */

#if 0
  for (scalar i in {s,w,s2})
    i.refine = i.prolongation = refine_biquadratic;
#endif

  /**
  We initialise the field with a function containing low-frequency and
  localised high-frequency components.  */
  
  foreach()
    s[] = sin(2.*pi*x) + 0.4*sin(15*pi*x)*max(sin(2.*pi*x), 0);

  /**
  The *w* field contains the wavelet transform of *s*. */

  wavelet (s, w);  

  /**
  We check that the inverse wavelet transform recovers the initial
  signal (almost) exactly. */

  inverse_wavelet (s2, w);
  foreach()
    assert (fabs (s2[] - s[]) < 1e-12);
  write_level ({s,w}, "", stderr);

  /**
  ## The original signal and its wavelet transform.
     
  ~~~gnuplot Original signal representation on each level
  maxlevel = 6
  spacing = 1.5
  set ytics -maxlevel,1,0
  set xlabel 'x'
  set ylabel 'Level'
  set grid ytics ls 0
  max(a,b) = a > b ? a : b
  set samples 500
  f(x)=sin(2.*pi*x)+0.4*sin(15.*pi*x)*max(sin(2.*pi*x), 0);
  plot for [i=0:maxlevel] '< grep ^'.i.' log' u 2:($3/spacing - i) \
                             w lp t '' lt 7,			   \
       for [i=0:maxlevel] f(x)/spacing - i t '' lt 1
  ~~~

  ~~~gnuplot Corresponding wavelet decomposition
  plot for [i=0:maxlevel] '< grep ^'.i.' log' u 2:($4/spacing - i) \
                             w lp t '' lt 7,			   \
       for [i=0:maxlevel] f(x)/spacing - i t '' lt 1
  ~~~
  
  ## Localised low-pass filtering
  
  We want to filter the high-frequency components of the signal, but
  only locally i.e. for $x > 1$. To do this we cancel the
  high-frequency (i.e. high-level) wavelet coefficients and recover
  the corresponding filtered signal *s2* with the inverse wavelet
  transform. */

  for (int l = 5; l <= 6; l++) {
    foreach_level (l)
      w[] *= (x <= 1);
    boundary_level ({w}, l);
  }
  inverse_wavelet (s2, w);
  write_level ({s2,w}, "a", stdout);

  /**
  ~~~gnuplot Low-pass filtered reconstruction (for $x > 1$)
  plot for [i=0:maxlevel] '< grep ^a'.i.' out' u 2:($3/spacing - i) \
                             w lp t '' lt 7,			   \
       for [i=0:maxlevel] f(x)/spacing - i t '' lt 1
  ~~~

  ## High-pass filtering

  Here we filter out the low-frequency components everywhere. */

  wavelet (s, w);
  for (int l = 0; l < 5; l++) {
    foreach_level (l)
      w[] = 0.;
    boundary_level ({w}, l);
  }
  inverse_wavelet (s2, w);
  write_level ({s2,w}, "b", stdout);

  /**
  ~~~gnuplot High-pass filtered reconstruction
  plot for [i=0:maxlevel] '< grep ^b'.i.' out' u 2:($3/spacing - i) \
                             w lp t '' lt 7,			   \
       for [i=0:maxlevel] f(x)/spacing - i t '' lt 1
  ~~~

  ## Mesh adaptation

  We can use the wavelet coefficients to decide where to coarsen the mesh. */

  wavelet (s, w);
  unrefine (fabs(w[]) < 0.1);
  write_level ({s,w}, "c", stdout);
}

/**
~~~gnuplot Adapted mesh
plot for [i=0:maxlevel-1] '< grep ^c'.i.' out' u 2:($3/spacing - i)	\
        w lp t '' lt 7,							\
      '< grep ^c6 out' u 2:($3/spacing - 6) w p t '' lt 7,		\
      for [i=0:maxlevel] f(x)/spacing - i t '' lt 1
~~~
  
## Discussion
  
* Note that the restriction operators assume that the discrete
  quantity is the average of the function over the cell, not its
  point value. While this makes no difference for second-order
  (i.e. linear) prolongation operators, it matters for higher-order
  prolongation but is not taken into account in this example.  
* Boundary conditions for the wavelet transform are important and
  non-trivial. This is not an issue here since we use a periodic
  domain.
* The link between the transform implemented here and formal wavelets 
  is not entirely clear. It is probably a kind of 
  [second-generation wavelet](#sweldens1998).
* The filters applied in the example could be formulated in a more 
  generic way using a filtering function $\phi(x,l)$ and a code 
  looking like:

~~~literatec
  foreach_cell()
    w[] *= phi(x,level);
~~~

combined with proper application of boundary conditions.

## References

~~~bib
@article{sweldens1998,
  title={The lifting scheme: A construction of second generation wavelets},
  author={Sweldens, Wim},
  journal={SIAM journal on mathematical analysis},
  volume={29},
  number={2},
  pages={511--546},
  year={1998},
  publisher={SIAM}
}
~~~
*/
