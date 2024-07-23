/**
# Breakup of a rectangular perturbation into a train of solitons

We reproduce the study of [Madsen et al,
2008](/src/references.bib#madsen2008). An initial
rectangular perturbation of width $2b$ and height $a$ is superposed on
an ocean of constant depth $h_0$. 

We solve the one-dimensional problem with an adaptive "bitree". We
compute both the Saint-Venant solution and the Green-Naghdi
solutions. */

#include "grid/bitree.h"
#if SAINT_VENANT
# include "saint-venant.h"
#else
# include "green-naghdi.h"
#endif

/**
We need a very high level of refinement to accurately reproduce the
results of Madsen et al. The depth $h_0$ and gravity are both set to
unity (Madsen et al use these non-dimensional units). The amplitude is
set to 0.1 and half-width to 12.2 as done in section 3.2.1 of [Madsen
et al, 2008](/src/references.bib#madsen2008). */

#define MAXLEVEL 15
double h0 = 1, a = 0.1, b = 12.2;

int main() {
  N = 1 << MAXLEVEL;
  L0 = 4000. [1];
  G = 1.;
#if !SAINT_VENANT
  gradient = NULL;
#endif
  run();
}

/**
The initial conditions are a (half)rectangle sitting on the axis of
symmetry at $x=0$. */

event init (i = 0) {
  foreach() {
    h[] = h0 + a*(x < b);
    zb[] = - h0;
  }
}

/**
We save the free-surface elevation at times $t\sqrt{g/h_0}=$ 35, 90,
700 and 2000 as in Figure 1 of [Madsen et al,
2008](/src/references.bib#madsen2008). We use a frame of reference
travelling with the wave i.e. we use $(x-b)/h_0-t\sqrt{g/h_0}$ as
coordinate. */

event plot (t = {35, 90, 700, 2000}) {
  char name[80];
  sprintf (name, "t-%g", t);
  FILE * fp = fopen (name, "w");
  foreach (serial)
    fprintf (fp, "%g %g\n", x - b - t*sqrt(G*h0), eta[]);
  fclose (fp);
}

/**
We need a low error threshold for adaptation on $\eta$. */

#if TREE
event adapt (i++) {
  astats s = adapt_wavelet ({eta}, (double[]){1e-6}, MAXLEVEL);
  fprintf (stderr, "%g refined %d cells, coarsened %d cells\n",
	   t, s.nf, s.nc);
}
#endif

/**
Finally we compare the two results (for Green-Naghdi and Saint-Venant)
with the results of [Madsen et al,
2008](/src/references.bib#madsen2008), Figure 1. The GN solutions
reproduces the breakup of the initial perturbation into a train of
solitary waves at long times, while the Saint-Venant solution becomes
very inaccurate.

The agreement between Madsen et al and the GN solution is quite good
given that the model of Madsen et al is a higher-order Boussinesq
expansion than the GN model ([Madsen et al, 
2006](/src/references.bib#madsen2006)).

~~~gnuplot Evolution of surface elevation
set term svg enhanced size 640,480 font ",9"
set multiplot layout 4,1 scale 1.,1.
set ytics 0,0.05,0.1
set yrange [-0.02:0.1]
set xtics -50,10,50
plot [-50:50]'../madf1a.plot' w l t 'Madsen et al, 2008', \
             '< sort -n -k1,2 t-35' w l t 'SGN', \
             '< sort -n -k1,2 ../madsen-sv/t-35' w l t 'Saint-Venant'
set xtics -100,20,100
set xrange [-100:100]
unset key
plot '../madf1b.plot' w l, \
     '< sort -n -k1,2 t-90' w l, \
     '< sort -n -k1,2 ../madsen-sv/t-90' w l
plot '../madf1c.plot' w l, \
     '< sort -n -k1,2 t-700' w l, \
     '< sort -n -k1,2 ../madsen-sv/t-700' w l
plot '../madf1d.plot' w l, \
     '< sort -n -k1,2 t-2000' w l, \
     '< sort -n -k1,2 ../madsen-sv/t-2000' w l
unset multiplot
~~~
*/
