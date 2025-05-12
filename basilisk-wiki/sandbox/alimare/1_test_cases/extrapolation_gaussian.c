/**
# 2D extrapolation 

We assess here the order of the extrapolation function
*embed_extrapolate* and compare it with a former version that I kept from
Ghigo's sandbox, which I
called *embed_extrapolate2*.
This case is similar to extrapolation_ghigo, this time with constant resolution
and an initial exponential (Gaussian) field. */

#include "grid/multigrid.h"
#include "../../ghigo/src/myembed.h"
#include "../embed_extrapolate_3.h"
#include "utils.h"
#include "view.h"
/**
## Exact solution

We define here the exact solution, evaluated at the center of each
cell. The solution is linear so the extrapolation should be exact
(since it is 2nd order and the restriction/prolongation functions are
second-order). */

static double exact (double x, double y)
{
  double theta = atan2(y, x), r = sqrt (sq(x) + sq(y));
  double phi = (0.3 + 0.15*cos(6.*theta)) - r;
  return exp(-sq(phi));
}

/**
We also define the shape of the domain. */

#if GEOM == 1 // concave
void p_shape (scalar c, face vector f)
{
  vertex scalar phi[];
  foreach_vertex() {
    double theta = atan2(y, x), r = sqrt (sq(x) + sq(y));
    phi[] = r - (0.25 + 0.05*cos(6.*theta));
  }
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f,
         smin = 1.e-14, cmin = 1.e-14);
}
#elif GEOM == 2 // convex
void p_shape (scalar c, face vector f)
{
  vertex scalar phi[];
  foreach_vertex() {
    double theta = atan2(y, x), r = sqrt (sq(x) + sq(y));
    phi[] = (0.25 + 0.05*cos(6.*theta)) - r;
  }
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f,
         smin = 1.e-14, cmin = 1.e-14);  
}
#elif GEOM == 3 // concave
void p_shape (scalar c, face vector f)
{
  vertex scalar phi[];
  foreach_vertex() {
    double theta = atan2(y, x), r = sqrt (sq(x) + sq(y));
    phi[] = r - (0.3 + 0.15*cos(6.*theta));
  }
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f,
         smin = 1.e-14, cmin = 1.e-14);  
}
#else // convex
void p_shape (scalar c, face vector f)
{
  vertex scalar phi[];
  foreach_vertex() {
    double theta = atan2(y, x), r = sqrt (sq(x) + sq(y));
    phi[] = (0.3 + 0.15*cos(6.*theta)) - r;
  }
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f,
         smin = 1.e-14, cmin = 1.e-14);  
}
#endif // GEOM


#if TREE
/**
When using *TREE*, we try to increase the accuracy of the restriction
operation in pathological cases by defining the gradient of *a* at the
center of the cell. */

void a_embed_gradient (Point point, scalar s, coord * g)
{
  g->x = -2*x*exp(-sq(x))+1;
  g->y = -2.+2*y;
}
#endif // TREE

/**
## Setup */

int lvl;

double myboundary_condition(Point point){
  coord n, p;
  embed_geometry (point, &p, &n);
  return exact(x+Delta*p.x,y+Delta*p.y);
}

int main()
{
  /**
  The domain is $1\times 1$. */
  
  origin (-L0/2., -L0/2.);

  for (lvl = 5; lvl <= 11; lvl++) { // minlevel = 4
    fprintf(stderr, "### %d\n", lvl);
    /**
    We initialize the grid. */

    N = 1 << (lvl);
    init_grid (N);
    
    /**
    We initialize the embedded boundary. */

#if TREE
    /**
    When using *TREE*, it is necessary to modify the refinement and
    prolongation of the volume and face fractions to account for
    embedded boundaries. In this case, using *embed_fraction_refine*
    for prolongation (as well as refine) seems to overall improve the
    accuracy of the results. */

    cs.refine = embed_fraction_refine;
    cs.prolongation = fraction_refine;

    foreach_dimension()
      fs.x.prolongation = embed_face_fraction_refine_x;
#endif // TREE
    
    cm   = cs;
    fm   = fs;
    csm1 = cs;

#if TREE
    /**
    When using *TREE*, we refine the mesh around the embedded
    boundary. */

    astats ss;
    int ic = 0;
    do {
      ic++;
      p_shape (cs, fs);
      ss = adapt_wavelet ({cs}, (double[]) {1.e-30},
        maxlevel = (lvl), minlevel = (lvl) - 2);
    } while ((ss.nf || ss.nc) && ic < 50);
#endif // TREE
    
    p_shape (cs, fs);

    /**
    We define *a* and set its boundary conditions. */
    
    scalar a[],a2[];
    foreach(){
      a[] = cs[] > 0. ? exact (x, y) : nodata;
      a2[] = cs[] > 0. ? exact (x, y) : nodata;
    }

#if GEOM == 1 || GEOM == 3
    a[left]   = exact (x - Delta/2., y);
    a[right]  = exact (x + Delta/2., y);
    a[top]    = exact (x, y + Delta/2.);
    a[bottom] = exact (x, y - Delta/2.);

    a2[left]   = exact (x - Delta/2., y);
    a2[right]  = exact (x + Delta/2., y);
    a2[top]    = exact (x, y + Delta/2.);
    a2[bottom] = exact (x, y - Delta/2.);
#endif // GEOM == 1 || GEOM == 3

    a[embed] = dirichlet  (exact(x,y));
    a2[embed] = dirichlet (exact(x,y));

#if TREE
    /**
    On *TREE*, we also modify the prolongation an restriction
    operators for *a*. */
    
    a.refine = a.prolongation = refine_embed_linear;
    a.restriction = restriction_embed_linear;
    a2.refine = a2.prolongation = refine_embed_linear;
    a2.restriction = restriction_embed_linear;
#endif // TREE
    
    /**
    We use "third-order" [face flux interpolation](/src/embed.h). */

#if ORDER2
    a.third = a2.third = false;
#else
    a.third = a2.third = true;
#endif // ORDER2

#if TREE
    /**
    We also define the gradient of *a* at the full cell center of
    cut-cells. */

    a.embed_gradient = a2.embed_gradient = a_embed_gradient;
#endif // TREE
            
    boundary ({a,a2});

    /**
    We now compute *a* in each cut-cell using the *embed_extrapolate
    function*. */

    foreach() {
      if (cs[] > 0. && cs[] < 1.) {
        coord b, nn;
        embed_geometry (point, &b, &nn);
        bool dirichlet;
        double ab = (a.boundary[embed] (point, point, a, &dirichlet));
        double ab2 = myboundary_condition(point);
        assert (dirichlet);
        a[] = embed_extrapolate (point, a, cs, nn, ab);
        a2[] = embed_extrapolate2 (point, a2, nn, b, ab2);
      }
    }
    boundary ({a,a2});
        
    /**
    The total (*e*), partial cells (*ep*) and full cells (*ef*) errors
    fields and their norms are computed. */
    
    scalar e[], ep[], ef[];
    scalar e2[], ep2[], ef2[];
    foreach() {
      if (cs[] == 0.)
  ep[] = ef[] = e[] = ep2[] = ef2[] = e2[] = nodata;
      else {
        e[]   = fabs (a[] - exact (x, y));
        e2[]  = fabs (a2[] - exact (x, y));
        ep[]  = cs[] < 1. ? e[] : nodata;
        ef[]  = cs[] >= 1. ? e[] : nodata;
        ep2[] = cs[] < 1. ? e2[] : nodata;
        ef2[] = cs[] >= 1. ? e2[] : nodata;
      }
    }
    norm n = normf (e), np = normf (ep), nf = normf (ef);
    norm n2 = normf (e2), np2 = normf (ep2), nf2 = normf (ef2);
    fprintf (stderr, "%d %.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3g\n",
       N, n.avg, n.max, np.avg, np.max, nf.avg, nf.max,
       n2.avg, n2.max, np2.avg, np2.max, nf2.avg, nf2.max);
    fflush (stderr);
    
    /**
    The solution and the error are displayed using bview. */

    if (lvl == 7) {
      clear ();
      view (fov = 18);
      
      squares ("a", spread = -1);
      draw_vof ("cs", "fs");
      save ("a.png");
      
      squares ("a", spread = -1);
      draw_vof ("cs", "fs");
      save ("a2.png");

      squares ("e", spread = -1);
      draw_vof ("cs", "fs");
      save ("e.png");

      squares ("e2", spread = -1);
      draw_vof ("cs", "fs");
      save ("e2.png");

      cells ();
      draw_vof ("cs", "fs");
      save ("mesh.png");
    }
  }
}

/**
## Results for $\Omega_4$

![Solution for *l=7*](extrapolation_gaussian/a.png)

![Error for *l=7*](extrapolation_gaussian/e.png)

![Error for *l=7*](extrapolation_gaussian/e2.png)

![Mesh for *l=7*](extrapolation_gaussian/mesh.png)

#### Errors

~~~gnuplot Average error convergence
set terminal svg font ",16"
set key bottom left spacing 1.1
set xtics 32,4,4096
set grid ytics
set ytics format "%.0e" 1.e-20,10000,1.e2 
set xlabel 'N'
set ylabel '||error||_{1}'
set xrange [32:4096]
set yrange [1.e-20:*]
set logscale
plot 'log' u 1:4 w p ps 1.25 pt 7 lc rgb "black" t 'cut-cells', \
     ''    u 1:6 w p ps 1.25 pt 5 lc rgb "blue" t 'full cells', \
     ''    u 1:2 w p ps 1.25 pt 2 lc rgb "red" t 'all cells',\
     ''    u 1:10 w p ps 1.25 pt 7 lc rgb "green" t 'cut-cells quad', \
     ''    u 1:12 w p ps 1.25 pt 5 lc rgb "cyan" t 'full cells quad', \
     ''    u 1:8 w p ps 1.25 pt 2 lc rgb "magenta" t 'all cells quad'
~~~

~~~gnuplot Maximum error convergence
set ylabel '||error||_{inf}'
plot '' u 1:5 w p ps 1.25 pt 7 lc rgb "black" t 'cut-cells', \
     '' u 1:7 w p ps 1.25 pt 5 lc rgb "blue" t 'full cells', \
     '' u 1:3 w p ps 1.25 pt 2 lc rgb "red" t 'all cells',\
     '' u 1:11 w p ps 1.25 pt 7 lc rgb "green" t 'cut-cells quad', \
     '' u 1:13 w p ps 1.25 pt 5 lc rgb "cyan" t 'full cells quad', \
     '' u 1:9 w p ps 1.25 pt 2 lc rgb "magenta" t 'all cells quad'
~~~

#### Order of convergence

The moment of truth...

~~~gnuplot Order of convergence of the  average error
reset
set terminal svg font ",16"
set key bottom right spacing 1.1
set xtics 32,4,4096
set grid ytics
set xlabel 'N'
set ylabel 'Order of ||error||_{1}'
set xrange [8:4096]
set yrange [*:*]
set ytics format "%.0e" 1.e-10,10000,1.e2 
set logscale

# Average order of converge

ftitle(b) = sprintf(", avg order n=%4.2f", -b);
ftitle2(b) = sprintf(", avg order quad n=%4.2f", -b);

f1(x) = a1 + b1*x; # cut-cells
f3(x) = a3 + b3*x; # all cells
f2(x) = a2 + b2*x; # cut-cells
f4(x) = a4 + b4*x; # all cells


fit [*:*][*:*] f1(x) '< sort -k 1,1n log' u (log($1)):(log($4)) via a1,b1; # cut-cells 
fit [*:*][*:*] f3(x) '< sort -k 1,1n log' u (log($1)):(log($2)) via a3,b3; # all cells
fit [*:*][*:*] f2(x) '< sort -k 1,1n log' u (log($1)):(log($10)) via a2,b2; # cut-cells 
fit [*:*][*:*] f4(x) '< sort -k 1,1n log' u (log($1)):(log($8)) via a4,b4; # all cells


plot '' u 1:4 w lp ps 1.25 pt 7 lc rgb "black" t 'cut-cells'.ftitle(b1), \
     '' u 1:2 w lp ps 1.25 pt 2 lc rgb "red" t 'all cells'.ftitle(b3),\
     '' u 1:10 w lp ps 1.25 pt 7 lc rgb "blue" t 'cut-cells'.ftitle2(b2), \
     '' u 1:8 w lp ps 1.25 pt 2 lc rgb "green" t 'all cells'.ftitle2(b4)
~~~

~~~gnuplot Order of convergence of the maximum error
set ylabel 'Order of ||error||_{inf}'

# Average order of convergence

fit [*:*][*:*] f1(x) '< sort -k 1,1n log' u (log($1)):(log($5)) via a1,b1; # cut-cells 
fit [*:*][*:*] f3(x) '< sort -k 1,1n log' u (log($1)):(log($3)) via a3,b3; # all cells
fit [*:*][*:*] f2(x) '< sort -k 1,1n log' u (log($1)):(log($11)) via a2,b2; # cut-cells 
fit [*:*][*:*] f4(x) '< sort -k 1,1n log' u (log($1)):(log($9)) via a4,b4; # all cells

plot '' u 1:5 w lp ps 1.25 pt 7 lc rgb "black" t 'cut-cells'.ftitle(b1), \
     '' u 1:3 w lp ps 1.25 pt 2 lc rgb "red" t 'all cells'.ftitle(b3),\
     '' u 1:11 w lp ps 1.25 pt 7 lc rgb "blue" t 'cut-cells'.ftitle2(b2), \
     '' u 1:9 w lp ps 1.25 pt 2 lc rgb "green" t 'all cells'.ftitle2(b4)
~~~

Gnuplot has spoken, my method gives better result in this specific case, which
is similar to my dendritic growth simulation.
It means that there is a strong requirement of local constant resolution in the
vicinity of the interface.

## References

~~~bib
@article{johansen1998,
  title={A Cartesian grid embedded boundary method for Poisson's
  equation on irregular domains},
  author={Johansen, Hans and Colella, Phillip},
  journal={Journal of Computational Physics},
  volume={147},
  number={1},
  pages={60--85},
  year={1998},
  publisher={Elsevier},
  url={https://pdfs.semanticscholar.org/17cd/babecd054d58da05c2ba009cccb3c687f58f.pdf}
}
~~~
*/
