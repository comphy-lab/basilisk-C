/**
# Unsymmetrical explosion 

We solve the Euler equations for a compressible gas under a
nonsymmetric initial pressure (density) distribution.  This problem
has been recently aborded by <a
href="https://arxiv.org/pdf/1701.00532.pdf">Eggers
et al. (2017)</a>. These authors showed that a selfsimilar structure exists
in the vicinity of the emerging shock wave. The shock spreads in the
transversal direction as $|t_o −t|^{1/2}$ and along the direction of
propagation as $|t_o − t|^{3/2}$. $t_o$ is the singularity time. */

#include "compressible.h"

#define LEVEL 8

int main() {

  /**
  We set Neumann conditions for normal velocity at all the boundaries. */

  foreach_dimension() {
    w.n[right] = neumann(0);
    w.n[left]  = neumann(0);
  }
  
  /**
  The domain spans $[-2:2]\times[-2:2]$. */

  origin (-2, -2);
  size (4.);
  DT = HUGE [0];
  init_grid (1 << (LEVEL));
  run(); 
}

/**
The initial density distribution is 
$$
\rho(x,y,0) = 0.2 + e^{-4(x^4 + y^2)}
$$
*/

event init (t = 0)
{  

  /**
  The momentum, $\mathbf{w} = \rho \mathbf{u}$, is initially null
  while the total energy is $E = \rho \mathbf{u}^2/2 + p/(\gamma-1)$.
  We assume that, initially, the isentropic relationship,
  $p/\rho^\gamma = cte$, holds. */

  foreach() {
    rho[] = 0.2 + exp(-4*(sq(y) + sq(x)*sq(x)));
    foreach_dimension()
      w.x[] = 0.;
    E[] = pow(rho[], gammao)/(gammao - 1.);
  }
}

/**
We wish to plot density distributions at different times. */

event plot (t = {0., 0.4, 0.511, 0.55}) {
  char name[80];
  sprintf (name, "density_%g", t);
  FILE * fp = fopen (name, "w");
  output_matrix(rho, fp, 1000); 
  fclose (fp);
}

/**
We will plot the time evolution of maximum entropy and the inverse of
the maximum gradient of the density. */ 

event graphs (i++) { 
  scalar rhograd[]; //gradient modulus 
  scalar S[]; //Entropy 
  foreach() { 
    double a = 0., w2 = 0.;
    foreach_dimension() { 
      a += sq(rho[1] -rho[-1]); 
      w2 += sq(w.x[]); 
    }
    rhograd[] = sqrt(a)/(2.*Delta); 
    double pres = (E[]-0.5*w2/rho[])*(gammao-1.); //pressure
    S[] = log(pres/(pow(rho[],(gammao)))); 
  } 
  stats sr = statsf (rhograd);
  stats ss = statsf (S);
  
  fprintf (stderr, "%g %g %g\n", t, 1./sr.max, ss.max); 
}

/**
We adapt the mesh by controlling the error on the density field. */

#if TREE
event adapt (i++) {
  adapt_wavelet ({rho}, (double[]){1e-4}, LEVEL + 3);
}
#endif

/**

## Results

We reproduce here figure 2 and 5 (sublots a and b) of the works of <a
href="https://arxiv.org/pdf/1701.00532.pdf">Eggers et
al. (2017)</a>. Note that figures of the paper were computed with a
largely finer mesh.

~~~gnuplot Density distribution for instants t = 0, 0.4, 0.511 and 0.55.
set pm3d
#set hidden
set cntrparam levels 7
set view 55,40
unset key
unset map
unset colorbox
set palette color
unset tics 
set xyplane relative -1.2

set contours
set multiplot layout 2,2
splot [0:2.][-1:1][0:1] 'density_0' binary u 1:2:3 with pm3d
splot [0:2.][-1:1][0:1] 'density_0.4' binary u 1:2:3 with pm3d
splot [0:2.][-1:1][0:1] 'density_0.511' binary u 1:2:3 with pm3d
splot [0:2.][-1:1][0:1] 'density_0.55' binary u 1:2:3 with pm3d
unset multiplot
~~~

~~~gnuplot Time evolution oth the maximum of the entropy and of the inverse of the density gradient.
set key left
set multiplot layout 2,1
set ylabel '|grad^{-1} ({/Symbol r})|_{max}'
set xtics 0.02
set ytics 0.5
plot [0.4:0.56][0:0.6] 'log' u 1:2 w l t 'Basilisk', 3.85*(0.511-x) t '|t - t_o| law'
set ylabel 'S_{max}'
set xlabel 't'
set ytics 0.005
plot [0.4:0.56][0:0.009] 'log' u 1:3 w l t 'Basilisk', 0.7*(x-0.511)**1.5 t '|t - t_o|^{3/2} law'
unset multiplot
~~~

*/
