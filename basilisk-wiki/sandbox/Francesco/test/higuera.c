/**
# Viscous hydraulic Jump

We want to reproduce the hydraulic jump of Higuera (1994). */

#include "grid/cartesian1D.h"
#include "saint-venant.h"

#define h1    0.1
#define x1    0.705*h1
#define Hw    1.8138
#define Uw    0.8964
#define h2    0.00001

double S;

int main() {
  X0 = x1;
  L0 = 1.-x1;
  nl = 30;
  nu = 1.;
  N = 256;

  for (S = 0.5; S <= 2; S *= 2) {
    G  = S;
    run();
  }
}

/**
We impose boundary condition for $h$ and $\eta$. */

h[left] = dirichlet (h1);
eta[left] = dirichlet (h1);

h[right] = dirichlet (h2*1);
eta[right] = dirichlet (h2*1);

/**
## Initialization

We define a scalar field *hc* to check for convergence on $h$. */

scalar hc[];

event init (i = 0) {

  /**
  We set a constant velocity at the inlet and a free outlet. */

  for (vector u in ul) {
    u.n[left] = h[left] ? 1./h[left]*1.58977503178463 : 0.; // some magic
    u.n[right] = neumann(0.);
  }
  
  /**
  We initialize *h* and *hc*. */
  
  foreach() 
    hc[] = h[] = h1;
}

/**
We check for convergence. */

#if 1
event logfile (t += 0.1; i <= 1000000) {
  double dh = change (h, hc);
  printf ("%g %g\n", t, dh);
  if (i > 0 && dh < 1e-4)
    return 1;
}
#endif

/**
Uncomment this part if you want on-the-fly animation. */

#if 0
event output (t += 0.01) {
  static FILE * fp = popen ("gnuplot", "w");
  fprintf (fp,
           "set title 'N=%d, nl=%d, t=%f'\n"
           "set xl 'x'\nset yl 'h'\nset key top left\n"
           "plot [][] '-' u 1:2 w l t 'h', 1.8136*x w l t 'watson'\n",
	   N, nl, t); 
  foreach()
    fprintf (fp, "%g %g\n", x, h[]);
  fprintf (fp, "e\n");
  fflush (fp);
}
#endif

/**
## Output

We print *h* and the velocity component. */

event output (t = end) {
  char name[80];
  sprintf (name,"htau-%d-%d-%g",N,nl,S);
  FILE *fp = fopen(name,"w");
  vector u0 = ul[0];
  foreach() 
    fprintf (fp,"%g %g %g\n", x, h[], 2*nl*u0.x[]/h[]);
  fclose(fp);
  // velocity field
  char name2[80];
  sprintf (name2,"vel-%d-%d-%g",N,nl,S);
  FILE *fp1 = fopen(name2,"w");
  foreach() {
    int l = 0;
    double sumh = 0.;
    for (vector u in ul) {
      sumh += h[]*layer[l];
      double z = sumh - h[]*layer[l++]*0.5;
      fprintf (fp1,"%g %g %g\n",x,z/h[],u.x[]);
    }
    fprintf (fp1,"\n\n");
  }
  fclose(fp1);
}

#if 0
event output (i+=10) {
  // fprintf(stdout,"step: %d t: %g dt: %g\n",i,t,dt);
  double sum = 0;
  int j = 0;
  foreach() {
    if (j++ == N/8) {
      int l = 0;
      for (vector u in ul) {
	sum += u.x[]*h[]*layer[l++];
      }
      fprintf(stdout,"step: %d sum: %g\n",j,sum);
    }
  }
}
#endif

/**
## Results

~~~gnuplot Comparison with Figure 2 of Higuera (1994).
X0 = 169
X1 = 604
Y0 = 222.24
Y1 = 528
unset tics
plot [0:][0:605] 'higuera.png' binary filetype=png with rgbimage not, \
  'htau-256-30-0.5' u (X0+$1*(X1-X0)):($2/2*(Y1-Y0)+Y0) t 'h S = 0.5' w l, \
  'htau-256-30-0.5' u (X0+$1*(X1-X0)):($3/15*(Y1-Y0)+Y0) t 'tau S = 0.5' w l, \
  'htau-256-30-1' u (X0+$1*(X1-X0)):($2/2*(Y1-Y0)+Y0) t 'h S = 1' w l, \
  'htau-256-30-1' u (X0+$1*(X1-X0)):($3/15*(Y1-Y0)+Y0) t 'tau S = 1' w l, \
  'htau-256-30-2' u (X0+$1*(X1-X0)):($2/2*(Y1-Y0)+Y0) t 'h S = 2' w l, \
  'htau-256-30-2' u (X0+$1*(X1-X0)):($3/15*(Y1-Y0)+Y0) t 'tau S = 2' w l
~~~

## Bibliography

* Higuera, F. 1994. [The hydraulic jump in a viscous laminar
flow](http://dx.doi.org/10.1017/S0022112094002041). J. Fluid
Mech. 274, 69–92.
*/
