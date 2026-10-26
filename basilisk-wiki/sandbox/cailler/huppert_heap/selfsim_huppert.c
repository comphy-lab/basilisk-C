/** 
# A self-similar solver for the viscous collapse of a heap 

This is an adaptation of the 
[code developed by *P.-Y. Lagrée*](/sandbox/M1EMN/Exemples/column_viscous.c) 
for exhibiting the self-similar behaviour of the collapse of 
a viscous heap for a newtonian fluid.

Here, we want to know if a [*self-similar solver*](/sandbox/cailler/self_sim_DNS/README) 
can be adaptated to this anisotropic problem, as self-similar coordinates'scalings 
differ in the two directions (see the problem transformed in self-similar coordinates 
in the below `selfsim_centered_huppert.h` file). 


~~~gnuplot Self-similar simulation of a Huppert's heap
reset session

set style line 1  lc rgb '#0025ad' lt 1 lw 1.5 # --- blue
set style line 2  lc rgb '#0060ad' lt 1 lw 1.5 #      .
set style line 3  lc rgb '#0099ad' lt 1 lw 1.5 #      .
set style line 4  lc rgb '#00ad88' lt 1 lw 1.5 #      .
set style line 5  lc rgb '#00ad4e' lt 1 lw 1.5 #      .
set style line 6  lc rgb '#09ad00' lt 1 lw 1.5 # --- green

set term pngcairo enhanced size 500,500
set output 'plot_selfsim_huppert.png'

set size ratio 1
set xrange [0:2]
set yrange [0:2]

p (9./10*(1.28338-x*x))**(1/3.) w l lw 2. lt -1 t "Analytic", \
  for [i=0:6] sprintf('out%d-0', i) u 1:2 w l ls i+1 t sprintf("τ = %d", i)
~~~
*/

#include "grid/multigrid.h"
#include "selfsim_centered_huppert.h"
#include "contact.h"
#include "selfsim_two-phase_huppert.h"
#include "tension.h"


// Domain extent
#define LDOMAIN 4.000001
//passive fluid
#define RHOF 1e-3
#define mug  1e-3
double mub;
double Bi=0; // *Bingham* number $Bi=(\frac{\tau_y L}{\mu_N U_0})$
double tmax;
double etamx;

#define LEVEL 7 // Maximum refinement 4./2**8 = 0.03125
#define LEVELmin 4 // Minimum refinement 2./2**4 = 0.125

scalar eta[];
scalar omega[]; // vorticity field
vector he[]; // height functions vector


/**
## Boundary conditions
*/

// Contact Angles
he.t[bottom] = contact_angle (pi/2.); 
he.t[left] = contact_angle (pi/2.);

// No-slip on the bottom
u.t[bottom] = dirichlet(0);
u.n[bottom] = dirichlet(0);

// Symmetry conditions for the half-heap
u.n[left] = dirichlet(0.);
u.t[left] = neumann(0.);

/** 
### Pressure boundary conditions

We shall ensure that
the normal component of the velocity is zero after projection. This
means that, at the boundary, the acceleration $\mathbf{\overline{a}}$ must be
balanced by the pressure gradient. But self-similar equations introduce new terms: 
a source term and a new velocity advection depending on the direction. 

#### Top and bottom boundaries

Taking care of boundary orientation
and staggering of $\mathbf{\overline{a}}$ and those new terms, this can be written
for top and bottom boundaries by solving this balance equation at the walls coming 
from the $\boldsymbol{e}_\eta$ projection of the momentum equation:
$$
\partial_\eta \overline{p} 
= 
- \dfrac{1}{\alpha_\eta(\tau)} 
\left(
  \dfrac{1}{5} \eta \, \partial_\eta \overline{u}_\eta
  - \dfrac{6}{5} \overline{u}_\eta
  + \mathrm{e}^{11 \tau / 5}
\right)
\quad ; \quad 
\alpha_\eta(\tau) = \mathrm{e}^{11 \tau / 5} / \rho
$$
First contribution of the right-hand side member corresponds to the non-zero 
advection term at the wall in the normal direction; the second contribution is 
the source term coming from the unsteady *Navier-Stokes* term.
Thir one is the gravity contribution.
*/

p[top] = neumann (
    neumann_pressure (ghost)
    - fm.n[ghost]/alpha.n[ghost] * (exp(2.2*t)
                                     + .2*y*(u.n[0,1] - u.n[])/Delta
                                     - 1.2*u.n[ghost])
); // alpha.n[ghost] already contains 'exp(2.2*t)*fm.y[]/rho(ff)'

p[bottom] = neumann (
    - neumann_pressure (0)
    + fm.n[0]/alpha.n[0] * (exp(2.2*t)
                            + .2*y*(u.n[] - u.n[0,-1])/Delta
                            - 1.2*u.n[0])
); // alpha.n[0] already contains 'exp(2.2*t)*fm.y[]/rho(ff)'


/** 
#### Left and right boundaries

In the $\boldsymbol{e}_\xi$ direction, the same principle is applied. 
It means that we need to verify:
$$
\partial_\xi \overline{p} 
= 
\dfrac{1}{\alpha_\xi(\tau)} 
\left(
  \dfrac{1}{5} \xi \, \partial_\xi \overline{u}_\xi
  + \dfrac{4}{5} \overline{u}_\xi
\right)
\quad ; \quad 
\alpha_\xi(\tau) = \mathrm{e}^{7 \tau / 5} / \rho
$$
*/

p[left] = neumann (
    - neumann_pressure(0) 
    - fm.n[0]/alpha.n[0] * (.2*x*(u.n[]-u.n[-1])/Delta
                             + 0.8*u.n[0])
); // alpha.n[0] already contains 'exp(1.4*t)*fm.x[]/rho(ff)'

p[right] = neumann (
    neumann_pressure(ghost) 
    + fm.n[ghost]/alpha.n[ghost] * (.2*x*(u.n[1]-u.n[])/Delta
                                     + 0.8*u.n[ghost])
); // alpha.n[ghost] already contains 'exp(1.4*t)*fm.x[]/rho(ff)'


/**
## Main Parameters
*/
int main() {
  L0 = LDOMAIN;
  N = 1 << LEVEL; // number of grid points
  DT = 0.025/4.; // maximum timestep
  TOLERANCE = 1e-3;
  
  /**
   **Numerical values**
    
  The problem is a heap of viscoplastic material of height = length
    released on a flat plate.
    We may define a cone for Abrams slup test (concrete).
  */
  // Initial conditions
#define Ltas 0.9999999999
#define Ls 0.9999999999
#define Htas 0.99999999999

/**
The fluid has a density $\rho_1$.
We define a quantity, $\mu_b$ which is in fact the inverse of
$Re_1=\frac{\rho_0 U_0 L}{\mu_N }$ and which is here equal to 1.
*/
    
  mub= 1;
  etamx=1000.;
  tmax = 6;
  DT = 0.025;

  // Works quite well up to tau = 6 with:
  // rho2=1.e-1 or 1.e-2, and mu2=1.e-2 or 1.e-3 
  // Best combo was rho2 = 1.e-1, mu2 = 1.e-3 !!!!!!!!!!  

  // Fails completely for rho=1.e-3.

  rho1 = 1., rho2 = 1.e-1; //1.e-3 ;
  mu1 = 1., mu2 = 1.e-3; //1.e-3 ;

  run();   
}

/**
## Main Events
*/
event init (t = 0) {
// Defining an initial trapezoidal shape:
  fraction (f, min(Htas - y, Ltas - (fabs(x)+ y * (Ltas-Ls))));
}

event vorti (i++){
  vorticity (u, omega);
}

#if QUADTREE
// if no #include "grid/multigrid.h" then adapt
event adapt (i++) {
  adapt_wavelet({f},(double[]){0.001}, maxlevel = LEVEL, minlevel = LEVELmin);
}
#endif

event stop (t = tmax);


event monitoring (i++){
  printf("time = %g ; iteration #%i\n", t, i);
  fflush(stdout); // empty buffer
}


/**
## Outputs

**Convergence outputs:**
*/
void mg_print (mgstats mg)
{
#if 0
    if (mg.i > 0 && mg.resa > 0.)
        fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
                 exp (log (mg.resb/mg.resa)/mg.i));
#else
  if (mg.i > 0 && mg.resa > 0.)
    fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
              mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0); 
#endif
}

event logfile (i++) {
  stats s = statsf (f);
  fprintf (stderr, "%g %d %g %g %g %g\n", t, i, dt, s.sum, s.min, s.max - 1.);
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
  fflush (stderr);

  double max_u = 0.;
  double max_p = 0.;
  int huge_u = 0;
  int huge_p = 0;

  foreach (reduction(max:max_u)
          reduction(max:max_p)
          reduction(+:huge_u)
          reduction(+:huge_p)) {

    if (isfinite(u.x[]) && isfinite(u.y[])) {
      double umag = sqrt(sq(u.x[]) + sq(u.y[]));
      if (umag > max_u)
        max_u = umag;
      if (umag > 1e10)
        huge_u++;
    }

    if (isfinite(p[])) {
      if (fabs(p[]) > max_p)
        max_p = fabs(p[]);
      if (fabs(p[]) > 1e10)
        huge_p++;
    }
  }

  fprintf(stderr,
          "tau=%g max_u=%g max_p=%g huge_u=%d huge_p=%d\n",
          t, max_u, max_p, huge_u, huge_p);
}

/** 
## Plots
 **Some interface profiles for comparison:**
*/

event pij0 (t = 0.) {
  char s[80];
  sprintf (s, "out0-%g",Bi);
  FILE * fp = fopen (s, "w");
  output_facets (f, fp);
  fclose (fp);
}
event pij1 (t = 1) {
  char s[80];
  sprintf (s, "out1-%g",Bi);
  FILE * fp = fopen (s, "w");
  output_facets (f, fp);
  fclose (fp);
}
event pij2 (t = 2) {
  char s[80];
  sprintf (s, "out2-%g",Bi);
  FILE * fp = fopen (s, "w");
    output_facets (f, fp);
  fclose (fp);
}
event pij3 (t = 3){
  char s[80];
  sprintf (s, "out3-%g",Bi);
  FILE * fp = fopen (s, "w");
  output_facets (f, fp);
  fclose (fp);
}
event pij4 (t = 4){
  char s[80];
  sprintf (s, "out4-%g",Bi);
  FILE * fp = fopen (s, "w");
  output_facets (f, fp);
  fclose (fp);
}
event pij5 (t = 5){
  char s[80];
  sprintf (s, "out5-%g",Bi);
  FILE * fp = fopen (s, "w");
  output_facets (f, fp);
  fclose (fp);
}
event pijtmax (t = tmax){
  char s[80];
  sprintf (s, "out6-%g",Bi);
  FILE * fp = fopen (s, "w");
  output_facets (f, fp);
  fclose (fp); 
}


/**
 **Saving the top position and runout position for slump measurements:**
*/
vector h[];
event timeseries (t += 0.1 ) {
  heights (f, h);
  double maxy = - HUGE,maxx = - HUGE;;
  foreach()
  if ((h.y[] != nodata) && (h.x[] != nodata)) {
    double yi = y + height(h.y[])*Delta;
    double xi = x + height(h.x[])*Delta;
    if (yi > maxy)
      maxy = yi;
    if (xi > maxx)
      maxx = xi;
  }
  char s[80];
  sprintf (s, "hmax-%g",Bi);
  static FILE * fp0 = fopen (s, "w");
  fprintf (fp0, "%g %g %g\n", t, maxx, maxy);
  fflush (fp0);
}

/** 
## MOVIE

We make a movie of the self-similar simulation:
*/
#include "view.h"
#include "cmaps.h"

event mov_p_omega (t+=1e-2) {
  char legend[100];
  sprintf(legend, "ln(t) = %0.2g", t);
  view (tx = 0, ty =-0.4, width = 1600., height = 800., fov = 18.); 
  box();
  draw_vof ("f", lw = 1.5);
  squares("omega", spread=-1, linear=true, map=viridis);
  vectors("u", scale=0.05, lc = {1, 1, 1});
  mirror ({1,0}) {
    squares ("p", spread=-1, linear=true, map=cool_warm);
    draw_vof ("f", lw = 1.5);
  }
  draw_string(legend, 2, size = 20., lw = 2.1);
  save ("selfsim_huppert_movie.mp4");
}

/**
![Collapsing Huppert's heap simulated in the self-similar space | *left* pressure | *right:* vorticity field | $N_{max} = 7$](selfsim_huppert/selfsim_huppert_movie.mp4)(width="40%")
*/


/** 
## Bibliography
 

* [Huppert H.](http://www.itg.cam.ac.uk/people/heh/Paper47.pdf)
 ”The propagation of two-dimensional and axisymmetric viscous gravity currents over a rigid horizontal surface” J . Fluid Mech. (1982), vol. 121, p p . 43-58
 * Lagrée  [M1EMN
Master 1 Ecoulements en Milieu Naturel](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
 * Lagrée  [M2EMN
Master 2 Ecoulements en Milieu Naturel](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/mainM2EMN.pdf)


 
*/