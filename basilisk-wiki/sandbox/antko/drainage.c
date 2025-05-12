/**
# Viscous drainage of a liquid film into a Plateau border 

Foams and emulsions involve bubbles and drops, but also liquid films that delimit these unit cells. As time flows, these liquid films drain, and burst, yielding foam/emulsion ageing and ripening.

We here track the evolution of such a thinning viscous liquid film. Typically, these films are connected to meniscii or *Plateau borders*. Because of the negative curvature of the interface there, these Plateau borders are low-pressure regions that draw out the liquid contained in the films, as illustrated in this sequence:

~~~gnuplot A liquid film drains into a low-pressure meniscus. This mechanism is key in foam ageing via film thinning/Plateau borders thickening.
set term svg enhanced size 750,250 font ",10" 
set multiplot layout 1,3 
set xrange [-1000:500]
set yrange [-10:10]
unset border 
unset xtics
unset ytics
unset key
plot 'log' index 'time = 0.01' u 1:2:(-$2) w filledcu lt rgb "#993498DB", '' index 'time = 0.01' u 1:2 w l lw 2 lt rgb "#3498DB", '' index 'time = 0.01' u 1:(-$2) w l lw 2 lt rgb "#3498DB"
plot 'log' index 'time = 1' u 1:2:(-$2) w filledcu lt rgb "#993498DB", '' index 'time = 1' u 1:2 w l lw 2 lt rgb "#3498DB", '' index 'time = 1' u 1:(-$2) w l lw 2 lt rgb "#3498DB"
plot 'log' index 'time = 10' u 1:2:(-$2) w filledcu lt rgb "#993498DB", '' index 'time = 10' u 1:2 w l lw 2 lt rgb "#3498DB", '' index 'time = 10' u 1:(-$2) w l lw 2 lt rgb "#3498DB"
unset multiplot
~~~

## Flowrate as a function of thickness

[Breward & Howell, 2002](#Breward2002) investigated this problem theoretically in the limit of thin films dominated by viscosity and capillarity. In this limit, the velocity scale $U$ in the film is
$$
U = \frac{\gamma}{\mu}\sqrt{\frac{h_0}{R_c/2}}.
$$
Here $\gamma$ stands for the surface tension, $\mu$ for the dynamic viscosity, $h_0$ for the film thickness and $R_c$ for the curvature radius of the meniscus. From this relation it is readily apparent that the flowrate leaving the liquid film will scale as $h^{3/2}$. Breward & Howell managed to obtain (after an infamous division by $h^{3/2}$ of the equation...) an expression for the non-dimensional flowrate:
$$
Q = \frac{3\sqrt{2} h^{3/2}}{16}
$$
indicating a "quenching" of the drainage as the film thins out. We recover this relation with Basilisk:
~~~gnuplot Relation between (nondimensional) flowrate and film thickness
reset
set xlabel 'Thickness'
set ylabel 'Flowrate'
plot [0:1] 'drainage.thickness' u 2:3 w l t 'Multilayer N = 256 nl = 1', 1*3./8./sqrt(2.)*x**(1.5) t 'Asymptotic theory'
~~~

## Evolution of the film thickness
With the liquid flux out of the liquid lamella, we can track the evolution of the film thickness with time:
$$
h = \left(1 + \frac{3\sqrt{2}}{32} t \right)^{-2}
$$
~~~gnuplot Evolution of the film thickness
set xlabel 'Viscous time'
set ylabel 'Thickness'
plot 'drainage.thickness' u 1:2 w l t 'Multilayer N = 256 nl = 1', 'drainage.thickness' u 1:4 w l t 'Asymptotic theory'
~~~

## Junction profile
The junction between the film and the Plateau border adopts a universal shape, function of the rescaled variable $\xi = \frac{x}{\sqrt{R_c h / 2}}$:
~~~gnuplot Junction profile
set xlabel "xi"
set ylabel "h/hl"
set xrange [-7:3]
set key top left
plot 'drainage.junction' index 'time = 1' u 1:2 t 't = 1', '' index 'time = 2' u 1:2 t 't = 2', '' index 'time = 5' u 1:2 t 't = 5', '' index 'time = 10' u 1:2 t 't = 10', 'drainage.junction.theory' u 1:2 t 'Asymptotic theory' w l dt 2 lw 2 lt rgb "gray50"
~~~

## References

~~~bib
@article{Breward2002,
	author = {C. J. W. Breward and P. D. Howell},
	journal = {J. Fluid Mech.},
	pages = {379--406},
	title = {The drainage of a foam lamella},
	volume = {458},
	year = {2002}
}
~~~

# Basilisk code

We use the `hydro-tension` version of the implicit multilayer solver.
*/
#include "grid/multigrid1D.h"
#include "layered/hydro-tension.h"
#include "layered/implicit.h"
#include "layered/remap.h"

/**
## Geometry 

Typically, in such drainage problems, there is a very large scale separation between the film thickness, its horizontal extent and the curvature radius. Here we consider a domain of size 10000.
Note also that because of the plug-flow velocity profile expected here, we only use a single layer.
*/
#define L 10000.   // The domain size is 10000 as wide as the film is thick
#define LEVEL 8    // horizontal resolution
#define layers 1  // vertical resolution ("mono-layer solver"!)
#define Rc 10000.  // Plateau border curvature radius
#define sqrtRc 70.7106781187
#define fl 8000.   // The film length is 4/5 of the domain
#define Oh 10.     // Non-dimensional viscosity (Ohnesorge number)
#define visctime (fl * Oh * sqrtRc [0,1])
#define tend (10. * visctime)

FILE * fpjunction  = NULL;
FILE * fpthickness = NULL;

/**
## Main function 
*/
int main()
{
  L0 = L;
  origin (-fl);
  N = 1 << LEVEL;
  CFL_H = 3;
  nl = layers;
  // gradient = NULL ; // unnecessary: slopes are kept small thanks to scale separation
  lambda_b[] = {HUGE}; // bottom free slip
  G = 0;
  nu = Oh;

  run();
}

double circular_arc (double pt, double radius) {
  double value = radius + 0.5 - sqrt( sq(radius) - sq(pt) );
  return value;
}

/**
## Asymptotic junction profile
The asymptotic form of the junction profile $h(\xi)$ is given implicitly as:
$$
2\sqrt{h} - \frac{2}{\sqrt{3}}\arctan\left(\frac{1+2\sqrt{h}}{\sqrt{3}}\right) + \frac{1}{3} \log \left(\frac{1-2\sqrt{h}+h}{1+\sqrt{h}+h}\right) = \sqrt{2} \xi.
$$
Here is a helper function that finds the value of $h$ for a given $\xi$ by a relaxation (Newton) method.
*/
double breward_junction (double xi) {
  double value = 2.;
  double residue = 1.;
  double relax = 0.1;

  while (fabs (residue) > 1.e-6) {
    residue = - 1./(value / (pow (value, 3./2.) - 1.)) *
      (2. * sqrt (value) - sqrt (2.) * xi 
      - 2. * atan ((1. + 2. * sqrt (value)) / sqrt (3.)) / sqrt (3.)
      + 1. / 3. * log ((1. - 2. * sqrt (value) + value) / (1. + sqrt (value) + value)));
    value += relax * residue;
  }

  return value;
}

/**
## Boundary conditions and initialization 

To enforce the boundary conditions at the domain boundaries we add this hook for the `half_advection` event:
*/
event half_advection (i++)
{
  foreach_dimension() {
    ha.n[left] = 0.;
    ha.n[right] = 0.;
    hu.n[left] = 0.;
    hu.n[right] = neumann(0.);
  }
}

event init (i = 0)
{
/** 
We here set the boundary conditions. Symmetry on the left; imposed curvature and constant flux ($h u = $cst) on the right.
*/
  eta[left] = neumann(0.);
  u.n[left] = dirichlet(0.);
  u.n[right] = dirichlet (eta[]*u.n[]/(2*eta[] - eta[-1] + 1./Rc * sq(Delta) * pow(1. + sq(eta[] - eta[-1])/sq(L0/N), 3./2.)));
  eta[right] = dirichlet ((3.*eta[] - eta[-1] + 1./Rc * sq(Delta) * pow(1. + sq(eta[] - eta[-1])/sq(L0/N), 3./2.))/2.);
 
/** 
The initial condition consists in a flat film on the left connected to a constant curvature profile on the right. Height and slopes are continuous whereas there is a jump in curvature.
*/
  foreach(){
    double H = (x < 0. ? 0.5 : circular_arc (x, Rc) );
    foreach_layer() {
      h[] = H/nl;
      u.x[] = 0.;
    }
    eta[] = H;
  }
}
/**
We output the film profile at three different times to create the thumbnails at the top of the page.
*/
event thumbnails (t = {0.01*visctime, 1. * visctime, 10. * visctime}) {
  fprintf (stderr, "\n\n# time = %g\n", t/(visctime));
  foreach ()
    fprintf (stderr, "%g %g\n", x, eta[]);
}

/**
We also output the shape of the junction in rescaled coordinates to compare it with the asymptotic prediction.
*/
event junctionscopy (t = {1*visctime, 2. * visctime, 5. * visctime, 10. * visctime}) { 
  static FILE * fpjunction = fopen ("drainage.junction", "w");
  double thickness_lamella = 0.;
  double length_lamella = 0.;
  double xnm1 = -100.;
  double xn = 100.;
  double xnp1;
  double residue = 1;
  thickness_lamella = 2. * interpolate (eta, X0+2.*L/5.); // interpolation mid-lamella
  while (fabs (residue) > 1e-3) {
    xnp1 = xn - (2. * interpolate (eta, xn) - thickness_lamella * breward_junction (0.)) * (xn - xnm1) / (2. * (interpolate (eta, xn) - interpolate (eta, xnm1)));
    residue = 2. * interpolate (eta, xnp1) - thickness_lamella * breward_junction (0.);
    xnm1 = xn;
    xn = xnp1;
  }
  length_lamella = xnp1 - X0;
  fprintf (fpjunction, "\n\n# time = %g\n", t/(visctime));
  foreach ()
    fprintf (fpjunction, "%g %g\n", ((x - X0) - length_lamella) / (sqrtRc * sqrt (thickness_lamella)), 2. * eta[] / thickness_lamella);
}

event junctionasymptotics (t = visctime) {
  FILE * fpjunctionasympt = fopen ("drainage.junction.theory", "w");
  double xi = -7.;
  while (xi < 3.) {
    fprintf (fpjunctionasympt, "%g %g\n", xi, breward_junction (xi));
    xi += 0.1;
  }
  fclose (fpjunctionasympt);
}

/**
## Trouton viscosity
We here enforce brutally the Trouton viscosity by setting a factor 4 in the longitudinal diffusion. Worse, the numerical treatment is implicit.
*/
event viscous_term (i++) {
  foreach()
    foreach_layer()
      foreach_dimension()
        u.x[] += dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
  // fixme: suppress the factor 4 and use `nh.h` instead
#if 0 
  horizontal_diffusion ((scalar*){u}, 4.*nu, dt);
#else
  NITERMIN = 2; // necessary for convergence
  implicit_horizontal_diffusion ((scalar*){u}, 4.*nu, dt);
  NITERMIN = 1; // back to default
#endif
  foreach()
    foreach_layer()
      foreach_dimension()
        u.x[] -= dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
}

/**
## Logs */

event logthickness (i += 500) {
  static FILE * fpthickness = fopen ("drainage.thickness", "w");
  scalar flowrate[];
  double thickness_lamella = 0;
  double err_thickness = 0.;
  double breward_prediction = 0.;
  foreach() {
    flowrate[] = 0.;
    foreach_layer() {
      flowrate[] += 2. * u.x[] * h[] / (fl / visctime);
    } 
  }
  thickness_lamella = 2. * interpolate (eta, X0+2.*L/5.); // interpolation mid-lamella
  breward_prediction = pow((1. + 3. * sqrt(2.) / 32. * t / visctime), -2);
  err_thickness = fabs (thickness_lamella - breward_prediction) / breward_prediction;
  fprintf (fpthickness, "%g %g %g %g %g\n", t / visctime, thickness_lamella, interpolate (flowrate, 0.), breward_prediction, err_thickness);
  fflush (fpthickness);
}
