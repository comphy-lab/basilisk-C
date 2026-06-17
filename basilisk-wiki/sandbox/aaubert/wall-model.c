/**
# Turbulent flow in a channel with a "wall model"

This is a variant of this [example](loglayer.c) which solves the
entire turbulent flow profile in a channel with a viscous sublayer and
Prandtl's "mixing length" turbulence model.

The goal is to exclude the small scales and strong gradients
associated with the viscous sublayer and its transition toward the log
profile, to drastically reduce the necessary mesh refinement.

The configuration is the same but the left boundary of the channel is
reduced by a small amount $y_c$ on which a special boundary condition
is applied: a "wall model". */

#include "grid/bitree.h"
#include "run.h"
#include "diffusion.h"

/**
The parameters are the acceleration $G$, the fluid viscosity $\nu$ and
the width of the channel $W$ (including the viscous sublayer etc.). */

const double G = 1., nu = 1e-4, W = 1. [1];
const double kappa = 0.41 [0]; // von Karman constant

/**
The left boundary is set at $y_c$, a small value compared to $W$ but
sufficient to exclude the viscous sublayer, since $y_+ =
y_c\nu/u_\star = 10$ here. */

double yc = 1e-3;

/**
The numerical parameters are the minimum and maximum level of
refinement. */

int minlevel = 3, maxlevel = 11;

/**
## Wall model

To construct the left boundary condition at $y = y_c$, we start from
the full analytical solution (see [Lagrée, 2026](#lagree2026) p.34)
$$
u_+ = \frac{\frac{1}{y_+} - \frac{\sqrt{4y^2_+\kappa^2+1}}{y_+} + 
      2\kappa\sinh^{-1}(2y_+\kappa)}{2\kappa^2}
$$
*/

double u_plus (double y_plus)
{
  return (1./y_plus - sqrt(4.*sq(y_plus*kappa) + 1.)/y_plus +
          2.*kappa*asinh(2.*y_plus*kappa))/(2.*sq(kappa));
}

/**
The gradient of this function is
$$
\partial_{y_+} u_+ = \frac{\sqrt{4\kappa^2y_+^2 + 1} - 1}{2\kappa^2y_+^2}
$$
*/

double du_plusdy_plus (double y_plus) {
  return (sqrt(4.*sq(kappa*y_plus)+1.)-1.)/(2.*sq(kappa*y_plus));
}

/**
This two formulas give two ways to obtain the friction velcoity $u_\star$. It can either be obtain from the velocity at a point $y_0$. This requires to solve the non-linear equation $u(y_0)=u_\star u_+\left(\frac{y_0u_\star}{\nu}\right)$. This equation is solve with Newton's method.*/

double f_newton(double y0, double uy, double utau) {
  return uy-utau*u_plus(y0*utau/nu);
}

double f_newton_derivative(double y0, double uy, double utau) {
  double yplus=y0*utau/nu;
  return -1./(2.*sq(kappa))*(1./yplus-sqrt(4.*sq(kappa)*sq(yplus)+1.)/yplus+2.*kappa*asinh(2.*yplus*kappa))-utau*y0/nu*(sqrt(4.*sq(kappa)*sq(yplus)+1.)-1.)/(2.*sq(kappa)*sq(yplus));
}

void ustar_compute_velocity(double y0, double uy, double * u_star) {
  if (fabs(uy)>0.) {
    double utau=*u_star?*u_star:1.;
    double utauprec=utau+1.;
    while (fabs(utau-utauprec)/utau>1.e-5) {
      utauprec=utau;
      utau=utauprec-f_newton(y0,fabs(uy),utauprec)/f_newton_derivative(y0,fabs(uy),utauprec);
    }
    *u_star=utau;
    fprintf(stderr,"%g\n",utau);
  }
  else {
    *u_star=0.;
  }
}

/**
Or the friction velcoity can come from the derivative at a point $y_0$. Assuming that $U$ follows the analytical solution above. By definition we have
$$
\partial_y U|_{y_0} = \frac{u_\star}{y_\star} \partial_{y_+} u_+|_{y_0}
$$
with $u_\star$ the "friction velocity", $y_\star = \nu/u_\star$ and
$y_+ = y/y_\star$. Using the derivative above we can write
$$
\partial_y U|_{y_0} = \frac{u^2_\star}{\nu}  \frac{\sqrt{4 \kappa^2 
  (y_0 u_\star / \nu)^2 + 1} - 1}{2 \kappa^2  (y_0 u_\star / \nu)^2}
$$
which gives after some manipulations
$$
  u_\star = \frac{\nu}{2 \kappa y_0}  \sqrt{\left( 1 + \frac{2 \kappa^2
  y_c^2}{\nu} \partial_y U|_{y_0} \right)^2 - 1}
$$*/
void ustar_compute_derivative(double y0, double dudy, double * u_star) {
  
*u_star = nu/(2.*kappa*y0)*sqrt(sq(1. + 2.*sq(kappa*y0)/nu*fabs(dudy)) - 1.);
}

/**
We are looking for the "correct" boundary condition that need to be applied at the $y=y_c$ to account for the fact that the mesh is unresolved here. We present here three methods.

The first one is based on a Dirichlet boundary condition with the velocity at $y=y_c$ given by the formula $u(y_c)=u_\star u_+(y_c^+)$. The friction velocity is obtained thanks to the derivative at $y=y_c$. It will be named "Dirichlet" in the following.

The second one is based on a Neumann boundary condition with the normal derivative linked to the formula for the derivative at $y=y_c$ : $\frac{\partial u}{\partial n}=-\frac{u_\star^2}{\nu}\left.\frac{du^+}{\partial y^+}\right|_{y_c}$. In this method, the friction velocity is obtained from the velocity at the first point of the mesh. It will be named "Neumann" in the following.

The last one is based on a Navier boundary condition. We want to construct a slip length $\lambda$ relating $U$ and its gradient
$\partial_yU$ at location $y_c$ i.e.
$$
U|_{y_c} = \lambda\partial_yU|_{y_c}
$$
*/

double slip (double yc, double dudy, double utau)
{
  /**
  Knowing $u_\star$ we can now compute $y_{c+} = y_c/y_\star$ and
  $$
  \lambda = \frac{u_\star u_+ (y_{c +})}{\partial_y U|_{y_c}}
  $$
  */
  
  double yc_plus = yc*(utau)/nu;
  return fabs(dudy) > 0. ? (utau)*u_plus (yc_plus)/fabs(dudy) : 0.;
}
/**
   This method will be named "Navier" in the following. */


bool navier_bc=false;
bool neumann_bc=false;
bool dirichlet_bc=false;

static FILE * log_dirichlet;
static FILE * log_neumann;
static FILE * log_navier;

static FILE * velocity_dirichlet;
static FILE * velocity_neumann;
static FILE * velocity_navier;

/**
## Setup

The timestep cannot be too large. The tolerance controls the accuracy
of the implicit diffusion solution. The initial mesh has
$2^\text{minlevel}$ grid points. */

int main() {
  DT = 0.1 [0,1];
  N = 1 << minlevel;
  
  log_dirichlet=fopen("log_dirichlet","w");
  log_neumann=fopen("log_neumann","w");
  log_navier=fopen("log_navier","w");

  velocity_dirichlet=fopen("velocity_dirichlet","w");
  velocity_neumann=fopen("velocity_neumann","w");
  velocity_navier=fopen("velocity_navier","w");

  
  for (maxlevel=4;maxlevel<10;maxlevel++) {
    for (int j=0;j<3;j++) {

      dirichlet_bc=(j<1);
      neumann_bc=(j%2==1);
      navier_bc=(j>1);
      
      /**
	 The domain size is $W - y_c$ and the left boundary is at $x =
	 y_c$. */
  
      size (W - yc);
      origin (yc);
      run();
      
    }
  }
}

/**
## Boundary conditions

The only unknown is the velocity $U$. On the right boundary the
default Neumann zero condition is used. On the left boundary we apply different boundary condition depending on the model used. */

scalar U[];
double lambda = 0;
double u_star=1.;

event init(i=0) {
  
  foreach() {
    U[]=1.;   //if u_star is found from the value at yc, U should not be zero
  }
  
  if (dirichlet_bc)
    U[left]=dirichlet(u_star*u_plus(yc*u_star/nu));
  if (neumann_bc)
    U[left]=neumann(-sq(u_star)/nu*du_plusdy_plus(yc*u_star/nu));
  if (navier_bc)
    U[left] = navier (0., lambda);

}

/**
## Time integration 

In this event we solve the diffusion equation
$$
\partial_t U = \partial_x(D\partial_x U) + G
$$
*/

mgstats mg;

event integration (i++)
{

  /**
  The timestep is computed taking account events if necessary and with
  a maximum value bounded by `DT`. */
  
  dt = dtnext (DT);

  /**
  Field `a` is just the acceleration. */
  
  scalar a[];
  foreach()
    a[] = G;

  /**
  We update the friction velocity $u_\star$ and the slip length
  $\lambda$ using only the value of the velocity gradient on the left
  boundary. */
  
  foreach_boundary (left) {
    //We first compute $u^*$ with the correct method for each boundary condition
    if (neumann_bc)
      ustar_compute_velocity(yc+Delta/2.,U[],&u_star);
    else
      ustar_compute_derivative(yc,(U[]-U[-1])/Delta,&u_star);

    
    if (navier_bc)
      lambda = slip (yc, (U[] - U[-1])/Delta, u_star);
  }

  /**
  To initialize the diffusion coefficient $D$, we need a model for the
  turbulent viscosity $\nu_t$. We use Prandtl's "mixing length" model
  which assumes that an important characteristic length of the flow is
  the distance to the boundary. The mixing length $l$ is then defined
  as $\kappa x$ and the turbulent viscosity is
  $$
  \nu_t = l^2 \partial_x U
  $$
  The diffusion coefficient is the sum of the fluid and the turbulent
  viscosity. */
  
  face vector D[];
  foreach_face(x) {
    double l = kappa*x;
    double nu_t = sq(l)*fabs (U[] - U[-1])/Delta;
    D.x[] = nu + nu_t;
  }
  mg = diffusion (U, dt, D, r = a);
}

/**
## Mesh adaptation

We do not adapt at every timestep to minimize adaptation noise. */

event adapt (i += 10)
  adapt_wavelet ({U}, {1e-2}, maxlevel, minlevel);  

/**
## Convergence

1000 timesteps are necessary/sufficient to reach a stationary
solution. */

event logfile (i++; i <= 1000) {
  if (i == 0)
    printf ("\n# yc = %g\n", yc);
  printf ("%g %d %g %g %g\n", t, mg.i, perf.t, lambda, u_star);
}

/**
## Results

We output the final velocity profile. The (analytical) characteristic
"friction velocity" is $u_\star = \sqrt{GW}$. */

event profile (t = end) {
  double u_star = sqrt(G*W), y_star = nu/u_star;
  FILE * output_v;
  if (neumann_bc)
    output_v=velocity_neumann;
  if (dirichlet_bc)
    output_v=velocity_dirichlet;
  if (navier_bc)
    output_v=velocity_navier;
  
  fprintf (output_v,"\n");
  fprintf (output_v, "\n# yc = %g level = %d\n", yc,maxlevel);
  foreach (serial)
    fprintf (output_v, "%g %g %g %d\n", x/y_star, U[]/u_star, Delta/y_star, level);
  
  FILE * output_err;
  if (neumann_bc)
    output_err=log_neumann;
  if (dirichlet_bc)
    output_err=log_dirichlet;
  if (navier_bc)
    output_err=log_navier;
  
  double err=0.;
  double n_point=0.;
  foreach (serial)
    if (x*u_star/nu<800) {  //only compare the value for $y^+<800 to be in the log zone
      err+=fabs(U[]-u_star*u_plus(x*u_star/nu));
      n_point+=1.;
    }
  if (n_point)
    fprintf(output_err,"%d %g\n",maxlevel,err/n_point);
}

/**
The solutions with $y_c = 10^{-3}$ corresponding to $y_{c+} = 10$ agree closely with the analytical solution for the highest refinment level and all the models.

~~~gnuplot Numerical and analytical velocity profiles for level 8
set xlabel 'y_+'
set ylabel 'u_+'
kappa = 0.41
set logscale x
u_plus(y) = (1./y - sqrt(4.*(y*kappa)**2 + 1.)/y + 2.*kappa*asinh(2.*y*kappa))/(2.*kappa**2)
u_log(y) = (log(y) + log(kappa) + log(4.) - 1.)/kappa
set key bottom right
set grid
plot [1:][0:21]\
'velocity_dirichlet' index 'yc = 0.001 level = 8' u 1:2 pt 5 ps 1.05 t 'level=8 Dirichlet', \
'velocity_navier' index 'yc = 0.001 level = 8' u 1:2 pt 5 t 'level=8 Navier', \
     'velocity_neumann' index 'yc = 0.001 level = 8' u 1:2 pt 5 t 'level=8 Neumann', \
     u_plus(x) lt 2 t 'Full analytical solution', \
     x lt 3 t 'Viscous sublayer', \
     u_log(x) lt 4 t 'Log layer'
~~~


But some difference appears at lower mesh refinment for the same $y_c$.

~~~gnuplot Numerical and analytical velocity profiles for level 5
set xlabel 'y_+'
set ylabel 'u_+'
kappa = 0.41
set logscale x
u_plus(y) = (1./y - sqrt(4.*(y*kappa)**2 + 1.)/y + 2.*kappa*asinh(2.*y*kappa))/(2.*kappa**2)
u_log(y) = (log(y) + log(kappa) + log(4.) - 1.)/kappa
set key top left
set grid
plot [1:][0:61]\
'velocity_dirichlet' index 'yc = 0.001 level = 5' u 1:2 pt 5 ps 1.05 t 'level=5 Dirichlet', \
     'velocity_navier' index 'yc = 0.001 level = 5' u 1:2 pt 5 t 'level=5 Navier', \
     'velocity_neumann' index 'yc = 0.001 level = 5' u 1:2 pt 5 t 'level=5 Neumann', \
     u_plus(x) lt 2 t 'Full analytical solution', \
     u_log(x) lt 4 t 'Log layer'
~~~


To highlight this difference between the models, we can plot the rate of convergence for this $y_c$ for the three models. It shows a higher rate of convergence for the Navier and Dirichlet models but the error with the Neumann model is always smaller. Some oscillations appear for the Dirichlet model at high refinment.

~~~gnuplot Rate of convergence
set xlabel 'Level'
set ylabel 'error'
unset logscale
set logscale y
set key bottom right
plot 'log_dirichlet' u 1:2 w l lw 1.1 t 'Dirichlet', \
'log_navier' u 1:2 w l t 'Navier', \
'log_neumann' u 1:2 w l t 'Neumann'
~~~

## See also

* [Viscous sublayer and log layer for turbulent flow in a channel](loglayer.c)

## References

~~~bib
@misc{lagree2026,
  author = {Pierre-Yves Lagrée},
  title = {{Equations de Saint Venant et application aux mouvements de fonds érodables. 
            ``MU4MEF04 - Ondes et Ecoulements en milieu naturel", M1 SU}},
  year = 2026,
  url = {http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf}
}
~~~
*/