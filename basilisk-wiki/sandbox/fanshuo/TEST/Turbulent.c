/**
We use multi-layer solver (hydro.h) to simulate a simple $2D$ stationary turbulent chute flow with Prandtl mixing length model.
*/

/**
Taking 2D Navier-Stokes equations ($u$ represtents the real streamwise velocity along $x$, $v$ represents the real transwise velocity along $y$ ), we apply thin layer approximations and the Reynolds decomposition ($u=U+u'$, $v = V+v'$). One can then obtain the Reynolds Averaged Navier-Stokes equations (RANS). When fully developped, it writes
*/
/**
$$
0 = \rho g \alpha + \frac{\partial}{\partial z}\tau,
\text{where } \tau = \rho(\nu\frac{\partial U}{\partial y}-<u'v'>).
$$
*/
/**
The shear stress of turbulent flow $\tau$ then depends on the laminar viscosity $\nu$ and the term $<u'v'>$, the velocity pertublations. To model this contribution, we choose to define a turbulent viscosity $\nu_t$ (viscosité de Boussinesq), such as 
$$
-<u'v'> = \nu_t\frac{\partial U}{\partial y}.
$$
Now the unknow of this modelisation becomes $\nu_t$. Prandtl provided an expression for this turbulent viscosity :
$$
\nu_t = l^2\frac{\partial U}{\partial y}
$$
According to his definition, this turbulent viscotcity results from flow shear ($\partial U/\partial y$) and characteristic length of eddy ($l$) presented in the turbulent flow. This characteristic length $l$, called "mixing length", can adopt a simple form:
$$
l = \kappa y.
$$
It states that the largest size of eddy increases linearly with $y$, with a prefactor $\kappa\approx 0,41$.

Then the shear stress of turbulent flow can be expressed as
$$
\tau = \rho(\nu+l^2\frac{\partial U}{\partial y})\frac{\partial U}{\partial y}=\rho(\nu + \kappa^2 y^2\frac{\partial U}{\partial y})\frac{\partial U}{\partial y}
$$
For the flows above lisse bottom (aerodynamics..), the value of $\nu$ implicates in the thickness of boundary layer. As we will be interested in flows in rver with rough bottom, one neglect commenly the laminar viscocity $\nu$, considering that the bottom roughness is more important than the boudary layer thickness
*/
/**
Now we can derive the turbulent velocity profil $U(y)$ with the Prandtl mixing length model with $\nu=0$. We choose to integrate Eq. from $z$ to $h$ as we are interested in the flow shape far from the rough bottom, this yields (we neglect the shear at top)
$$
\alpha gh(\frac{y}{h}-1) + \kappa^2y^2(\frac{\partial U}{\partial y})^2=0
$$
If we measure $y$ with $d$, the characteristic bottom roughness, and $U$ with $u^*$ whose expression remains to be determined, the above equation can be rewritten with $\tilde{y}=y/d$ and $\tilde{u}=U/u^*$
$$
\alpha gh(\frac{d}{h}\tilde{y}-1)+\kappa^2\tilde{y}^2 (u^*)^2(\frac{\partial \tilde{u}}{\partial \tilde{y}})^2=0
$$
Considering $d/h<<1$, then according the dominant balance ("principe de moindre dégénérescence"), one gets the characteristic velocity
$$
u^* = \sqrt{\alpha g h},
$$
and the following non-dimension equation is obtained
$$
\kappa^2\tilde{y}^2(\frac{\partial \tilde{u}}{\partial \tilde{y}})^2=1.
$$
Solve it and one can get the logarithmic velocity profile with dimensions
$$
U(y) = \frac{u^*}{\kappa}log(\frac{y}{d})+B
$$
Experimental measurements show the parameter $B$ is a function of $\nu, u^*,d$. When the bottom is very rough, one gets $B=8,5$. With approximation $8,5\approx log(30)/\kappa$ with $\kappa\approx 0,41$, one may get
$$
\frac{U(y)}{u^*} = \frac{1}{\kappa}log(\frac{y}{y_0}), \textit{with }y_0=d/30
$$
Here the $y_0$ represents the vitual bed level where $U(y_0)=0$. With coordinate shift $\tilde{y} = y-y_0$, one can get, where $U(0)=0$
$$
\frac{U(y)}{u^*} = \frac{1}{\kappa}log(1+\frac{y}{y_0})
$$
To be noticed, in order to obtain the above velocity profile, the mixing length $l$ should be according changed, due to the coordinate shift. To do so, we substitute the expression of $U$ into $\tau$ and solve $l$ with Eq, one can get
$$
l = \kappa (y_0+y)\sqrt{1-\frac{y}{h}}
$$
The objective of this test is to obtain the above velocity profile with this modified mixing length, using muli-layer solver. 
*/



#define ML 1
#define HYDRO 1
#define MUI 1
#define TURBULENT  1

#include "grid/multigrid1D.h"
#if !ML
# include "saint-venant.h"
#else // ML
#include "hydroMT.h"

# define phi q
# if !HYDRO
#   include "layered/nh.h"
# endif
# include "layered/remap.h"
#endif // ML


const double NU = 0.1, T0 = 10000, HR = 1.;
scalar uold[];

double Ut(Point point)
{
  double zc = 0.;
  for (int l = - point.l; l < nl - point.l; l++) {
    if (l < 0)
      zc += h[0,0,l];
  }
  zc += h[]/2.;

 
  return sqrt(HR*G*slope)/kappa*log(1+zc/rugo);
}

double Ut2(double z)
{
  
  return sqrt(HR*G*sin(slope))/kappa*log(1+z/rugo);
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - */
int main()
{
  periodic(right);
 
  slope = 0.2618;  // ~ 15 degrees
  L0 = 1;
  G = 1.;
  N = 32; 
  
  rugo = 0.01;
  nu = NU;
  nl = 256; 


  run();
}

/**
We initialise the topography and the initial thickness of each layer *h*. */

event init (i = 0)
{
  foreach() {
    zb[] = 0.;
#if !ML
    h[] = HR - zb[];
#else
    foreach_layer()
      h[] = (HR)/nl;
#endif
 	eta[] = HR;  
  }

//Mettre une vitesse initiale dans le domaine 
  foreach () {
    foreach_layer()
      u.x[] = Ut(point);
    uold[] = 0;
  }
}

event acc(i++){
  foreach () 
    foreach_layer()
      u.x[] = u.x[] + G*sin(slope)*dt;
}


// /**exit criteria  */
event logfile (t += 1; t<=T0) {

  double du = change (u.x,uold);
    if (i > 0 && du < 1e-6)
      return 1; /* stop */
}
/**
## Outputs

We save profiles at regular intervals. */
#if 0
event profiles (t += 5;t<=T0)
{
  foreach (serial) {
#if !ML
    double H = h[];
#else
    double H = 0.;
    foreach_layer()
      H += h[];
#endif
    fprintf (stderr, "%g %g %g\n", x, zb[] + H, zb[]);
  }
  fprintf (stderr, "\n\n");
}
#endif

/**
For the hydrostatic case, we compute a diagnostic vertical velocity
field `w`. Note that this needs to be done within this event because
it relies on the fluxes `hu` and face heights `hf`, which are only
defined temporarily in the [multilayer solver](hydro.h#update_eta). */

#if HYDRO
scalar w = {-1};

event update_eta (i++)
{
  if (w.i < 0)
    w = new scalar[nl];
  vertical_velocity (w, hu, hf);

  /**
  The layer interface values are averaged at the center of each
  layer. */
  
  foreach() {
    double wm = 0.;
    foreach_layer() {
      double w1 = w[];
      w[] = (w1 + wm)/2.;
      wm = w1;
    }
  }
}
#endif // HYDRO

event output (t = end)
{
 
    foreach () {
      double z = zb[];
    foreach_layer() {
      z += h[];
      fprintf (stdout, "%g %g %g %g\n", x, z, u.x[],Ut2(z));
    }
      fprintf (stdout,"\n \n");   
  }
#if HYDRO
  delete ({w});
#endif
 
}



/**

~~~gnuplot Velocity and shear profiles for Turbulent flow
 set logscale x
 set xlabel "y"
 set ylabel "u"
 p "out" u 2:3 w p t'U computed', "" u 2:4 w l t'U ana', 
~~~
*/ 























