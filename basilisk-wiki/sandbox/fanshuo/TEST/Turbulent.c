/**
# Purpose

We use *modifed* multi-layer solver (hydroMT.h *?*) to simulate a simple $2D$ stationary turbulent chute flow with Prandtl mixing length model.
*/

/**
# Turbulent flow, the Prandtl mixing length model

Taking 2D Navier-Stokes equations ($u$ represtents the real streamwise velocity along $x$, $v$ represents the real transwise velocity along $y$ ), we apply the Reynolds decomposition ($u=U+u'$, $v = V+v'$, $<>$ repressents a mean value, $<u>=U$, $<u'>=0$ etc). One can then obtain the Reynolds Averaged Navier-Stokes equations (RANS). When fully developped, and using thin layer approximations,   it writes
 
$$
\rho (\frac{\partial U }{\partial t}+
U \frac{\partial U }{\partial x}+
W \frac{\partial U}{\partial z}) = 
-\rho g \frac{\partial h }{\partial x}+
\rho g \alpha + \frac{\partial}{\partial z}\tau,
\text{where } \tau = \rho(\nu\frac{\partial U}{\partial y}-<u'v'>).
$$
 
The *total* shear stress of turbulent flow $\tau$ then depends on the laminar viscosity $\nu$ and the term $<u'v'>$, the *mean product* of velocity pertublations. To model this contribution, we choose to define a turbulent viscosity $\nu_t$ (Boussinesq viscosity), such as 
$$
-<u'v'> = \nu_t\frac{\partial U}{\partial y}.
$$
Now the unknow of this modelisation becomes $\nu_t$. Prandtl $proposed$ an expression for this turbulent viscosity :
$$
\nu_t = l^2\frac{\partial U}{\partial y}
$$
According to his definition, this turbulent viscosity results from flow shear ($\partial U/\partial y$) and characteristic length of eddy ($l$) presented in the turbulent flow. This characteristic length $l$, called "mixing length", can adopt a simple form:
$$
l = \kappa y.
$$
It states that the largest size of eddy increases linearly with $y$, with a prefactor $\kappa\approx 0,41$.

Then the  stress of turbulent flow can be expressed as
$$
\tau = \rho(\nu+l^2\frac{\partial U}{\partial y})\frac{\partial U}{\partial y}=\rho(\nu + \kappa^2 y^2\frac{\partial U}{\partial y})\frac{\partial U}{\partial y}
$$
For flows above *smooth* bottom (aerodynamics..), the value of $\nu$ implicates *a "viscous" sublayer* where $U/u^* \simeq u^* y/\nu$. As we will be interested in flows in river with rough bottom, one neglects commonly the laminar viscocity $\nu$, considering that the bottom roughness is more important than the *viscous sublayer* thickness (*which is* $O(\nu/u^*$))
*/
/**
Anyway we can derive the turbulent velocity profile $U(y)$ with the Prandtl mixing length model with $y \gg \nu/u^*$ and $y \ll h$. We choose to integrate Eq. from $z$ to $h$ as we are interested in the flow shape far from the rough bottom, this yields (we neglect the shear at top)
$$
\alpha gh(\frac{y}{h}-1) + \kappa^2y^2(\frac{\partial U}{\partial y})^2=0
$$
If we measure $y$ with $d$ (this $d$d is any length, smaller than $h$,
classicaly one takes $\nu/u^*$, but we will take at the end the characteristic bottom roughness), and $U$ with $u^*$ whose expression remains to be determined, the above equation can be rewritten with $\tilde{y}=y/d$ and $\tilde{u}=U/u^*$
$$
\alpha gh(\frac{d}{h}\tilde{y}-1)+\kappa^2\tilde{y}^2 (u^*)^2(\frac{\partial \tilde{u}}{\partial \tilde{y}})^2=0
$$
Considering $d/h \ll 1$, then according the dominant balance ("principe de moindre dégénérescence"), one gets the characteristic velocity (what ever is the chosen length))
$$
u^* = \sqrt{\alpha g h},
$$
and the following non-dimension equation is obtained
$$
\kappa^2\tilde{y}^2(\frac{\partial \tilde{u}}{\partial \tilde{y}})^2 \simeq 1.
$$
Solve it and one can get the logarithmic velocity profile with dimensions
$$
U(y) = \frac{u^*}{\kappa}log(\frac{y}{d})+B
$$
Experimental measurements show the parameter $B$ is a function of $\nu, u^*,d$. When the bottom is very rough, one gets $B=8,5$, and indeed use $d$ the saize of grain is used. With approximation $8,5\approx log(30)/\kappa$ with $\kappa\approx 0,41$, one may get
$$
\frac{U(y)}{u^*} = \frac{1}{\kappa}log(\frac{y}{y_0}), \textit{with }y_0=d/30
$$
Here the $y_0$ represents the virtual bed level where $U(y_0)=0$. With coordinate shift $\tilde{y} = y-y_0$, one can get, where $U(0)=0$
$$
\frac{U(y)}{u^*} = \frac{1}{\kappa}log(1+\frac{y}{y_0})
$$
To be noticed, *this profile is "given" in the litterature in order to obtain the above velocity profile*, the mixing length $l$ should be according changed, due to the coordinate shift. To do so, we substitute the expression of $U$ into $\tau$ and solve $l$ with Eq, one can get a mixing length *so that the profile is exactly the logaritmic one*. 
$$
l = \kappa (y_0+y)\sqrt{1-\frac{y}{h}}
$$
The objective of this test is to obtain the above velocity profile with this modified mixing length, using modified muli-layer solver. 



# Code
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

/** exact solution, without dimension (for layers)
*/
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
/** exact solution, without dimension as a function
*/
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

/**
Simple acceleration : projection of weight along the slope
*/

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
# Results

Plots of velicity profile in cartesian and semilog 

~~~gnuplot Velocity and shear profiles for Turbulent flow
 set xlabel "u"
 set ylabel "y"
 p "out" u 3:2 w p t'U computed','' u (1.272*log($2/0.01+1)):2 t'u* /kappa*log(1+y/y0)' w l
~~~


~~~gnuplot Velocity  r profiles for Turbulent flow in semi log
 set logscale x
 set xlabel "y"
 set ylabel "u"
 p "out" u 2:3 w p t'U computed',1.272*log(x/0.01+1) t'u* /kappa*log(1+y/y0)'
~~~
*/ 























