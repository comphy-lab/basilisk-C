#include "run.h"
#include "timestep.h"

/**
We perform a 4th order Runge-Kutta integration for a double using Basilisk events to integrate a system 
$$ \frac{d\,Y}{dt} = f(t,Y) $$
where $Y= (y, \dot{y})^T$ is the state vector.

*/

#include "RK4.h" // one step integration

/**
Some non-dimensionnal parameters
*/

// double pinf=1.;   // pressure at infinity, used to adimensionalize
// double rhoL=1.;   // liquid density, used to adimensionalize
// double R0=1.;     // initial bubble radius, used to adimensionalize

double pg0=1./6.;    // initial gas pressure    !! pg0 = pl0 + 2/We
double dotR0=0.;     // initial bubble velocity
double We=1e6;       // Weber number
double Re=1e4;       // Reynolds number
double Ma=0.007;     // Mach number

double gammaG=1.4;   // gas polytropic coeficient


/**
We define the function to integrate the Keller-Miksis, cf [Keller and Miksis (1980) (eq.(3.4) and their is a typo in it)] or [prosperetti (1984) (eq.(15))] or [gaitan et al (1992) (eq.(1))] or [Lauteborn and Kurz (2010) (eq.(11))]) for a bubble initially at equilibrium subjected to a suddent high pressure increase at infinity:
$$ \left( 1 - \dfrac{\dot{R}}{c}\right) R \ddot{R} + \dfrac{3}{2} (\dot{R})^2 \left(1-\dfrac{\dot{R}}{3c}\right) = \left( 1 + \dfrac{\dot{R}}{c}\right) \dfrac{p_l - p_\infty}{\rho_l} + \dfrac{R}{\rho_l c}\dfrac{\text{d} p_l}{\text{d}t} $$
where we dropped the retarded time term (usefull for sinusoidal infinity pressure), and with $\rho_l$ the liquid density, $c$ the sound celerity in liquid, $R$ the bubble radius, $p_\infty$ the pressure at infinity and $p_l$ the liquid pressure at the bubble interface:
$$ p_l = p_g - 2\dfrac{\sigma}{R} - 4\mu\dfrac{\dot{R}}{R} $$
with $\sigma$ the surface tension, $\mu$ the viscosity, $p_g$ the gas pressure which can be expressed as
$$ p_g = p_{g,0} \left(\dfrac{R_0}{R}\right)^{3k} $$
where we dropped the saturated vapor pressure and
with $p_{g,0}$ the initial gas pressure, $R_0$ the initial bubble radius and $k$ the polytropic coefficient.
The pressure time derivative is thus
$$ \dfrac{\text{d} p_l}{\text{d}t} = -3k p_{g,0} \left(\dfrac{R_0}{R}\right)^{3k} \dfrac{\dot{R}}{R} + \dfrac{2\sigma\dot{R}}{R^2} - \dfrac{4\mu\ddot{R}}{R} + \dfrac{4\mu(\dot{R})^2}{R^2} $$

Using $\rho_l$, $p_\infty$, $R_0$ to create non-dimensional variables, we obtain 
$$ \left( 1 - M_a\dot{R}^*\right) R^* \ddot{R}^* + \dfrac{3}{2} (\dot{R}^*)^2 \left(1-\dfrac{M_a}{3}\dot{R}^*\right) = \left( 1 + M_a\dot{R}^*\right) (p_l^* - 1) + M_aR^*\dfrac{\text{d} p_l^*}{\text{d}t^*} $$
and
$$ p_l^* = p_{g,0}^* \left(\dfrac{1}{R^*}\right)^{3k} - \dfrac{2}{W_e} \dfrac{1}{R^*} - \dfrac{4}{R_e} \dfrac{\dot{R}^*}{R^*} $$
and
$$ \dfrac{\text{d} p_l^*}{\text{d}t^*} = -3k p_{g,0}^* \left(\dfrac{1}{R^*}\right)^{3k} \dfrac{\dot{R}^*}{R^*} + \dfrac{2}{W_e} \dfrac{\dot{R}^*}{(R^*)^2} - \dfrac{4}{R_e} \dfrac{\ddot{R}^*}{R^*} + \dfrac{4}{R_e}\dfrac{(\dot{R}^*)^2}{(R^*)^2} $$
with $W_e = \rho_l u_c^2 R_0/\sigma$, $R_e = \rho_l u_c R_0/\mu_l$, $M_a = u_c/c = 1/c^*$ and $u_c = \sqrt{p_\infty/\rho_l}$.

As the pressure time derivative contains $\ddot{R}$, we rewrite the equation as
$$ \left( 1 - M_a\dot{R}^*\right) R^* \ddot{R}^* + 4\dfrac{M_a}{R_e}\ddot{R}^* + \dfrac{3}{2} (\dot{R}^*)^2 \left(1-\dfrac{M_a}{3}\dot{R}^*\right) = \left( 1 + M_a\dot{R}^*\right) (p_l^* - 1) + M_aR^*S $$
with 
$$ S = -3k p_{g,0}^* \left(\dfrac{1}{R^*}\right)^{3k} \dfrac{\dot{R}^*}{R^*} + \dfrac{2}{W_e} \dfrac{\dot{R}^*}{(R^*)^2} + \dfrac{4}{R_e}\dfrac{(\dot{R}^*)^2}{(R^*)^2} $$


For time integration, the equation rewrites as
$$ \ddot{R}^* =\dfrac{ - \dfrac{3}{2} (\dot{R}^*)^2 \left(1-\dfrac{M_a}{3}\dot{R}^*\right) + \left( 1 + M_a\dot{R}^*\right) (p_l^* - 1) + M_aR^*S  }{\left( 1 - M_a\dot{R}^*\right) R^* + 4\dfrac{M_a}{R_e} } $$

*/

my_vec f_rk4(double t, my_vec Y){   // dY/dt = f_rk4(t,Y)

  double pl = pg0*pow(Y.y0,-3.*gammaG) - 2./We/Y.y0 - 4./Re*Y.y1/Y.y0;
  double RS = -3.*gammaG*pg0*pow(Y.y0,-3.*gammaG)*Y.y1 + 2./We*Y.y1/Y.y0 + 4./Re*sq(Y.y1)/Y.y0;  // R*S
  
  my_vec Y2;
  Y2.y0 = Y.y1;
  Y2.y1 = ( -3./2.*sq(Y.y1)*(1.-Ma/3.*Y.y1) + (1.+Ma*Y.y1)*(pl-1.) + Ma*RS )/( (1.-Ma*Y.y1)*Y.y0 + 4.*Ma/Re  ); // K-M eq.
  return Y2;
}

/**
We initialize the state vector.
It contains the (non-dimensional) bubble radius and its first derivative.

*/

my_vec Y; 

event init (i=0) {
  Y.y0 = 1.; // R, bubble radius
  Y.y1 = dotR0; // dot(R)
}

/**
We integrate at each time step
*/

event integ(i++,last){
  Y = RK4(t, dt, Y, f_rk4);  // 4th order Runge-Kutta
}


/**
We include and define some useful things, such as the timestep and the output of the result.
*/

double dtmax;
event set_dtmax (i++) dtmax = 0.001;

event stability (i++) {
  dt = dtnext(dtmax);   // dtnext() updates the internal variables of Basilisk to compute each iterations
}

int main(){
  init_grid(1<<1);
  run();
  free_grid();
}


/**
Output: non-dimensional time $t^* = t/t_c$ with $t_c = \sqrt{p_\infty/\rho_l}/R_0$, $R^*(t)$, $\dot{R}^*(t)$, and 
 $$p_g^*(t) = p_{g,0}^* \left(\frac{1}{R^*}\right)^{3k}$$
*/

event logfilef(i++){
  fprintf(stderr,"%g %g %g %g\n", t, Y.y0, Y.y1, pg0*pow(1./Y.y0,3.*gammaG) );
}

event endsimu(t=9){}



/**
We plot the result

~~~gnuplot Bubble radius (non-dimensional) evolution predicted by the Rayleigh-Plesset equation and RK4 integration.
set term pngcairo enhanced size 500,500
set output 'radius.png'

set xlabel 't^*'
set ylabel 'R^*'

p "log" u 1:2 w p t 'RK4'
~~~

~~~gnuplot Bubble pressure (non-dimensional) evolution using the Rayleigh-Plesset equation and RK4 integration.
set term pngcairo enhanced size 500,500
set output 'pressure.png'

set xlabel 't^*'
set ylabel 'p^*'

p "log" u 1:4 w p t 'RK4'
~~~
*/



