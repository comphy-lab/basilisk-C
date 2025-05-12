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

double pv=0.;        // vapor pressure
double pg0=1./3.;    // initial gas pressure    !! pg0 = pl0 + 2/We
double dotR0=0.;     // initial bubble velocity
double We=1e1;       // Weber number
double Re=1e4;       // Reynolds number

double gammaG=1.4;   // gas polytropic coeficient


/**
We define the function to integrate the Rayleigh-Plesset equation, cf [Brennen (2013), Cavitation and bubble dynamics, eq. (2.27), paragraph 2.4] for a bubble initially at equilibrium subjected to a suddent high pressure increase at infinity:
$$ \dfrac{p_v-p_\infty}{\rho_l} + \dfrac{p_{g,0}}{\rho_l}\left(\dfrac{R_0}{R}\right)^{3k} = R\ddot{R} + \dfrac{3}{2}\left(\dot{R}\right)^2 + 4\nu_l\dfrac{\dot{R}}{R} + \dfrac{2\sigma}{\rho_l R} $$ 
with $p_v$ the saturated vapor pressure, $p_\infty$ the pressure at infinity, $\rho_l$ the bulk (liquid) density, $p_{g,0}$ the initial (gas) bubble pressure, $R = R(t)$ the bubble radius and $R_0$ the initial bubble radius, $k$ the polytropic coefficient ($k=1$: constant bubble temperature, $k=\gamma$: adiabatic), $\nu_l$ the dynamic viscosity and $\sigma$ the surface tension.

Using $\rho_l$, $p_\infty$, $R_0$ to create non-dimensional variables, we obtain
$$ p_v^*-1 + p_{g,0}^*\left(\dfrac{1}{R^*}\right)^{3k} = R^*\ddot{R}^* + \dfrac{3}{2}\left(\dot{R}^*\right)^2 + \dfrac{4}{R_e}\dfrac{\dot{R}^*}{R^*} + \dfrac{2}{W_e}\dfrac{1}{R^*} $$
which gives
$$ \ddot{R}^* = \dfrac{p_v^*-1}{R^*} + \dfrac{p_{g,0}^*}{R^*}\left(\dfrac{1}{R^*}\right)^{3k} - \dfrac{3}{2}\dfrac{\left(\dot{R}^*\right)^2}{R^*} - \dfrac{4}{R_e}\dfrac{\dot{R}^*}{\left(R^*\right)^2} - \dfrac{2}{W_e}\dfrac{1}{\left(R^*\right)^2} $$
with $p_{g,0}^* = p_{g,0}/p_\infty$, $p_v^* = p_v/p_\infty$, $R^* = R/R_0$, $\dot{R}^* = \dot{R} \sqrt{\rho_l/p_\infty}$, $\ddot{R}^* = \ddot{R} \rho_l R_0/p_\infty$, $R_e = \rho_l u_c R_0/\mu_l$ the Reynolds number with $u_c = \sqrt{p_\infty/\rho_l}$, $W_e = \rho_l u_c^2 R_0/\sigma$ the Weber number.

Note: we consider $p_\infty$ independant of the time.
*/

my_vec f_rk4(double t, my_vec Y){   // dY/dt = f_rk4(t,Y) 
  my_vec Y2;
  Y2.y0 = Y.y1;
  Y2.y1 = (pv-1.)/(Y.y0) + pg0/(Y.y0)*pow(1./Y.y0,3.*gammaG) - 3./2.*sq(Y.y1)/Y.y0 - 4./Re*Y.y1/sq(Y.y0) - 2./We*1./sq(Y.y0); // R-P eq.
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



