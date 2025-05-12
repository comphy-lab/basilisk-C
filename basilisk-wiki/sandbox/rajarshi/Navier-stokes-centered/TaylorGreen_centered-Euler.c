/**
# Taylor-Green vortices --> Centered NS - WENO5 + Poisson4 + Euler.
*/

#include "grid/multigrid.h"
#include "centered_weno5_Poisson4-Euler.h"

scalar omega[],ke[];

int main() {
 
  Prolongation_Weight_Initialization();
  origin (-0.5,-0.5);
  foreach_dimension()
    periodic (right);
  
  for (N = 64; N <= 64; N *= 2)
    run();
}

event init (i = 0) {
 
  foreach() {
    u.x[] =   (sin(2.*pi*(x+Delta/2.))-sin(2.*pi*(x-Delta/2.)))*(cos(2.*pi*(y+Delta/2.))-cos(2.*pi*(y-Delta/2.)))/(sq(2.*pi*Delta));
    u.y[] = - (cos(2.*pi*(x+Delta/2.))-cos(2.*pi*(x-Delta/2.)))*(sin(2.*pi*(y+Delta/2.))-sin(2.*pi*(y-Delta/2.)))/(sq(2.*pi*Delta));
    p[]   = - (sin(4.*pi*(x+Delta/2.))-sin(4.*pi*(x-Delta/2.)) + sin(4.*pi*(y+Delta/2.))-sin(4.*pi*(y-Delta/2.)))/(16.*pi);
  }

  boundary ({p});
  boundary ((scalar *){u});
  foreach()
    foreach_dimension()
      g.x[] = - (p[-2] -8.*p[-1] + 8.*p[1] - p[2])/(48.*Delta);
  boundary ((scalar *){g});
}

event animation (i++,t<=1){

   foreach(){
       omega[] =  (u.y[1,0]-u.y[-1,0])/(2.*Delta) - (u.x[0,1]-u.x[0,-1])/(2.*Delta);
       ke[]    =  (sq(u.x[]) + sq(u.y[]) )/2.;
    }
  
    output_ppm(omega,min=-12.5,max=12.5,file="Omega.mp4");
    output_ppm(p,min=-0.5,max=0.5,file="Pressure.mp4");
    output_ppm(ke,min=0,max=0.5,file="Kinetic-Energy.mp4");
}

/**
We produce animations of the vorticity, pressure and kinetic-energy fields 

![Animation of the vorticity field.](TaylorGreen_centered-Euler/Omega.mp4)

![Animation of the pressure field.](TaylorGreen_centered-Euler/Pressure.mp4)

![Animation of the kinetic energy field.](TaylorGreen_centered-Euler/Kinetic-Energy.mp4) 

*/

event error (t = 1) {
 
  scalar e[];
  foreach() {
    double u0 = - cos(2.*pi*x)*sin(2.*pi*y);
    double v0 =   sin(2.*pi*x)*cos(2.*pi*y);
    e[] = norm(u) - sqrt(sq(u0) + sq(v0));
  }
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
}
