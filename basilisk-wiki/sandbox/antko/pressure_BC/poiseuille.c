/**
# Changing a Dirichlet condition for the velocity into a Neumann condition for the pressure

In this example we look at changing a Dirichlet condition for the velocity into a Neumann condition for the pressure.
The setup used for this investigation is a classic Poiseuille entrance flow. A fluid is injected in a 2D pipe with a uniform velocity profile. Upon growth and merging of the vortical boundary layers, the velocity profile ultimately becomes parabolic. 
Here, we look specifically at the entrance boundary condition.

## How to change a velocity boundary condition into a pressure one?

Looking at the projection step used in basilisk, we see that
$$
u^{n+1} = u^\star - \frac{\Delta t}{\rho} \nabla p
$$
where $u^\star$ is formally:
$$
u^\star = u^n + \Delta t \left( \frac{\mu}{\rho} \Delta u^n - \frac{1}{\rho} \left. (u \cdot \nabla) u \right|^n + f^n \right)
$$
(actually the computation of $u^\star$ is more involved -- implicit computation of the viscous term, second order scheme).
*/

#include "navier-stokes/centered.h"

int main() {
  L0 = 8.;
  origin (-0.5, -L0/2.);
  N = 512;
  run(); 
}

#ifndef NEUMANN
u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
#else
p[left]  = neumann(-(a.x[]-1./dt)*fm.x[]/alpha.x[]);
pf[left]  = neumann(-(a.x[]-1./dt)*fm.x[]/alpha.x[]);
#endif

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);

event init (t = 0) {

  mask (y >  0.5 ? top :
	y < -0.5 ? bottom :
	none);
  
  const face vector muc[] = {0.25,0.25};
  mu = muc;
  foreach() {
    u.x[] = 1.;
  }
}


#if 1
event movies (i += 1; t <= 0.1) {
  #ifdef NEUMANN
  static FILE * fpgif = popen ("ppm2gif > pressure_neum.gif", "w");
  #else
  static FILE * fpgif = popen ("ppm2gif > pressure.gif", "w");
  #endif
  output_ppm (p, fpgif, n = 128, box = {{-0.5,-0.5},{0.5,0.5}},
	      min = 0, max = 80, linear = true);
}
#endif

event output_velocity_profiles (t += 0.01; t<=0.1) {
   #ifdef NEUMANN
   static FILE * fp = fopen ("poiseuille_profiles_neum.dat", "w");
   #else
   static FILE * fp = fopen ("poiseuille_profiles.dat", "w");
   #endif
   for (double y = -0.5; y <= 0.5; y += 0.0156249) {
     fprintf (fp, "%.4f %.4f\n", y, interpolate (u.x, 0., y));
   }
   fprintf (fp, "\n\n");
}

event output_pressure_profiles (t = end) {
   #ifdef NEUMANN
   static FILE * fppres = fopen ("poiseuille_pressure_profile_neum.dat", "w");
   #else
   static FILE * fppres = fopen ("poiseuille_pressure_profile.dat", "w");
   #endif
   for (double x = -0.5; x <= 7.5; x += 1e-2) {
     fprintf (fppres, "%.4f %.4f\n", x, interpolate (p, x, 0.));
   }
   fprintf (fppres, "\n\n");
}

#if 0
event gfsview (i += 1) {
  static FILE * fp = popen ("gfsview2D poiseuille_velocity.gfv", "w");
  output_gfs (fp, t = t);
}
#endif

/**
![Velocity profiles](/sandbox/antko/pressure_BC/pics/profiles.png)
![Pressure profiles](/sandbox/antko/pressure_BC/pics/pressure_profiles.png)
![Pressure field at entrance with a Dirichlet condition for the velocity](/sandbox/antko/pressure_BC/pics/pressure.gif)
![Pressure field at entrance with a Neumann condition for the pressure](/sandbox/antko/pressure_BC/pics/pressure_neum.gif)
*/