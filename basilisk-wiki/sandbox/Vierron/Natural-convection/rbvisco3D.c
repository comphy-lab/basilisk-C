/**
This code simulate Rayleigh-Bénard instability in 3D with Strongly Temperature-
Dependent Viscosity. Convection patterns changes for hexagonal or squares rolls.
*/

#include "grid/multigrid3D.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "navier-stokes/perfs.h"
#include "view.h"
//#include "display.h"
#include "lambda2.h"
//#include "tracer-particles.h"
//#include "scatter2.h"

scalar T[];
scalar * tracers = {T};
face vector av[], muc[];
//Particles flow;

double Ttop = 20;
double Tbot = 26;
#define T0 ((Tbot + Ttop)/2.)
#define deltaT (Tbot - Ttop)

double Ra, Pr;
double EndTime= 20.;
int ii;

int main() {
  size (npe()/2); 
  dimensions(nx = npe()/2, ny = npe()/2, nz = npe()/4);
  TOLERANCE = 1e-3;
  init_grid (64);
  mu = muc; //viscosity
  a = av; // acceleration
  Ra = 3000; Pr = 1000.;
  run();
}

/**

## Boundary conditions 

no-slip walls and fixed temperature at the top and the bottom. In addition Neumman condition for the temperature on the side of the domain.

*/
T[front] = dirichlet(-0.5);
T[left] = neumann(0.);
T[right] = neumann(0.);
T[back] = dirichlet(0.5);
T[top] = neumann(0.);
T[bottom] = neumann(0.);

u.n[top] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.n[right] = dirichlet(0.);
u.n[left] = dirichlet(0.);
u.n[back] = dirichlet(0.);
u.n[front] = dirichlet(0.);
u.t[top] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.t[left] = dirichlet(0.);
u.t[back] = dirichlet(0.);
u.t[front] = dirichlet(0.);

u.r[top] = dirichlet(0.);
u.r[bottom] = dirichlet(0.);
u.r[right] = dirichlet(0.);
u.r[left] = dirichlet(0.);
u.r[back] = dirichlet(0.);
u.r[front] = dirichlet(0.);

/**

## Initialization

Initialization of the fluid layer and linear temperature profile with no motion.
*/
event init (t=0) {
  foreach(){
    T[] = - 4/npe()*z + 0.5;
    foreach_dimension()
      u.x[] = 0.01*noise();
  }
  DT = 0.001;
  dtnext(DT);
  TOLERANCE=10E-6;
  //boundary ({T,u});
  //if(pid()==0)
    //flow = init_tp_square (20);
}

/**

## Strongly temperature-dependent viscosity (glycerol)

*/

#define mu0 (exp(4.549 - 0.12309*T0 + 9.1129e-4 * pow(T0,2) -4.7562e-4 * pow(T0,3) + 1.3296e-8 * pow(T0,4)))
#define mu(T) ((exp(4.549 - 0.12309*((T*deltaT+T0)) + 9.1129e-4 * pow((T*deltaT+T0),2) -4.7562e-4 * pow((T*deltaT+T0),3) + 1.3296e-8 * pow((T*deltaT+T0),4)))/mu0)

event properties (i++) {
  foreach_face(){
    double ff = face_value(T,0);
    muc.x[] = fm.x[]*Pr/sqrt(Ra)*(mu(ff));
  }
  //boundary ((scalar*){muc});
}
/**

## Thermal diffusion

*/
event tracer_diffusion (i++) {
  face vector D[];
  foreach_face()
    D.x[] = fm.x[]*1./sqrt(Ra);
  //boundary ((scalar*){D});
  diffusion (T, dt, D);
  //boundary ({T});
}


/**
## Boussinesq approximation
*/

event acceleration (i++) {
  ii++;
  foreach_face(z)
    av.z[] += Pr*((T[] + T[0,-1])/2.);
    
  if ((i==10||ii==10)){
    DT=0.1;
    dtnext(DT);
    TOLERANCE=10e-4;
  }
}

/**
##  Outputs
*/

scalar l2[];
event movies(t+=0.1; t<=EndTime){
  clear();
  view (quat = {0.243, 0.493, 0.717, 0.430},
      fov = 30, near = 0.01, far = 1000,
      tx = -0.187, ty = 0.081, tz = -3.134,
      width = 1550, height = 940);
  box ();
  lambda2 (u, l2);
  isosurface (f = "l2", v = -0.04, color = "T");
  save("l2.mp4");
}

event video2(t+=0.1; t<=EndTime){
  clear();
  view (quat = {0.243, 0.493, 0.717, 0.430},
      fov = 30, near = 0.01, far = 1000,
      tx = -0.187, ty = 0.081, tz = -3.134,
      width = 1550, height = 940);
  box ();
  isosurface (f = "T", v = 0, color = "u.z", linear = true, min = -0.2, max = 0.2);
  //scatter(flow, s = 5);
  save("test.mp4");
}

#if 1
event shadow(t = EndTime){
  scalar shadow[];
  foreach()
    shadow[] = 0.;
  foreach(){
    if(z == 31 * Delta) 
      shadow[] += ((T[1,0,0] - 2.*T[] + T[-1,0,0])/sq(Delta) + (T[0,1,0] - 2.*T[] + T[0,-1,0])/sq(Delta)) * Delta;
  }
  stats s = statsf(shadow);
  printf("min = %g, max = %g\n", s.min, s.max);
  output_ppm(shadow, file="shadow.png", n = 1280, spread = 1, linear = true, map = gray, min = s.min, max = s.max);
  dump();
}
#endif

/**
## Results
![vortical structures](rbvisco3D/l2.mp4)
![Isosurface of temperature field](rbvisco3D/test.mp4)
![Numerical shadowgraph](rbvisco3D/shadow.png)


## References
D. S. Oliver & J. R. Booker (1983) Planform of convection with
strongly temperature-dependent viscosity, Geophysical & Astrophysical Fluid Dynamics,
27:1-2, 73-85, DOI: 10.1080/03091928308210121
*/