/**
This code simulate Rayleigh-Bénard instability with Strongly Temperature-
Dependent Viscosity. Convection patterns changes for hexagonal or squares rolls.
*/
#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "navier-stokes/perfs.h"
//#include "view.h"
//#include "display.h"


scalar T[];
scalar * tracers = {T};
face vector av[], muc[];

double Ttop = 20;
double Tbot = 35;
#define T0 ((Tbot + Ttop)/2.)
#define deltaT (Tbot - Ttop)

double Ra, Pr;
double EndTime= 50.;
int ii;

int main() {
  size (npe()); 
  Y0 = -0.5;
  dimensions(nx = npe(), ny = 1);
  TOLERANCE = 1e-3;
  init_grid (256);
  mu = muc; //viscosity
  a = av; // acceleration
  Ra = 3000; Pr = 1000.;
  run();
}

/**

## Boundary conditions 

no-slip walls and fixed temperature at the top and the bottom. In addition Neumman condition for the temperature on the side of the domain.

*/
T[top] = dirichlet(-0.5);
T[left] = neumann(0.);
T[right] = neumann(0.);
T[bottom] = dirichlet(0.5);

u.n[top] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.n[right] = dirichlet(0.);
u.n[left] = dirichlet(0.);
u.t[top] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.t[left] = dirichlet(0.);

/**
## Initialization
*/

event init (t=0) {
  foreach(){
    T[] = - y;
    foreach_dimension()
      u.x[] = 0.01*noise();
  }
  DT = 0.001;
  dtnext(DT);
  TOLERANCE=10E-6;
  //boundary ({T,u});
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
  foreach_face(y)
    av.y[] += Pr*((T[] + T[0,-1])/2.);
    
  if ((i==10||ii==10)){
    DT=0.1;
    dtnext(DT);
    TOLERANCE=10e-4;
  }
}

/**
##  Outputs
*/

event movies(t+=1.0; t<=EndTime){
  output_ppm (T, file="temperature.mp4", n = 1280, min = -0.48, max = 0.48, linear = true, box = {{0,-0.5},{L0, 0.5}});
}

event shadowgraph(t = EndTime){
  scalar shadow[];
  foreach()
    shadow[] = (T[1,0] -2*T[] + T[-1,0])/sq(Delta) + (T[0,1] - 2*T[] + T[0,-1])/sq(Delta);
  output_ppm (shadow, file="shadow.png", n = 1280, spread = 1, linear = true, map = gray, box = {{0,-0.5},{L0, 0.5}});
}

/**
## Results
![Temperature field.](rbvisco/temperature.mp4)
![Numerical shadowgraph.](rbvisco/shadow.png)


## References
D. S. Oliver & J. R. Booker (1983) Planform of convection with
strongly temperature-dependent viscosity, Geophysical & Astrophysical Fluid Dynamics,
27:1-2, 73-85, DOI: 10.1080/03091928308210121
*/