#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "embed.h"
#include "view.h"

face vector av[];
face vector muv[];

int maxlevel = 8; //maximum refinement level
double be = 1e-2, ue = 3e-2; //max. discretisation error

double nu = 1e-3; //kinematic viscosity
double U0 = 1.; //intial flow velocity
double end = 25.; //simulation end time

int main() {
  L0 = 8; //domain size
  X0 = Y0 = -L0/2.;
  N = 256;
  mu = muv;
  periodic(left); //periodic left-right boundary condition
  run();
}

//no-slip boundary conditions for lower boundary
u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);

event init (t = 0)
{
  foreach()
    u.x[] = U0; //initialise fluid velocity
}

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*nu; //set fluid viscosity
}

event adapt (i++)
  adapt_wavelet ({u}, (double[]){be,ue}, maxlevel, 4);

// output vorticity movie with cells
event movies (i += 4; t <= end)
{
  scalar omega[];
  vorticity (u, omega);
  view (camera = "front", width = 1024, height = 1024);
  squares("omega", linear = true);  
  box();         
  cells();           
  save ("vort.mp4");
}

// output full velocity field (u,v) for further analysis
event turbulence(t = 0, t <= end, t+=5) {
  output_matrix (u.x, fp = fopen("u", "a"), linear=true);
  output_matrix (u.y, fp = fopen("v", "a"), linear=true);
}