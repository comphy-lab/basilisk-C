#include "grid/quadtree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "view.h"

face vector av[];
face vector muc[];

int maxlevel = 8; //maximum refinement level
double be = 1e-5, ue = 1e-2; //max. discretisation errors

double nu = 1e-5; //kinematic viscosity
double rhograss = 1.; //grass density per unit distance
double hgrass = 1.; //grass height 
double U0 = 1.; //fluid velocity
double end = 10.; //simulation end time

// this function creates the "grass"
void grass (scalar cs, face vector fs)
{
  vertex scalar phi[];
  vertex scalar test[];
  foreach_vertex() {
    phi[] =  y - 0.3;
    for (int xp = 0; xp < L0*rhograss; xp++) {
      double d = xp/rhograss;
      phi[] = intersection(phi[], union(fabs(x-d) > 0.05, y > (hgrass + 0.3)));
    }
  }
  fractions (phi, cs, fs);
}

int main() {
  L0 = 8; //domain size
  X0 = Y0 = 0;
  N = 512;
  mu = muc;
  a = av;
  periodic(left); //periodic left-right boundary condition
  run();
}

//set fluid viscosity
event properties (i++)
{
  foreach_face()
    muc.y[] = fm.y[]*nu;
  boundary((scalar *){muc});
}

//embedded (grass) and ground no-slip boundary conditions
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);

//initialise fluid velocity and embedded grass
event init (t = 0)
{
  grass (cs, fs);
  foreach(){
    u.x[] = cs[] ? U0 : 0.;
  }
}

event adapt (i++)
  adapt_wavelet ((scalar *){cs, u}, (double[]){be, ue, ue}, maxlevel, 4);

//output vorticity movie and grid cells with a colorbar
event movies(i += 4, t <= end)
{
  scalar omega[];
  vorticity (u, omega);
  view (tx = -0.5, ty = -0.5, camera = "front", width = 512, height = 512, fov = 20); 
  squares("omega", linear = true, cbar = true, mid = true, horizontal = true, 
          border = true, levels = 10, pos = {-0.5, 0.7}, format = "%9.3f", size = 10);                           
  box();                    
  save ("vort.mp4");
}

//output the full velocity field at t=1s timesteps
event outputfile(t = 0, t <= end, t+=1) {
  output_matrix (u.x, fp = fopen("u", "a"), linear=true);
  output_matrix (u.y, fp = fopen("v", "a"), linear=true);
}