#include "grid/octree.h"
#include "navier-stokes/centered.h"
#define TWO_WAY 1
#include "stokes-particles.h"
#include "view.h"
#include "scatter2.h"  

u.n[front] = dirichlet (0.7);
u.r[front]  = neumann(0);
u.r[front]  = neumann(0);
p[front]   = neumann(0);
pf[front]  = neumann(0);


u.n[back] = neumann(0);
u.t[back] = neumann(0);
u.r[back] = neumann(0);
p[back]   = dirichlet(0);

Particles drop;

#define MUZ 1.8e-5
#define U0 27.6 //m/s
#define rout 1e-4 //m.

int main() 
{
  const face vector muc[] = {MUZ, MUZ, MUZ}; //mu Water
  mu = muc;
  G.y = -9.81;
  L0 = 1.1;
  X0 = Y0 = Z0 = -L0/2;
  mu = muc;
  N = 64;
  run();
  DT = 0.00006;
}

// Initially, the flow is steady and seeded with particles near the top of the domain.

event init (t = 0) 
{
  const scalar rhof[] = 1.2;
  rho = rhof;
  drop = new_inertial_particles (0);
  //add one particle in the system
  particle p = { .x = 0, .y = 0, .z = 0.15, 
        .u.x = 0.5*U0*noise()*tan(pi/6), .u.y = 0.5*U0*noise()*tan(pi/6), .u.z = -U0,
        .u2.x = 1000, .u2.y = 20e-6, .tag = 1*t};
  add_particle (p, drop);
}

//remove particle if goes out of box
event remover (i++) 
{
  remove_particles (drop, z > L0/2 || y > L0/2 || x > L0/2 || z < -L0/2 || y < -L0/2 || x < -L0/2);
}


event adapt (i++) 
  adapt_wavelet ((scalar*){u}, (double[]){0.01, 0.01, 0.01}, 9, 5);


//create movie
event mov (t += .02) {
  view (theta = -1.5);
  scatter (drop, pc = {1, 1, 1});
  box();
  translate (z = -L0/2)
    cells();
  save ("locs.mp4");
}

event end (t = 10);

