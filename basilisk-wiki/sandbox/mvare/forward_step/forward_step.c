/**
# Flows over a flat plate and over a forward step

A uniform flow of constant velocity over a flat plate or a forward step is 
computed using embed boundary conditions for the solid surfaces and 
adaptive mesh refinement (AMR). The prolongation function of the AMR is 
modified to take into account a no-slip condition on the solid boundary.

*/

#include "embed2.h"
#include "navier-stokes/centered.h"
#include "view.h"

double Reynolds = 160.;
int maxlevel = 9;
face vector muv[];

double eps = 1.e-6;
double mygeom(double x, double y){
  /* 3 geometries are considered:
  case 1: flow over a flat plate not coinciding with the grid
  case 2: flow over a flat plate coinciding with the grid
  case 3: flow over a forward facing step not coinciding with the grid*/
  //return(y+L0/2.-L0/20.-eps);   // case 1
  return(y+L0/2.-L0/16.-eps); // case 2
  /*if(x<L0/20.-eps)                // case 3 
  {
    return(y+L0/2.-L0/20.-eps);
  }
  else
  {
    return(y+L0/2.-2.*L0/20.-eps);
  }*/

}

int main() {
  L0 = 8.;
  origin (-L0/2., -L0/2.);
  N = 512;
  mu = muv;
  run(); 

}

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*0.125/Reynolds;
}

// Boundary conditions

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);

u.n[bottom] = dirichlet(0.); 
u.t[bottom] = dirichlet(0.);

p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

event init (t = 0) // initialization
{
  vertex scalar distn[];
 
  foreach_vertex() {
    distn[] = mygeom(x,y);
    }
  boundary({distn});
  fractions(distn,cs,fs);
  fractions_cleanup(cs,fs);
  
  foreach()
    u.x[] = cs[] ? 1. : 0.;

  foreach_dimension()
    u.x.prolongation = refine_embed_prolong_wall; // change of the prolongation functions
}

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

event adapt (i++) {
  adapt_wavelet ({cs,u}, (double[]){1e-2,3e-2,3e-2}, maxlevel, 4);
}

// Drawing of pictures

void movies_core(int i){
 
  scalar omega[],m[];
  vorticity(u,omega);
  foreach()
    m[] = cs[] - 0.5;
  
  draw_vof ("cs", "fs", filled = -1);
  squares("omega", linear = true, min = -10, max = 10);
  char name[80];
  sprintf (name, "vort_%05d.png", i);
  save(name);  

  cells();
  squares("cs",linear = true, min = 0, max = 1);
  sprintf (name, "cells_%05d.png", i);
  save(name);

}

event movies ( i += 4; t<= 10){
  view(fov = 6, width = 1600, height = 500, ty = 0.35, tx = 0.000001);   
  movies_core(i);
}

