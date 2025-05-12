/**
# Navier-Stokes equation with an embedded geometry based on a STL model



Several unexplained issues:
- No vortices (timestep too big ?)
- Restart works but fills cells with crappy values
- Properties in centered.h modifies uf with incoherent values...
*/

#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "distance.h"
#include "utils.h"
#include "view.h"
#include "lambda2.h"
#define Reynolds 100

#define MAXLEVEL 9

u.n[bottom] = dirichlet(1.);
p[bottom]   = neumann(0.);
pf[bottom]  = neumann(0.);

/**
Outflow uses standard Neumann/Dirichlet conditions.  */
        
u.n[top]  = neumann(0.);
p[top]    = dirichlet(0.);
pf[top]   = dirichlet(0.);

u.t[back] = 0.; // ground
u.r[back] = 0.;


/**
Boundary conditions for the solid and fraction tracers. */

u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);
#if dimension == 3
u.r[embed] = dirichlet (0.);
#endif

double eps = 1.e-4; // precision on the geometry

face vector muc[];

double nu = 1./Reynolds; // 60 is the size of the blade of the Wind Turbine

int main()
{
  mu = muc;
  init_grid (32); 
  L0 = 800.;
  origin (-L0/2.,-L0/2.,0.); // centered in x,  centered in y, z is altitude
  DT = 1.; // (uref= 1/ lref = 100) * 1.e-2
  run();
}

event init(i=0){
  // if (!restore (file = "dump")) {

    scalar d[];

// /**
// This part is adapted from the distance.c example.
// */
    system ("test -f windTurbineGeom.stl || "
    "wget \"https://drive.google.com/file/d/12rucJ8LEoFZmXBpch7fB-1buX13ZBg3V\" -O windTurbineGeom.stl"); 
    coord * windTurbine = input_stl (fopen ("windTurbineGeom.stl", "r"));
    distance (d, windTurbine);
    fprintf(stderr, "%ld\n", grid->tn);
    for (int j = 0; j < 5; ++j)
    {
      adapt_wavelet ({d}, (double[]){eps*L0}, MAXLEVEL);
      fprintf(stderr, "%d %ld\n", j, grid->tn);
    }
    vertex scalar phi[];
    foreach_vertex()
      phi[] = -(d[] + d[-1] + d[0,-1] + d[-1,-1] +
         d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
    boundary ({phi});
    fractions (phi, cs, fs);
    fractions_cleanup(cs,fs);
/**
Simple flow in the y direction
*/
    foreach(){
      u.x[] = 0.;
      u.z[] = 0.;
      u.y[] = cs[] ? 1. : 0.;
    }
    boundary((scalar *){u});
    // dump();
    // exit(1);
  // }
}


event properties (i++) {
  foreach_face()
    muc.x[] = fs.x[]*nu;
  boundary ((scalar*){muc});
}

event stability(i++){
// the event properties in centered.h modifies uf with incoherent values...
  foreach_face(){
    if(fs.x[] == 0. && uf.x[] !=0.)uf.x[] =0.;
  }
}

event snapshot (i += 100)
{
  char name[80];
  sprintf (name, "dump-%d", i);
  dump (file = name);
}

event movie (t += 1.; t<=200){
  scalar l2[];

  view (fov = 10,quat = {0.575, -0.191, -0.261, 0.752},ty = -0.2);
  draw_vof ("cs","fs");
  lambda2 (u, l2);
  cells(n= {1,0,0});
  squares("u.y",n= {1,0,0}, min = -0.4, max = 1.25);
  isosurface ("l2", -0.002, color = "u.y", min = -0.4, max = 1.25);
  save ("l2.mp4");
}

event adapt(i++){
  double uemax = 2.e-2;
  adapt_wavelet ({cs,u},
     (double[]){1.e-3 ,uemax,uemax,uemax}, MAXLEVEL, 5);
  fractions_cleanup(cs,fs);
  fprintf(stderr, "%g %ld\n", t, grid->tn);
  // dump();
  // exit(1);
}
