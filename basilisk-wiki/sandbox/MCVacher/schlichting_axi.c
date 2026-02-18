/**
# Schlichting Jet (axisymmetric case)

We try to recover the famous problem of the laminar 2D jet in an immersed domain, also called the Schlichting jet. We look at the axisymmetric case. This analytical solution can be found in H. Schlichting's book "Boundary layer theory" (McGraw-Hill) in section XI.a.2. It was derived by H. Schlichting himself.

It is a monophasic problem depending only on the Reynolds number $Re$. Here we take $Re = 30$.

**Reminder** : *"For problems with a symmetry of revolution around the z-axis of a cylindrical coordinate system : The longitudinal coordinate (z-axis) is x and the radial coordinate (r-axis) is y"*.

Here the videos are rotated for convenience (using the most infamous manner I found ^^).

*/

/**
## Simulation
*/

#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "axistream.h"
#include "view.h"
#define BVIEW 1

scalar s[];
scalar * tracers = {s};

double U0;
double R_d;

double Re;

FILE * fpmax; 

face vector muv[];

int main() {

  Re=30;

  R_d=0.025; //Rayon du jet initial
  L0=1; //Taille de la boite  
  U0=1;

  TOLERANCE = 1e-3 [*];

  u.n[left] = dirichlet (U0*(y < R_d));
  u.t[left] = dirichlet(0.);
  s[left] = dirichlet (U0*(y < R_d));

  u.n[right] = neumann(0);
  p[right] = dirichlet(0.);

  u.n[top] = neumann(0);
  p[top] = dirichlet(0);

  N=128;
  origin (L0, 0);
  init_grid(N);

  mu=muv;

  fpmax =  fopen("log.dat", "w"); 

  run();
}


event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*U0*R_d/Re;
}


event logfile (i++) { 
  fprintf (stderr, "%d %g \n", i, t); 
  fprintf (fpmax, "%d %g\n", i, t);
}

//Ask Antoon what is going on, but it works !
#define vertex_value(b) ((b[] + b[-1] + b[0,-1] + b[-1,-1])/4.)

event movie_streamlines (t += 0.05; t <= 30) {
  scalar omega[];
  vertex scalar stream[];
  scalar psi[];

  psi[left] = dirichlet (0.);
  psi[right]    = neumann (0.);
  psi[top]   = neumann (0.);
  
  vorticity(u, omega);  
  axistream(u, psi);
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = vertex_value(psi);
  boundary ({omega, phi});

  clear();
  view(fov=24, 
       width = L0/2,
       height = L0/2,
       tx     = -1.45*L0,      
       ty     = -0.4*L0);
  
  isoline("phi", n = 30);

  // Draw velocity vectors optionally
  // vect("u", scale=0.1);

  box();
  save("streamlines.mp4");
}

event ppm_output (t = 0; t += 0.05; t <= 30) {

    char name1[80];
    sprintf (name1, "uX.mp4");
    output_ppm (u.x, file = name1, n = 512, min = -U0, max = +U0, linear = true);

    char name2[80];
    sprintf (name2, "s.mp4");
    output_ppm (s, file = name2, n = 512, min = 0., max = U0, linear = true);
}

event profile (t = end) {
  char name[80];
  sprintf (name, "res_end.txt");
  FILE * fpres = fopen(name, "w");
  foreach()
    fprintf (fpres, "%g %g %g %g \n", x, y, u.x[], u.y[]);
  fclose(fpres);
  
  printf ("-----END-----\n");
}


/**
<div style="display:flex; flex-direction:column; gap:60px; align-items:center;">

  <video src="schlichting_axi/s.mp4"
         controls
         loop
         muted
         playsinline
         style="height:400px; transform: rotate(90deg);">
  </video>

  <video src="schlichting_axi/uX.mp4"
         controls
         loop
         muted
         playsinline
         style="height:400px; transform: rotate(90deg);">
  </video>
  
  <video src="schlichting_axi/streamlines.mp4"
         controls
         loop
         muted
         playsinline
         style="height:400px; transform: rotate(90deg);">
  </video>

</div>

*/

/**
## Comparison with theory

*soon*

*/

