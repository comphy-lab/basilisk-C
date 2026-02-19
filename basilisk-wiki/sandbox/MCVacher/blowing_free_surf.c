/**
# Laminar axi. jet impinging on a free surface

Inspired from original code from the "We Make Waves" team [here](/sandbox/WMW/blow.c).
  
*/

/**
## Simulation
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "curvature.h"
#include "view.h"
#include "navier-stokes/perfs.h"
#include "tracer.h"

scalar s[];
scalar * tracers = {s};

//Dim. variables
double U0;
double R_d;
double h;

//Undim. variables
double Re;
double rhob;
double mub;
double box_ratio;

FILE * fpmax; 

int main()
{
  Re=700;                  // rho_2*U0*R_d/mu2
  rhob=1./1000.;           // air/water
  mub= 0.000018/0.001;     // air/water
  box_ratio=0.1;          // R_d/L0
    
  L0=1; 
  R_d=box_ratio*L0; 
  
  U0=10;
  
  u.n[top] = neumann(0);
  p[top] = dirichlet(0);

  u.n[right] = dirichlet (-U0*(1-f[])*(y < R_d));
  s[right] = dirichlet (U0*(y < R_d));
  
  origin (-L0/2.,0.);
  
  rho1 = 1.;
  rho2 = rhob*rho1;
  
  mu2 = rho2*U0*R_d/Re;
  mu1 = mu2/mub;
 
  G.x = -1;
  
  f.sigma=0.001;
  
  N = 128;
  init_grid(N);
  
  fpmax =  fopen("log.dat", "w"); 
  
  run();
}

event init (i = 0) {
  foreach()
    f[] = x < 0;
  boundary ({f});
}

event movie (t += 0.01; t <= 10) {
  view (fov = 21.0715, quat = {0,0,-0.707,0.707},
	tx = -0.000575189, ty = 0.0132879, bg = {1,1,1},
	width = 1208, height = 666, samples = 4);
  box (notics = true);
  draw_vof ("f", filled = 0, fc = {1,1,1});
  squares ("u.x", linear = true);
  mirror (n = {0,1}) {
    draw_vof ("f", filled = -1, fc = {1,1,1});
    squares ("u.y", linear = true);
    box (notics = true);
  }
  save ("movie.mp4");
}

event ppm_output (t = 0; t += 0.01; t <= 10) {

    char name1[80];
    sprintf (name1, "uX.mp4");
    output_ppm (u.x, file = name1, n = 512, min = -U0, max = +U0, linear = true);
    
    scalar s_other[];
    foreach()
      s_other[] = (1-f[])*s[];
  
    char name2[80];
    sprintf (name2, "s.mp4");
    output_ppm (s_other, file = name2, n = 512, min = 0, max = U0, linear = true);
  
    char name3[80];
    sprintf (name3, "f.mp4");
    output_ppm (f, file = name3, n = 512, min = 0., max = 1, linear = true);
}

event profile (t = end) { 
  printf ("-----END-----\n");
}

/**
![Animation of free surface, vertical velocity (left), 
  horizontal velocity (right).](blowing_free_surf/movie.mp4)(width=604,height=333)

<div style="display:flex; flex-direction:column; gap:60px; align-items:center;">
  
  <video src="blowing_free_surf/f.mp4"
         controls
         loop
         muted
         playsinline
         style="height:400px; transform: rotate(-90deg) scaleY(-1);">
  </video>

  <video src="blowing_free_surf/s.mp4"
         controls
         loop
         muted
         playsinline
         style="height:400px; transform: rotate(-90deg) scaleY(-1);">
  </video>

  <video src="blowing_free_surf/uX.mp4"
         controls
         loop
         muted
         playsinline
         style="height:400px; transform: rotate(-90deg) scaleY(-1) ;">
  </video>

</div>

*/


