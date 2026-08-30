/** # Simulation of the jet-edge interaction with embedded boundary

This code investigates the dynamics of a planar jet impinging on a shaped plate (wedge/edge) using 
the embedded boundary method (`embed.h`). 

*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

double Reynolds = 150; 
int maxlevel = 7;
double b = 1.; // Jet diameter
double U0 = 1.; // Reference velocity (top-hat profile)
double W = 5.; // Jet-to-plate distance 
double freq = 0.067; // Forcing frequency
double OMEG ; // Forcing angular frequency
double eps = 1e-3; // Perturbation amplitude
double delta_erf; // Shear layer thickness
double h, R;
double Tramp; // Ramp-up duration

scalar omega[];
scalar base[];
scalar un[];
face vector muv[];
double du;
double Tend = 500;

face vector av[];

FILE * fp;

int main(int argc, char * argv[])
{
  /** Interactive display controls configuration. */
  display_control(Reynolds,0,1000);
  display_control(maxlevel,0,15);

  
  L0 = 50; 
  origin (0., -L0/2.);
  N = 512; 
  mu = muv;
  a = av;
  
  delta_erf = b/10.;
  OMEG = 2*M_PI*freq;
  Tramp = 50; // Time ramp duration to avoid numerical shocks

  init_grid (N);

  fp = fopen("data.dat","w");
    
  dt = 0.1;
  DT = 0.1;
  run();
}

/** Viscosity and vorticity field calculation at each time step. */
event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*(b*U0)/Reynolds;
  vorticity (u, omega);
}

/** Inlet boundary conditions (left):
    Smoothed velocity profile (erf error functions) modulated by a temporal 
    harmonic perturbation and a progressive activation ramp. */
u.n[left] = dirichlet(
    (min(t/Tramp,1.))
    *
    (0.5*(erf((y+b/2.)/delta_erf)
          -erf((y-b/2.)/delta_erf))
     *(1.+eps*cos(OMEG*t)))
);

u.t[left] = dirichlet((min(t/Tramp,1.)) * eps * sin(OMEG*t));

/** Outlet and lateral confinement boundary conditions. */
u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);

/** No-slip boundary condition on the embedded solid interface. */
u.n[embed] = dirichlet(0.); 
u.t[embed] = dirichlet(0.);

/** Solid edge geometry initialization using the `solid()` function. */
event init (t = 0)
{
  h = 0.1;
  R = h/2;

  /** Analytical geometry definition of the solid obstacle. */
  solid (cs, fs,
         min(
             max(max(W - x, x - L0),
                 max(y - h/2., -y - h/2.)),

             max(sq(x - W) + sq(y) - sq(R),
                 x - W)
         ));
}

/** Time-series logging and extraction of velocity/pressure probe values. */
event logfile (i++, t<=Tend) {
  vorticity(u,omega);

  /** Standard error output (console). */
  fprintf (stderr, "%d %g %g %g %g %g %g %g %g\n", 
           i, t, dt, 
           interpolate(u.y, 0.1*W, 0.), 
           interpolate(u.y, 0.9*W, 0.),
           interpolate(u.x, 0.05*W, 0.),
           interpolate(p, 0.9*W, W/2),
           interpolate(p, 0.9*W, -W/2),
           interpolate(omega, 0.9*W, 0.));  

  /** Writing to external output file `data.dat`. */
  fprintf (fp, "%d %g %g %g %g %g %g %g %g\n", 
           i, t, dt, 
           interpolate(u.y, 0.1*W, 0.), 
           interpolate(u.y, 0.9*W, 0.),
           interpolate(u.x, 0.05*W, 0.),
           interpolate(p, 0.9*W, W/2),
           interpolate(p, 0.9*W, -W/2),
           interpolate(omega, 0.9*W, 0.));  
}

/** Movie generation for vorticity field dynamics. */
event movie (t += 2; t <= Tend)
{
  scalar m[];

  foreach()
    m[] = cs[];

  vorticity (u, omega);

  output_ppm (omega,
              file = "vorticity.mp4",
              n = N,
              min = -0.5, max = 0.5,
              linear = true,
              mask = m);
}

/** High-resolution final frame rendering with solid boundary line drawing. */
event final_images (t = Tend)
{
  vorticity (u, omega);

  clear();

  view (
      width = 1600,
      height = 900,
      fov = 18,
      tx = -0.5,
      ty = 0.
  );

  squares ("omega",
           min = -1,
           max = 1,
           linear = true);

  draw_vof ("cs", lw = 2);

  box();

  save ("omega_final.png");
}

/** Initial state snapshot of the solid geometry. */
event init_images (t = 0.)
{
  output_ppm (cs,
              file = "solid.png",
              n = 1024,
              box = {{W-0.2,-0.2},{W+0.2,0.2}},
              min = 0,
              max = 1,
              linear = false);
}


/** Spatial snapshot extraction for Dynamic Mode Decomposition (DMD) analysis. */
/*
event dmd_snapshots (t = 250; t += 0.5; t <= Tend) 
{
  vorticity (u, omega);

  char filename[80];
  sprintf (filename, "snapshot_t%.2f.dat", t);
  FILE * fp_dmd = fopen (filename, "w");

  int N_dmd = 256; 

  // Spatial extraction domain defined around the edge
  coord box[2] = {
    {0.0, -5.0},   // Lower left corner [x_min, y_min]
    {20,  5.0}    // Upper right corner [x_max, y_max]
  };

  output_field ({omega}, fp_dmd, N_dmd, box = box, linear = true);

  fclose (fp_dmd);
}*/





/**
# Visualisations


![Animation of the vorticity field](edge_tone/vorticity.mp4)(loop)


*/