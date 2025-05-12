/**
# Atomisation of a pulsed liquid jet

A dense cylindrical liquid jet is injected into a stagnant lighter
phase (density ratio 1/27.84). The inflow velocity is modulated
sinusoidally to promote the growth of primary shear
instabilities. Surface tension is included and ultimately controls the
characteristic scale of the smallest droplets.

We solve the two-phase Navier--Stokes equations with surface
tension. We need the *tag()* function to count the number of
droplets. We generate animations online using Basilisk View. */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "view.h"
#include "../thinning.h"

/**
We define the radius of the jet, the initial jet length, the Reynolds
number and the surface tension coefficient. */

#define radius 1./12.
#define length 0.025
#define Re 5800
#define SIGMA 3e-5

/**
The default maximum level of refinement is 10 and the error threshold
on velocity is 0.1. */

int maxlevel = 10;
double uemax = 0.05;

scalar skeleton[];
/**
To impose boundary conditions on a disk we use an auxilliary volume
fraction field *f0* which is one inside the cylinder and zero
outside. We then set an oscillating inflow velocity on the
left-hand-side and free outflow on the right-hand-side. */

scalar f0[];
u.n[left]  = dirichlet(f0[]*(1. + 0.05*sin (10.*2.*pi*t)));
u.t[left]  = dirichlet(0);
#if dimension > 2
u.r[left]  = dirichlet(0);
#endif
p[left]    = neumann(0);
f[left]    = f0[];

u.n[right] = neumann(0);
p[right]   = dirichlet(0);

/**
The program can take two optional command-line arguments: the maximum
level and the error threshold on velocity. */

int main (int argc, char * argv[])
{
  skeleton.restriction  =  myrestrict;
  skeleton.prolongation =  myprolongation;

  if (argc > 1)
    maxlevel = atoi (argv[1]);
  if (argc > 2)
    uemax = atof (argv[2]);

  /**
  The initial domain is discretised with $64^3$ grid points. We set
  the origin and domain size. */
  
  init_grid (64);
  origin (0, -1.5, -1.5);
  size (3.);

  /**
  We set the density and viscosity of each phase as well as the
  surface tension coefficient and start the simulation. */
  
  rho1 = 1., rho2 = 1./27.84;
  mu1 = 2.*radius/Re*rho1, mu2 = 2.*radius/Re*rho2;  
  f.sigma = SIGMA;

  run();
}

/**
## Initial conditions */

event init (t = 0) {
  if (!restore (file = "restart")) {

    /**
    We use a static refinement down to *maxlevel* in a cylinder 1.2
    times longer than the initial jet and twice the radius. */
    
    refine (x < 1.2*length && sq(y) + sq(z) < 2.*sq(radius) && level < maxlevel);
    
    /**
    We initialise the auxilliary volume fraction field for a cylinder of
    constant radius. */
    
    fraction (f0, sq(radius) - sq(y) - sq(z));
    f0.refine = f0.prolongation = fraction_refine;
    restriction ({f0}); // for boundary conditions on levels
    
    /**
    We then use this to define the initial jet and its velocity. */

    foreach() {
      f[] = f0[]*(x < length);
      u.x[] = f[];
    }
    boundary ({f,u.x});
  }
}

/**
## Outputs

We log some statistics on the solver. */

event logfile (i++) {
  if (i == 0)
    fprintf (ferr,
         "t dt mgp.i mgpf.i mgu.i grid->tn perf.t perf.speed\n");
  fprintf (ferr, "%g %g %d %d %d %ld %g %g\n", 
       t, dt, mgp.i, mgpf.i, mgu.i,
       grid->tn, perf.t, perf.speed);
}

event snapshot (t = 0.6; t += 0.6; t <= 1.8) {
  char name[80];
  sprintf (name, "snapshot-%g", t);
  scalar pid[];
  foreach()
    pid[] = fmod(pid()*(npe() + 37), npe());
  boundary ({pid});
  foreach(){
    if(f[]> 1.e-3){
      skeleton[] = 1;
    }
    else{
     skeleton[] = 0;
   }
  }
  boundary({skeleton});
  restriction({skeleton});

  thinning2D(skeleton);
  /**
  Simple method to select only cells near the interface. Here the idea is to
  remove the skeleton in the bulk and only select the film. Once we have the
  skeleton we can remove cells which are further than `n` cells away from the
  interface using a `n` pass algorithm and tagging functions.
  */
  int n = 3;
  scalar tag[], tag2[];
  /**
  First we tag interfacial cells.
  */
  foreach(){
    if(f[]> 1.e-3 && f[] < 0.999){
      tag[] = 1;
    }
    else 
      tag[] = 0;
  }
  boundary({tag});

  /**
  We propagate the tag value.
  */
  for (int i = 0; i < n; i++){
    foreach(){
      tag2[] = tag[];
    }
    boundary({tag2});
    
    foreach(){
      if(f[] > 1.e-3 && tag[] == 0){
        double maxval = 0;
        foreach_neighbor(1)
          maxval = max(maxval,tag2[]);
        tag[] = max(maxval,tag[]); 
      }
    }
    boundary({tag});
  }

  /**
  After the propagation we remove the skeleton in the bulk.
  */
  foreach(){
    if(tag[] == 0) skeleton[] =nodata;
  }
  boundary({skeleton});

  dump (name);
}




/**
We generate an animation using Basilisk View. */

event movie (t += 1e-2; t<2.)
{
  scalar omega[];
  vorticity (u, omega);
  view (fov = 7.1, tx = -0.66/2.7*t);
  clear();
  draw_vof ("f");
  squares ("omega", linear = true, spread = 10);
  box ();
  save ("movie.mp4");
  
  scalar skel[];
  foreach(){
    if(skeleton[] == 1 && point.level > maxlevel-2)skel[] = 1;
    else skel[] = nodata;
  }
  boundary({skel});

  draw_vof("f");
  squares("skel");
  save("skeleton.mp4");
}

  
/**
## Mesh adaptation

We adapt the mesh according to the error on the volume fraction field
and the velocity. */

event adapt (i++) {
  // let's build a first skeleton with no a priori estimates
  foreach(){
    if(f[]> 1.e-3){
      skeleton[] = 1;
    }
    else{
     skeleton[] = 0;
   }
  }
  boundary({skeleton});
  restriction({skeleton});

  thinning2D(skeleton);

  // now remove skeleton too far from the interface
  // foreach(){
  //   if(h.x[] == nodata && h.y[] == nodata) skeleton[] = 0.;
  // }
  // boundary({skeleton});
  // restriction({skeleton});

  // scalar m[];
  // foreach()
  //   m[] = f[] > 1e-3;
  // int n = tag (m);

  // *
  // Once each cell is tagged with a unique droplet index, we can easily
  // compute the volume *v* and position *b* of each droplet. Note that
  // we use *foreach_leaf()* rather than *foreach()* to avoid doing a
  // parallel traversal when using OpenMP. This is because we don't have
  // reduction operations for the *v* and *b* arrays (yet). 

  // double v[n];
  // for (int j = 0; j < n; j++)
  //   v[j] 0.;
  // foreach_leaf()
  //   if (m[] > 0) {
  //     int j = m[] - 1;
  //     v[j] += dv()*f[];
  // }
  // foreach(){

  // }
  int n = 3;
  scalar tag[], tag2[];
  foreach(){
    if(f[]> 1.e-3 && f[] < 0.999){
      tag[] = 1;
    }
    else 
      tag[] = 0;
  }
  boundary({tag,tag2});
  for (int i = 0; i < n; i++){
    foreach(){
      tag2[] = tag[];
    }
    boundary({tag2});
    foreach(){
      if(f[] > 1.e-3 && tag[] == 0){
        double maxval = 0;
        foreach_neighbor(1)
          maxval = max(maxval,tag2[]);
        tag[] = max(maxval,tag[]); 
      }
    }
    boundary({tag});
  }

  foreach(){
    if(tag[] == 0) skeleton[] =0;
  }
  boundary({skeleton});

  adapt_wavelet ({f,u,skeleton}, (double[])
    {0.005,uemax,uemax,uemax,0.1},
   maxlevel);
}

/**

## Results


![Atomisation of a pulsed liquid jet. 4096^3^ equivalent resolution.](atomisation_skeleton/movie.mp4)(loop)

![Associated skeleton](atomisation_skeleton/skeleton.mp4)(loop)

*/
