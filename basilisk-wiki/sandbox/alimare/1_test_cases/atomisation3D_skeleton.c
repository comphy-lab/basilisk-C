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

/**
 * Compilation :
CC99='h5pcc -std=c99' qcc -O2 -D_MPI=2 atomisation3D_skeleton.c -L$BASILISK/gl -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -lOSMesa -lGLU -lglutils -lfb_osmesa -lm -lhdf5 -lgsl -lgslcblas -o atomisation3D_skeleton/prog

 */
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "view.h"
#include "../thinning.h"
long unsigned int n_part;
long unsigned int out_skel = 0;
#include "../output_skeleton.h"
#include "../output_hdf.h"
#define QUADRATIC 1
#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))
#include "alex_functions.h"
#include "../LS_reinit.h"
int count_it = 0;

// #pragma autolink -lhdf5
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

int maxlevel = 11;
double uemax = 0.05;

scalar cskeleton[];
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

void bulkRemoval(scalar c){
  int filter = 1;
  scalar tag[], tag2[];
  /**
  First we tag interfacial cells.
  */
  foreach(){
    if(f[]> 1.e-4 && f[] < 0.999){
      tag[] = 1;
    }
    else 
      tag[] = 0;
  }
  boundary({tag});

  for (int i = 0; i < filter; i++){
    foreach(){
      tag2[] = tag[];
    }
    boundary({tag2});
    
    foreach(){
      if(f[] > 1.e-4 && tag[] == 0){
        double maxval = 0;
        foreach_neighbor(1)
          maxval = max(maxval,tag2[]);
        tag[] = max(maxval,tag[]); 
      }
    }
    boundary({tag});
  }

  foreach(){
    if(tag[] == 0) c[] =0;
  }
  boundary({c});
  restriction({c});
}

/**
The program can take two optional command-line arguments: the maximum
level and the error threshold on velocity. */

int main (int argc, char * argv[])
{
  cskeleton.restriction  =  myrestrict;
  cskeleton.prolongation =  myprolongation;

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

// event snapshot (t = 0.6; t += 0.6; t <= 1.2) {
//   char name[80];
//   sprintf (name, "snapshot-%g", t);
//   scalar pid[];
//   foreach()
//     pid[] = fmod(pid()*(npe() + 37), npe());
//   boundary ({pid});
//   foreach(){
//     if(f[]> 1.e-3){
//       cskeleton[] = 1;
//     }
//     else{
//      cskeleton[] = 0;
//    }
//   }
//   boundary({cskeleton});
//   restriction({cskeleton});

//   thinning3D(cskeleton);
//   dump (name);
// }




/**
We generate an animation using Basilisk View. */

// event movie (i+=10; i<1000)
// {
//   scalar omega[];
//   vorticity (u, omega);
//   view (fov = 7.1, tx = -0.66/2.7*t);
//   clear();
//   draw_vof ("f");
//   squares ("omega", linear = true, spread = 10);
//   box ();
//   save ("movie.mp4");
  
//   scalar skel[];
//   foreach(){
//     if(cskeleton[] == 1 && point.level > maxlevel-2)skel[] = 1;
//     else skel[] = nodata;
//   }
//   boundary({skel});

//   draw_vof("f");
//   squares("skel");
//   save("cskeleton.mp4");
// }

  
/**
## Mesh adaptation

We adapt the mesh according to the error on the volume fraction field
and the velocity. */

event build_skeleton (i++,last) {
  // let's build a first cskeleton with no a priori estimates
  foreach(){
    if(f[]> 1.e-5){
      cskeleton[] = 1;
    }
    else{
     cskeleton[] = 0;
   }
  }
  boundary({cskeleton});
  restriction({cskeleton});

  thinning3D(cskeleton);

  bulkRemoval(cskeleton);
}

event output(t+=0.01,last;t<4){
  n_part= 0;  
  int n = 0;
  foreach(){
    if(cskeleton[]){
      n_part++;
    }
  }
  fprintf(stderr, "N_PART %lu \n", n_part);
  char name[80];
  double * myfield = malloc (n_part*sizeof(double));
  coord * loc = malloc (n_part*sizeof(coord));
  foreach() {
    if(cskeleton[]){
      coord cc = {x, y, z};
      foreach_dimension()
        loc[n].x = cc.x;
      myfield[n] = pid();
      n++;
    }
  }

  sprintf(name, "skeleton_%ld", out_skel);

  output_pskeleton(loc,myfield,name);

/**We output the results in HDF5 format.*/

  scalar dist[];

  foreach()
   dist[] = (2*f[]-1)*L0/(1 << grid -> maxdepth)*0.75;

  boundary({dist});
  restriction({dist});
  LS_reinit(dist,dt = 0.1*L0/(1 << grid -> maxdepth), it_max = 30);

  scalar * list = {f,dist};
  vector * vlist = {u};
  
  char buf[100];
  char time[100];
  char it[100];
  
  snprintf(buf, sizeof(buf), "out_%05d.xmf", count_it);
  snprintf(it, sizeof(it), "%05d", count_it);
  snprintf(time, sizeof(time), "%06ld", out_skel);
  out_skel++;
  count_it++;
  
  FILE * fp = fopen(buf, "w");
  output_xmf_h5_foreach(list, vlist, 64, fp, it, time); 
  
  fclose(fp);

  if(i%100==0 && i!=0){
    dump();
    fprintf(stderr, "FILE DUMPED\n");
  }


/** interface output*/
  // sprintf(name, "%s-%.5d-CC.vtk", "out", i);
  // output_paraview_IF(name, t, 0.0, f);
}

event adapt(i++){
  adapt_wavelet ({f,u,cskeleton}, (double[])
    {0.005,uemax,uemax,uemax,0.1},maxlevel);

  // output_vtu ((scalar *) {f},(vector *) {NULL},  name);
}

/**

## Results


![Atomisation of a pulsed liquid jet. 4096^3^ equivalent resolution.](atomisation_skeleton/skeleton.mp4)(loop)

![Associated skeleton](atomisation_skeleton/skeleton.mp4)(loop)

*/
