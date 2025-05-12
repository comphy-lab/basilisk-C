/** 
# Stokes flow through a periodic array of spheres with a fictitious-domain method
This computes the dynamic of a flow through an infinite array of
spheres.  The concentration is set to $c=0.45$ and the expected drag
coefficient is $K=28.1$ (analytical solution by [Zick and Homsy
1982](#references)).  This is a setup file, it is not meant to be ran
on the wiki's server.  */
# define LEVEL 6
# include "grid/octree.h"
# define NPARTICLES 1
# define DLM_alpha_coupling 1

/* Stokes flow in an infinite array of spheres in a cubic domain */
# define diam 0.9508 
# define rhoval 1000. // density 
# define tval 60. // dynamical viscosity

/* output and numerical parameters */
# define Tc (diam) // convective velocity scale and we suppose Uc = 1
# define mydt  (Tc/50.)
# define maxtime (3*Tc)

/** Include the fictitious domain implementation */
# include "DLMFD_reverse_Uzawa.h"
# include "view.h"

/** File pointer to save the flowrate */
static FILE * flowpointer;

double deltau;
scalar un[];

int main() {
  L0 = 1;

  /* Stokes flow */
  stokes = true;
  
  /* set time step */
  DT = mydt;
   
  init_grid (1 << LEVEL);

  /* Tri-periodic cubic domain */ 
  foreach_dimension()
    periodic(top);

  /* Convergence criteria */
  TOLERANCE = 1e-5; 
  NITERMIN = 2;
  
  run();
}


event init (i = 0) {
  const scalar rhoc[] = rhoval;
  rho = rhoc;

  const face vector muc[] = {tval, tval, tval};
  mu = muc;

  origin (0., 0., 0.);

   /* Initialize acceleration (face) vectors for pressure gradient to drive the flow */
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
    boundary ((scalar *){a});
  }
  
  /* We set the initial fluid at rest. */
  foreach() {
    foreach_dimension()
      u.x[] = 0.;
    un[] = u.x[];
  }

  /* initial condition: particles/fictitious-domains positions */
  particle * p = particles;
  
  init_file_pointers (pdata, fdata, 0);

  /* Open file for storing the flowrate */
  flowpointer = fopen ("flowrate_right", "w");
  
  /**
The fictitious domain parameters are set here */
  for (int k = 0; k < NPARTICLES; k++) {

  /**
Sphere's geometric parameters are set here */
    GeomParameter gp = {0};
    gp.center.x = 0.5*L0;
    gp.center.y = 0.5*L0;
    gp.center.z = 0.5*L0;
    gp.radius = 0.5*diam;
    p[k].g = gp;
    
    p[k].iswall = 0;
    p[k].iscube = 0;
    p[k].pnum = k;
  }

}

/**
We log the number of iterations of the multigrid solver for pressure
and viscosity. */

event logfile (i++) {
  deltau = change (u.x, un);
  fprintf (stderr, "%d %g %d %d %g \n", i, t, mgp.i, mgu.i, deltau);
  if (t > maxtime) {
    return 1;
  }
}


/**
The acceleration event is used to impose a pressure gradient in
order to drive the flow. */

event acceleration(i++) {
  face vector av = a;
  coord dp = {10./L0/rhoval, 0, 0};
  foreach_face() {
    av.x[] = dp.x; 
  }
}

#if DLM_alpha_coupling
/**
 We overload the viscous_term event to the coupling of the two sub-problems */
event viscous_term (i++) {
  foreach() {
    foreach_dimension() {
      u.x[] += -dt*DLM_explicit.x[]/(dv()*rho[]);
    }
  }
  boundary((scalar*){u});
}
#endif

/**
We overload the adapt event to plug in the fictitious-domain problem */
event adapt (i++; t <= maxtime) {
  
  particle * pp = particles;
  DLMFD_subproblem (pp, i, rhoval);
   
  /* We save the forces acting on particles before adapting the mesh */
  sumLambda (pp, fdata, t, dt, flagfield, DLM_lambda, index_lambda, rhoval);

  /* We free particle's dlmfd points (we dont need them anymore) */
  free_particles (pp, NPARTICLES);
  
  /* We save all particles trajectories */
  particle_data(pp, t, i, pdata);

  /* We compute the flowrate at the "right" boundary */
  compute_flowrate_right(flowpointer, u, LEVEL);

}


event movie (t += Tc/24; t <= maxtime) {

  stats statsvelox;
  statsvelox = statsf (u.x);
  
  view (fov = 38.7802, quat = {-0.221938,0.254244,0.0441847,0.940295}, tx = -0.130528, ty = -0.0308596, bg = {0.3,0.4,0.6}, width = 724, height = 674, samples = 1);
  clear();
  squares("u.x", n = {0,0,1}, alpha = L0/2, map = cool_warm, min = statsvelox.min, max = statsvelox.max);
  squares("u.x", n = {1,0,0}, alpha = L0/2, map = cool_warm, min = statsvelox.min, max = statsvelox.max);
  cells(n = {1,0,0}, alpha = L0/2);
  cells(n = {0,0,1}, alpha = L0/2);
  box();
  save("movie.mp4")
}

event dump_last(t = end) {
  dump("dump");
}
 
/** 
# Solution and convergence rate

![Drag coefficient K vs dt](pics/zick_c0_45_K_vs_timestep.png)

![error vs dt](pics/zick_c0_45_error_vs_timestep.png)

*/

/**
# References 
~~~bib
@Article{zick1982,
  Title                    = {{Stokes flow through periodic arrays of spheres}},
  Author                   = {Zick, A.A. and Homsy, G.M.},
  Journal                  = {Journal of Fluid Mechanics},
  Year                     = {1982},
  Number                   = {1},
  Pages                    = {13--26},
  Volume                   = {115},

  File                     = {:JFM_115_13-26_1982.pdf:PDF},
  Publisher                = {Cambridge Univ Press},
  Timestamp                = {2012.08.20}
}
~~~
*/
