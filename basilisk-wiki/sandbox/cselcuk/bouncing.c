/**
# Bouncing motion of a sphere at Stokes = 58 with a fictitious-domain method
*/

/** 
We study the motion of a sphere falling under gravity in a viscous
fluid.  This is a setup/example file, it is not meant to be ran on the
wiki. We compare our results with respect to [Ardekani et al. 2008](#references) and [Gondret et al. 2002](#references).
*/
# define LEVEL 6
# include "grid/octree.h"
# define DLM_alpha_coupling 1
# define DLM_Moving_particle 0
# define NPARTICLES 1
# define adaptive 1

# if adaptive
# define MAXLEVEL (LEVEL + 1)
# endif

# define diam 1. 
# define rhoval 1.
# define tval 100.
# define Uc 1.
# define Tc (diam/Uc)

/* output and numerical parameters */
# define mydt  (0.00001)
# define maxtime (20.*Tc)
# include "DLMFD_reverse_Uzawa.h"

double deltau;
scalar un[];

int main() {
  L0 = 10*diam;
  stokes  = true;
  
  // set time step
  DT = mydt;
   
  init_grid(1 << LEVEL);

  /* boundary conditions */
  u.n[right] = dirichlet(0.);
  u.t[right] = neumann(0.);

  u.n[left] = dirichlet(0.);
  u.t[left] = neumann(0.);  

  u.n[top] = dirichlet(0.);
  u.t[top] = dirichlet(0.);

  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);

  u.r[bottom] = dirichlet(0.);
  u.r[top] = dirichlet(0.);
  u.r[left] = neumann(0.);
  u.r[right] = neumann(0.);
  
  u.n[front]  = dirichlet(0.);
  u.t[front]  = dirichlet(0.);
  u.r[front]  = dirichlet(0.);
  
  u.n[back]   = dirichlet(0.);
  u.t[back]   = dirichlet(0.);  
  u.r[back]   = dirichlet(0.);
  
  p[front]  = neumann(0.);
  p[back]   = neumann(0.);
  p[left]   = neumann(0.);
  p[right]  = neumann(0.);
  p[bottom] = neumann(0.);
  p[top]    = neumann(0.); 
  
    
  // Convergence criteria
  TOLERANCE = 1e-5; 
  
  run();
}


event init (i = 0) {
  particle * p = particles;
  
  const scalar rhoc[] = rhoval;
  rho = rhoc;

  const face vector muc[] = {tval, tval, tval};
  mu = muc;

  origin (0., 0., 0.);

  if (!restore (file = "dump")) {
    /* Initial condition: fluid at rest. */
    foreach() {
      foreach_dimension() {
	u.x[] = 0.;
      }
      un[] = u.y[];
    }

    /* Initial condition: particles/fictitious-domains positions */
    init_file_pointers(pdata, fdata, 0);

    for (int k = 0; k < NPARTICLES; k++) {
    
      GeomParameter gp;
      p[k].iswall = 0;
      gp.center.x = 0.5*L0;
      gp.center.y = 1.025*0.5*diam;
      gp.center.z = 0.5*L0;
      gp.radius = 0.5*diam;
      p[k].g = gp;
    
      p[k].pnum = k;
      (p[k]).imposedU.y = -Uc;
    }
  }
  else { 
    restore_particle (p);
    
    /* Initialize file poiters */
    init_file_pointers(pdata, fdata, 1);
    fprintf (ferr, "simulation restored: it has to go for t_end = %f\n", maxtime);
  }
}

/**
We log the number of iterations of the multigrid solver for pressure
and viscosity. */

event logfile (i++) {
  deltau = change (u.y, un);
  fprintf (stderr, "%d %g %d %d %g \n", i, t, mgp.i, mgu.i, deltau);
  if (t > maxtime){
    return 1;
  }
}


#if DLM_alpha_coupling
event viscous_term (i++) {
  foreach() {
    foreach_dimension() {
      u.x[] += -dt*DLM_explicit.x[]/(rho[]*dv());
    }
  }
  boundary((scalar*){u});
}
#endif

event adapt (i++) {
  
  particle * pp = particles;
  
  DLMFD_subproblem (pp, i, rhoval);
   
  /* We save the forces acting on particles before adapting the mesh */
  sumLambda (pp, fdata, t, dt, flagfield, DLM_lambda, index_lambda, rhoval);

  /* We free particle's dlmfd points (we dont need them anymore) */
  free_particles (pp, NPARTICLES);
  
  /* We save all particles trajectories */
  particle_data(pp, t, i, pdata);


#if adaptive
  astats s = adapt_wavelet ((scalar *){flagfield, u}, (double[]){1e-4, 1e-3,1e-3,1e-3}, maxlevel = MAXLEVEL);
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
#endif
}

event output_data (i += 100; t <= maxtime) {
  dump(file = "dump");
  particle * p = particles;
  dump_particle(p);
}
 
event lastdump (t = end) {
  dump(file = "dump");
  particle * p = particles;
  dump_particle(p);
}
/**
# Results

## Temporal convergence

![Sphere's center position versus time for various time-steps](pics/stokes_58_sphere_pos_vs_time.png)

![Sphere's velocity versus time for various time-steps](pics/stokes_58_sphere_velo_vs_time.png)

Cannot obtain a clear convergence on this case.

## Spatial convergence

![Sphere's center position versus time for various spatial resolutions](pics/stokes_58_sphere_pos_vs_time_spatial_convergence.png)

![Sphere's velocity versus time for various spatial resolutions](pics/stokes_58_sphere_velo_vs_time_spatial_convergence.png)

## Annimation
<video width="1280" height="720" controls>
<source src="http://www.basilisk.fr/sandbox/cselcuk/pics/movie_St_58_Re_63.mp4" type="video/mp4">
</video>

# References

~~~bib
@article{ardekani2008numerical,
  title={Numerical investigation of particle--particle and particle--wall collisions in a viscous fluid},
  author={Ardekani, AM and Rangel, RH},
  journal={Journal of fluid mechanics},
  volume={596},
  pages={437--466},
  year={2008},
  publisher={Cambridge University Press}
}

@article{gondret2002bouncing,
  title={Bouncing motion of spherical particles in fluids},
  author={Gondret, P and Lance, M and Petit, L},
  journal={Physics of fluids},
  volume={14},
  number={2},
  pages={643--652},
  year={2002},
  publisher={AIP}
}
~~~
*/
