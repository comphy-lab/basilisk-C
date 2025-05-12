/**
# Flow past a sphere with DLMFD
*/

/** Flow past a sphere with DLMFD */

# define LEVEL 6
# include "grid/octree.h"
# define DLM_Moving_particle 0
# define adaptive 1
# define NPARTICLES 1

# if adaptive
# define MAXLEVEL (LEVEL + 1)
# endif

/** 
# Physical parameters 
*/
# define diam (1.) 
# define Um 1
# define rhoval 1
# define Re 50
# define fs_density_ratio  2. // fluid solid density ratio
# define Ld_ratio 8. // box size-particle diameter ratio
# define Ldomain (Ld_ratio*diam)
# define rhosolid (fs_density_ratio*rhoval) //particle density
# define tval (Um*rhoval*diam)/Re

/** 
# Output and numerical parameters 
*/

# define Tc (diam/Um) // caracteristic time scale
# define mydt (Tc/100) // maximum time-step (the time step is also adaptive in time but it won't exceed this value)
# define maxtime (2.)
# define tsave (Tc/200.)

/** 
We include the ficitious-domain implementation 
*/

# include "DLMFD_reverse_Uzawa.h"
# include "view.h"
# include "lambda2.h"

 
double deltau;
scalar un[];

int main() {		
  L0 = Ldomain;

  // set time step
  DT = mydt;
  
  /* initialize grid */
  init_grid(1 << LEVEL);

  /* left boundary */
  u.n[left] = dirichlet(Um); 
  u.r[left] = dirichlet(0.); 
  u.t[left] = dirichlet(0.); 
  
  /* right boundary */
  u.n[right] = neumann(0.);   
  u.r[right] = dirichlet(0.);
  u.t[right] = dirichlet(0.);
 
  /* top boundary */
  u.n[top] = dirichlet(0.);
  u.r[top] = dirichlet(Um);  
  u.t[top] = dirichlet(0.);  

  /* bottom boundary */
  u.n[bottom] = dirichlet(0.); 
  u.r[bottom] = dirichlet(Um);  
  u.t[bottom] = dirichlet(0.);
  
  /* front boundary */
  u.n[front] = dirichlet(0.);  
  u.r[front] = dirichlet(0.);   
  u.t[front] = dirichlet(Um);   

  /* back boundary */
  u.n[back] = dirichlet(0.);  
  u.r[back] = dirichlet(0.);  
  u.t[back] = dirichlet(Um);   
  
  p[left]    = neumann(0.);
  p[right]   = dirichlet(0.);

  /* Convergence criteria */
   TOLERANCE = 1e-3; 
  
  run();
}


event init (i = 0) {
  /* set origin */
  origin (0, 0, 0);

  /* set dynamic viscosity */
  const face vector muc[] = {tval,tval,tval};
  mu = muc;

  /* set density of the flow */ 
  const scalar rhoc[] = rhoval;
  rho = rhoc;


  /* We set the initially horizontal velocity to unity everywhere.  */
  foreach(){
    u.x[] = Um;
    un[] = u.x[];
  }

  /* initial condition: particles position */
  particle * p = particles;
  
  init_file_pointers (pdata, fdata, 0);
  
  for (int k = 0; k < NPARTICLES; k++) {
    GeomParameter gp = {0};
    gp.center.x = L0/4;
    gp.center.y = L0/2;
    gp.center.z = L0/2;
    gp.radius = diam/2;

    p[k].g = gp;
   
    /* particle id */
    p[k].pnum = k;
    
    /* density rho_s of the particle */
    p[k].rho_s = rhosolid;

    /* Volume or surface of the particle (circle of sphere) */
    p[k].Vp = 4*pi*pow(gp.radius,3)/3;

    /* total weight of the particle */
    p[k].M = rhosolid*(p[k].Vp);

#if DLM_Moving_particle
    /* The inertia tensor is: */
    /* When in doupt check https://en.wikipedia.org/wiki/List_of_moments_of_inertia */
    /* Ip[0] = Ixx */
    /* Ip[1] = Iyy */
    /* Ip[2] = Izz */
    /* Ip[3] = Ixy */
    /* Ip[4] = Ixz */
    /* Ip[5] = Iyz */
    
    /* For a solid disk: */
    p[k].Ip[0] = 2*(p[k].M)*sq(gp.radius)/5;
    p[k].Ip[1] = p[k].Ip[0];
    p[k].Ip[2] = p[k].Ip[0]; 
    p[k].Ip[3] = 0.;
    p[k].Ip[4] = 0.;
    p[k].Ip[5] = 0.;
    

    /* initial condition: particle's velocity */
    coord c;
    foreach_dimension()
      c.x = 0;
#if TRANSLATION
    p[k].U = c;
#endif
#if ROTATION
    p[k].w = c;
#endif
   p[k].wished_ratio = 0.1;
   p[k].en = 1;
   p[k].vzero = 1;
   compute_wo (p);
#endif

  }
}




/**
We log the number of iterations of the multigrid solver for pressure
and viscosity. */

event logfile (i++) {
  deltau = change (u.x, un);
  fprintf (stderr, "%d %g %d %d %g \n", i, t, mgp.i, mgu.i, deltau);
}

#if DLM_alpha_coupling
event viscous_term (i++) {
  foreach() {
    foreach_dimension()
      u.x[] += dt*DLM_lambda.x[]/(rho[]*dv());
  }
}
#endif


event adapt (i++) {
  particle * pp = particles;
  int totalcell = totalcells ();
  if (pid() == 0) {
    printf ("total cells = %d\n", totalcell);
    printf ("starting dlmfd subproblem \n");
  }
  DLMFD_subproblem (pp, i, rhoval);
   
  /* Save forces acting on particles before the adapting the mesh */
  sumLambda (pp, fdata, t, dt, flagfield, DLM_lambda, index_lambda, rhoval);

  /* Free particle structures (we dont need them anymore) */
  free_particles (pp, NPARTICLES);
   
  /* Save particles trajectories */
  particle_data (pp, t, i, pdata);
   
#if adaptive
#if DLM_Moving_particle  
  astats s = adapt_wavelet ((scalar *){flagfield_mailleur, u}, (double[]){1e-4, 0.01,0.01,0.01}, maxlevel = MAXLEVEL);
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);

#elif DLM_Moving_particle == 0
  
  astats s = adapt_wavelet ((scalar *){flagfield, u}, (double[]){1e-4,1e-2,1e-2,1e-2}, maxlevel = MAXLEVEL);
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
#endif
#endif
}


event output_data (t+=maxtime/50; t < maxtime) {

  /* char name[80] = "flowsphere"; */
  /* scalar * list = {p, flagfield}; */
  /* vector * vlist = {u, index_lambda}; */
  /* save_data (list, vlist, i, t, name); */
  
  view (fov = 39.6758, quat = {-0.200394,0.310041,0.0701111,0.926715}, tx = -0.0703854, ty = -0.0351928, bg = {0.3,0.4,0.6}, width = 600, height = 600, samples = 1);

  stats statsvelox;
  statsvelox = statsf (u.x);
  clear ();
  box ();
  cells (n = {0,0,1}, alpha = L0/2);
  squares ("u.x", n = {0,0,1}, alpha = L0/2, map = cool_warm, min = statsvelox.min, max = statsvelox.max);
  squares ("u.x", n = {0,1,0}, alpha = L0/2, map = cool_warm, min = statsvelox.min, max = statsvelox.max);
  cells (n = {0,1,0}, alpha = L0/2);
  scalar l2[];
  lambda2 (u, l2);
  isosurface ("l2", -0.0002);
  save ("movie.mp4");
}
 
event lastdump (t = maxtime) {
  dump (file = "dump");
}


/**
# Animation

## Lambda_2 criterion and streamwise velocity u.x
<video width="600" height="600" controls>
<source src="sphere/movie.mp4" type="video/mp4">
</video>
*/
