/** 
# Sphere with initial rotation rate and subject to gravity in a box of L=120D
*/

# include "grid/octree.h"
# define DLM_Moving_particle 1
# define NPARTICLES 1
# define DLM_alpha_coupling 1



/** 
# Parameters related to the spatial resolution 
*/
# define LEVEL 4 
# define adaptive 1
# if adaptive
# define MAXLEVEL (LEVEL + 8)
# endif

/** 
# Physical parameters 
The characteristic velocity scale is $U_c = \sqrt(gD)$ with D the sphere's diameters (characteristic length scale).The Reynolds number based on $U_c$ is Re = 25.
*/
# define rhoval 1000. // fluid density
# define rhosolid 1140.
# define gravity_y -9.81 // gravity acceleration
# define Lc (0.01)  // characteristic length - particle diameter
# define tval 0.125284 // fluid dynamic viscosity, here Re = 25
# define Uc 0.31321 // sqrt(grav*Lc)
# define Tc (Lc/Uc) // characteristic time

/* output and numerical parameters */
# define mydt (Tc/50) //time-step
# define maxtime (5*Tc) // simulation time


/** 
# Fictitious domain implementation 
*/
# include "dlmfd-toygs.h"
# include "navier-stokes/perfs.h"


double deltau;
scalar un[];

int main () {
  
  origin (0., 0., 0.);
 
  L0 = 1.2;
   
  /* set time step */
  DT = mydt;

  /* initialise a uniform grid */
  init_grid (1 << LEVEL);

  /* boundary conditions */
  u.t[left]   = dirichlet(0.); //v
  u.r[left]   = dirichlet(0.); //w
  u.n[left]   = dirichlet(0.); //u
  uf.t[left]   = dirichlet(0.); //v
  uf.r[left]   = dirichlet(0.); //w
  uf.n[left]   = dirichlet(0.); //u

  u.t[right]  = dirichlet(0.); //v
  u.r[right]  = dirichlet(0.); //w
  u.n[right]  = dirichlet(0.); //u
  uf.t[right]  = dirichlet(0.); //v
  uf.r[right]  = dirichlet(0.); //w
  uf.n[right]  = dirichlet(0.); //u

  
  u.n[bottom] = dirichlet(0.); //v
  u.t[bottom] = dirichlet(0.); //w
  u.r[bottom] = dirichlet(0.); //u
  uf.n[bottom] = dirichlet(0.); //v
  uf.t[bottom] = dirichlet(0.); //w
  uf.r[bottom] = dirichlet(0.); //u
 
  u.n[top]    = dirichlet(0.); //v
  u.t[top]    = dirichlet(0.); //w
  u.r[top]    = dirichlet(0.); //u
  uf.n[top]    = dirichlet(0.); //v
  uf.t[top]    = dirichlet(0.); //w
  uf.r[top]    = dirichlet(0.); //u
  
  u.n[front]  = dirichlet(0.); //w
  u.t[front]  = dirichlet(0.); //u
  u.r[front]  = dirichlet(0.); //v
  uf.n[front]  = dirichlet(0.); //w
  uf.t[front]  = dirichlet(0.); //u
  uf.r[front]  = dirichlet(0.); //v

  u.n[back]   = dirichlet(0.); //w
  u.t[back]   = dirichlet(0.); //u
  u.r[back]   = dirichlet(0.); //v
  uf.n[back]   = dirichlet(0.); //w
  uf.t[back]   = dirichlet(0.); //u
  uf.r[back]   = dirichlet(0.); //v

  
  // Convergence criteria for the N-S solver
  TOLERANCE = 1e-4; 
  NITERMIN = 2;
  run();
}

event init (i = 0) {
  
  particle * p = particles;

  /* Initial condition for the fluid */
  foreach() {
    foreach_dimension() {
      u.x[] = 0.;
    }
    un[] = u.x[];
  }

  double diam = Lc;
  /* Set the particle types (spheres and cubes only for now) */
  for (int k = 0; k < NPARTICLES; k++) {   

    GeomParameter gp = {0};
    gp.center.x = L0/2;
    gp.center.y = L0-20*diam;
    gp.center.z = L0/2;
    gp.radius = diam/2.;
    
    p[k].g = gp;
    
    /* initial condition: particle's velocities */
    coord c = {0., 0., 0. };
    p[k].U = c;
    
    c.z = 1.;
    p[k].w = c; 
  }    
}


event output_data (t += Tc; t <= maxtime)
{

  dump(file = "dump");
  particle * p = particles;
  dump_particle(p);
  /* performances display/output */
  if (i > 10) output_dlmfd_perf (dlmfd_globaltiming, i, p);
  
}
 event last_output (t=end){
  p.nodump = false;
  dump("dump");
  
  /* performances display/output */
  particle * p = particles;
  if (i > 10) output_dlmfd_perf (dlmfd_globaltiming, i, p);
}


/**

# Results

~~~gnuplot dimensionless particle angular velocities $t_c\omega^p_x$, $t_c\omega^p_y$ , $t_c\omega^p_z$
set grid
show grid
Lc = 0.01
Uc = 0.31321
tc = Lc/Uc
set xlabel "t/t_c"
set ylabel "t_c*\omega"
plot "particle-data-0" u ($1/tc):($8*tc) w l title "t_c\omega^p_x", "particle-data-0" u ($1/tc):($9*tc) w l title "\omega^p_y t_c", "particle-data-0" u ($1/tc):($10*tc) w l title "\omega^p_z t_c"
~~~

~~~gnuplot dimensionless particle translational velocities $U_x/Uc$, $U_y/Uc$, $U_z/Uc$ 
reset
set grid
show grid
Lc = 0.01
Uc = 0.31321
tc = Lc/Uc
set xlabel "t/t_c"
set ylabel "u/Uc"
plot  "particle-data-0" u ($1/tc):($5/Uc) w l title "U_x/Uc","particle-data-0" u ($1/tc):($6/Uc) w l title "U_y/Uc", "particle-data-0" u ($1/tc):($7/Uc) w l title "U_z/Uc"
~~~
~~~gnuplot particle trajectory in x-y plane at $z=L0/2$
reset
diam = 0.01
show grid
set xlabel "x/diam"
set ylabel "y/diam"
plot "particle-data-0" u ($2/diam):($3/diam) w l title "trajectory"
~~~
*/
