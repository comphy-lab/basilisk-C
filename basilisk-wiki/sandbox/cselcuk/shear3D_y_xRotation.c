/**
# Rotation test x-direction: sphere freely rotating in a simple shear (Stokes) flow with DLMFD
*/

/** We measure the sphere's rotation rate in the x-direction by
imposition a simple shear along the y direction. The velocity field
reads $\vec{u} = \left(u,v,w\right) = \left(0,\dot{\gamma}z,0
\right)$, with $\dot{\gamma}$ the shear rate.  */ 

/** The surface velocity on the particle satisfies $\bm{\Omega}\cdot
\vec{x} = \vec{\omega^p}\times \vec{x}$, with $\bm{\Omega}$ being the
rotation rate tensor, $\vec{\omega^p}$ the (pseudo vector) `rotation
vector` and $\vec{x}$ the position vector of a point on the
sphere. Note that $\vec{\omega^p}$ is related to the vorticity
$\vec{\hat{\omega}}$ as $\vec{\omega^p} =
\frac{1}{2}\vec{\hat{\omega}}$.*/

/** Given $\vec{x} = \left(x,y,z\right)$, a point on the surface of the sphere and $\vec{\omega^p} = \left(\omega^p_x,\omega^p_y,\omega^p_z \right)$, it comes 

$$ \omega^p_y z = \omega^p_z y$$

$$\frac{1}{2}\dot\gamma z  = \omega^p_z x - \omega^p_x z$$ 

$$-\frac{1}{2}\dot\gamma y  = \omega^p_x y - \omega^p_y x$$ 

*/

/** for $\vec{x} = \left(D/2,0,0\right)$, it comes $\omega^p_z = \omega^p_y = 0$*/

/** for $\vec{x} = \left(0,D/2,0\right)$, it comes $$-\frac{1}{2}\dot{\gamma} = \omega^p_x$$ and $\omega^p_z = 0$ */

/** for $\vec{x} = \left(0,0,D/2\right)$, it comes $\omega^p_y = 0.$ which is redundent.   */

# define LEVEL 4
# include "grid/octree.h"
# define DLM_Moving_particle 1
# define TRANSLATION 0
# define ROTATION 1
# define DLM_alpha_coupling 1

# define NPARTICLES 1

# define adaptive 1
# define MAXLEVEL (LEVEL + 5)

/** # Physical parameters */
# define Uc 1. //caracteristic velocity
# define rhoval 1. // fluid density
# define diam (1.) // particle diameter
# define ReD  (0.001) // Reynolds number based on the particle's diameter to setupd the viscosity
# define Ldomain (20.)
# define rhosolid 2. //particle density
# define tval (rhoval*Uc*diam/ReD) // fluid density
 

/** # Output and numerical parameters */
# define Tc (diam/Uc) // caracteristic time scale
# define mydt (Tc/200.) // maximum time-step
# define maxtime (1.)
# define tsave (Tc/1.)

/** 
We include the ficitious-domain implementation 
*/

# include "dlmfd.h"
# include "view.h"

double deltau;
scalar un[];

int main() {
  L0 = Ldomain;

  stokes = true;
  
  /* set time step */
  DT = mydt;
     
  /* initialize grid */
  init_grid(1 << (LEVEL));

  /* boundary conditions */
  /** The shear rate is $\dot{\gamma} = 2Uc/L0$*/

  periodic(left);
  periodic(top);
  
  /* from boundary */
  uf.n[front] = dirichlet(0.);
  uf.r[front] = dirichlet(Uc);
  uf.t[front] = dirichlet(0.);

  u.n[front] = dirichlet(0.);
  u.r[front] = dirichlet(Uc);
  u.t[front] = dirichlet(0.);
  p[front] = neumann(0.);
  pf[front] = neumann(0.);
    
  /* back boundary */
  uf.n[back] = dirichlet(0.);
  uf.r[back] = dirichlet(-Uc);
  uf.t[back] = dirichlet(0.);
  u.n[back] = dirichlet(0.);
  u.r[back] = dirichlet(-Uc);
  u.t[back] = dirichlet(0.);
  p[back] = neumann(0.);
  pf[back] = neumann(0.);
  
  
  /* Convergence criteria */
  TOLERANCE = 1e-4;

  run();
}





/**
We initialize the fluid and particle variables. */

event init (i = 0) {

  /* set origin */
  origin (0., 0., 0.);

  
  if (!restore (file = "dump")) {
    /* fluid initial condition: */
    foreach() {
      u.y[] = -Uc + 2*z*Uc/L0;
      un[] = u.x[];
    }
    /* initial condition: particles position/velocity */
    particle * p = particles;

    for (int k = 0; k < NPARTICLES; k++) {
      GeomParameter gp = {0};
    
      p[k].iscube = 0;
      p[k].iswall = 0;
    
      gp.center.x = L0/2.;
      gp.center.y = L0/2.;
      gp.center.z = L0/2.;
      gp.radius   = diam/2.;
   
      p[k].g = gp;
   
      /* initial condition: particle's velocity */
      coord c = {0., 0., 0.};
      p[k].w = c;
    
    }
  }
  else { // restart run, the default init event will take care of it
  }
}

/**
We log the number of iterations of the multigrid solver for pressure
and viscosity. */

event logfile (i++) {
  deltau = change (u.x, un);
  fprintf (stderr, "log output %d %g %d %d %g %g %g %ld\n", i, t, mgp.i, mgu.i, mgp.resa, mgu.resa, deltau, grid->tn);

}



event output_data (t += 0.01; t < maxtime) {
  stats statsvelox;
  /* scalar omega[]; */
  view (fov = 22.3366, quat = {1,0,0,1}, tx = -0.465283, ty = -0.439056, bg = {1,1,1}, width = 890, height = 862, samples = 1);
  /* vorticity(u, omega); */
  statsvelox =  statsf (u.y);
  clear();
  squares ("u.y", n = {1,0,0}, alpha = L0/2, map = cool_warm, min = statsvelox.min, max = statsvelox.max); 
  cells(n = {1,0,0}, alpha = L0/2);
  save ("movie.mp4");
}

event last_output (t=end){
  p.nodump = false;
  dump("dump_bc");
}
/**

# Results

~~~gnuplot particle's angular velocities $\omega^p_z$, $\omega^p_y$
set grid
show grid
plot "particle-data-0" u 1:10 w l title "\omega^p_z", "particle-data-0" u 1:9 w l title "\omega^p_y", 0 title "analytical omega^p_z",0 title "analytical omega^p_y"
~~~

~~~gnuplot particle's angular velocity $\omega^p_x$ 
set grid
show grid
Uc = 1.;
L = 20.;
shear = 2.*Uc/L
set yrange [-0.055:-0.025]
plot  "particle-data-0" u 1:8 w l title "\omega^p_x",-shear*0.5 title "analytical omega^p_x"
~~~

~~~gnuplot Torque on the particle $T_x, T_y, T_z$
reset
set grid
show grid
plot "sum_lambda-0" u 1:5 w l title "T_x", "sum_lambda-0" u 1:6 w l title "T_y","sum_lambda-0" u 1:7 w l title "T_z"
~~~
*/

