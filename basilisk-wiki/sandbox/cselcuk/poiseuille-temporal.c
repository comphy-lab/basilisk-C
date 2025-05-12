/**
# Poiseuille flow using masking or DLMFD.
*/

#if DLMFD
# define NPARTICLES 2
# include "dlmfd.h"
#else
# include "navier-stokes/centered.h"
#endif

int main() {
  size (2.);
  origin (0, -L0/2.);
  stokes = true;
  TOLERANCE = 1e-5;
  
  const face vector g[] = {1.,0.};
  a = g;

  const face vector muc[] = {1.,1.};
  mu = muc;
  
  N = 64;
  DT = 1;
  
  for (int i = 1; i <= 4 ; i++) {
    DT /= 10;
    run();
  }
}

u.t[top] = dirichlet(0);
u.t[bottom] = dirichlet(0);
u.n[left] = neumann(0);
p[left] = dirichlet(0);
u.n[right] = neumann(0);
p[right] = dirichlet(0);


scalar un[];

event init (i = 0) {
#if !DLMFD
  mask (y > 0.5 ? top : y < -0.5 ? bottom : none);  
#endif
  
  foreach(){
    un[] = u.x[];
  }
#if DLMFD
  /* initial condition: particles position */
  particle * p = particles;
  
  init_file_pointers(pdata, fdata, 0);
  
  for (int k = 0; k < NPARTICLES; k++) {
    p[k].iswall = 1;
    coord wallmin, wallmax, wallpos;
    wallmin.x = 0;
    wallmax.x = L0;

    wallmax.y = -L0/4 + k*3*L0/4;
    wallmin.y = -L0/2 + k*3*L0/4;
    
    if (k == 0)
      wallpos.y = wallmax.y;
    if (k == 1)
      wallpos.y = wallmin.y;
    p[k].wallmin = wallmin;
    p[k].wallmax = wallmax;
    p[k].wallpos = wallpos;
  }
#endif
}



event logfile (t += DT*10; t <= 10) {
  double du = change (u.x, un);
  /* printf("#du = %10.8f, t = %f, dt = %f\n", du, t, dt); */
  if (i > 0 && du < 1e-6)
    return 1; /* stop */
}

event profile (t = end) {
  printf ("\n");
  scalar e[];
  foreach() {
    e[] = u.x[] - 0.5*(0.25 - y*y)*(fabs(y) < 0.5);
    printf ("%g %g %g %g %g %g %f\n", x, y, u.x[], u.y[], p[], e[], DT);//out
  }
  norm n = normf (e);
  fprintf (stderr, "%f %g %g %g\n", DT, n.avg, n.rms, n.max); // log

#if DLMFD
  /* performances display/output */
  /* particle * p = particles; */
  /* if (i > 10) output_dlmfd_perf (dlmfd_globaltiming, i, p); */
#endif
}

/** # Solution comparison
## Mask() function
~~~gnuplot Velocity profiles (mask)
plot [-0.5:0.5][0:0.14]                   \
     "< grep '0.100000$' out" u 2:3 t 'DT = 1e-1', \
     "< grep '0.010000$' out" u 2:3 t 'DT = 1e-2', \
     "< grep '0.001000$' out" u 2:3 t 'DT = 1e-3', \
     "< grep '0.000100$' out" u 2:3 t 'DT = 1e-4', \
     0.5*(0.25 - x*x)
~~~

~~~gnuplot Accuracy of the solution as a function of the level of refinement (mask)
set xlabel 'dt'
set ylabel 'Error norms'
set cbrange [1:1]
set logscale
ftitle(a,b) = sprintf("order %4.2f", b)
f2(x)=a2+b2*x
fit f2(x) 'log' u (log($1)):(log($3)) via a2,b2
fm(x)=am+bm*x
fit fm(x) 'log' u (log($1)):(log($4)) via am,bm
set xrange [0.0001:1]
plot exp (f2(log(x))) t ftitle(a2,b2), \
     exp (fm(log(x))) t ftitle(am,bm),  \
     'log' u 1:3 t '|e|_2' ps 1.5, \
     'log' u 1:4 t '|e|_{max}' ps 1.5 lc 0

~~~

## DLMFD
![Velocity profiles (DLMFD)](poiseuille-temporal-dlmfd/_plot0.png)

![Accuracy of the solution as a function the timestep (DLMFD)](poiseuille-temporal-dlmfd/_plot1.png)


## DLMFD-explicit forcing 
![Velocity profiles (DLMFD) with interior-oriented stencils](poiseuille-temporal-eforcing/_plot0.png)

![Accuracy of the solution as a function the timestep with explicit forcing (DLMFD)](poiseuille-temporal-eforcing/_plot1.png)

*/
