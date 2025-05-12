/**
# Purely rotating interface

VOF advection doesn't cope well with rotations. Here we test a correction on the
advection.

[Same example with Gerris](http://http://gerris.dalembert.upmc.fr/gerris/tests/tests/rotate.html)

~~~gnuplot error plot central cell
set logscale y
plot 'err_reference' w l t 'classic VOF', 'log' w l t 'rotation corrected'
~~~

~~~gnuplot Results
unset logscale y
set size ratio -1
set style line 2 dt 2 lw 1.5 lc rgb 'black'
f(x,t) = x/t
plot [-0.5:0.5][-0.5:0.5]'out' w l t 'interface', f(x,1) ls 2 t 'theory' , f(x,2) ls 2 t '', f(x,5) ls 2 t '',\
                         'VOF_classic' w l t 'std VOF' lc rgb 'blue'
~~~



*/
double capteur(Point point, scalar s){
  return s[];
}

@define str(x) #x

#include "advection.h"
#include "../vof_rotation.h"
// #include "vof.h"
#include "view.h"

scalar f[],f2[];
scalar * interfaces = {f}, * tracers = NULL;

#define Pi 3.14159265358979323846
#define T 5



char filename [100];
FILE * fp1;



int main()
{
  L0 = 2;
  origin (-L0/2., -L0/2.);
  // snprintf(filename, 100,  "reference");
  fp1 = fopen (filename,"w");
  /**
  We then run the simulation for different levels of refinement. */

  init_grid (1 << 4);
  CFL = 0.45;
  DT = CFL*L0/(1 << grid->maxdepth);

  run();
  fclose(fp1);
}


event init (i = 0)
{
  fraction (f, x);
  foreach_face(x)
    u.x[] = y;
  foreach_face(y)
    u.y[] = 0;
  u.t[top]    = dirichlet(y);
  u.t[bottom] = dirichlet(y);
  u.n[left]   = dirichlet(y);
  u.n[right]  = dirichlet(y);
  boundary ((scalar *){u});
  event("velocity");
  
}

event check_volume (i++,last) {
  fraction(f2, x/sqrt(1+t*t) - t*y/sqrt(1+t*t)); // normal to the interface
  double myDelta = L0/(1 << grid->maxdepth);
  Point p = locate(myDelta,0);
  double val1 = capteur(p,f);
  double val2 = capteur(p,f2);
  fprintf(stderr, "%g %g\n", t, fabs(val1- val2));
}


event movie(i++){
  draw_vof("f");
  save("f.mp4");
}

event shape (t ={1,2,5}) {
  output_facets (f,stdout); 
}

event logfile(t  = {0, T},last){
  fprintf(stderr, "##%d %g\n", i, t);
  if(t==T)exit(1);
}


/**

*/

