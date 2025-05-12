/**
# Mass flow 3D test 
*/
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "view.h"
#define tmax 1.*2.*M_PI
#define ADAPT 0 

/**
##The main function . */
int LEVEL = 4;
int MAX_LEVEL = 6;
double Pout;
double DEBIT;
double QQQ;

FILE * fp1;
FILE * fp2;
FILE * fp3;
scalar f0[];
int main() {
  //We set the domain geometry and initial refinement.
  size (1);
  origin (0.,-L0/2., -L0/2.);
  fp1 = fopen ("testpoint", "w");
  fp2 = fopen ("pressure.dat", "w");
  init_grid (32);
  run();
}

//normal pressure condition
pf[left]=neumann(0.);
p[left]=neumann(0.);

//normal flow at two side
u.n[left] = dirichlet(1.*f0[]*cos(t));
u.n[right] = neumann(0.);
u.t[right] =  dirichlet(0.);
u.t[left] = dirichlet(0.);
u.r[right] =  dirichlet(0.);
u.r[left] = dirichlet(0.);


event init (t = 0) {
  fraction(f0, sq(0.4) - sq(z) - sq(y));
  const face vector muc[] = {0.1, 0.1, 0.1};
  mu = muc;
  view(fov = 30.65,tx= - L0/2., width = 640, height = 640);	
  clear();
  box();
  draw_vof("f0", edges = true, lw = 0.5);
  save("f0.png");
}

event midpressure(t <= tmax; i += 1) {
  //BC : lateral velocity = 0
  foreach()
    foreach_dimension()
    u.x[] = f0[]*u.x[];
    QQQ =0.;
  DEBIT = 0.;
  // we calculate the mass flow 
  foreach_boundary(right)
    DEBIT += u.x[]*sq(Delta)*f0[];
  Pout = QQQ * 0.2;
  // idea 2 
  for (double yy = -0.5; yy < 0.5; yy += L0/100. ){
    for (double zz = -0.5; zz < 0.5; zz += L0/100. ){
      QQQ += interpolate(u.x, 0.9, yy, zz)*sq(L0/100.);
    }	
  }
  pf[right]=dirichlet(Pout);
  p[right]=dirichlet(Pout);
  //output video and datas  
  fprintf(fp1,"%g %g %g %g\n" , t, interpolate(u.x, L0/2, 0.,0.), DEBIT, QQQ);
  view (fov = 30.65,tx= - L0/2., width = 640, height = 640);
  box();
  squares("u.x",min=-1.5,max=1.5);
  cells();
  save("ux.mp4");
}


/**
We add the acceleration of gravity (unity) in the downward (-y)
direction. */
event tracer (t = end) {
  for (double xx = 0.; xx < 1.; xx += L0/100. ){
    fprintf(fp2,"%d %g %g %g %g\n", N, t, xx, interpolate(p , xx, 0.,0.),Pout);
  }
}


/**
We adapt the mesh by controlling the error on the volume fraction and
velocity field. */
#if ADAPT
event adapt (i++) {
  double uemax = 1e-2;
  //  adapt_wavelet ({f,u}, (double[]){5e-4,uemax,uemax,uemax}, MAX_LEVEL, 5);
  adapt_wavelet ({f,u}, (double[]){5e-3,uemax,uemax,uemax}, LEVEL, 2);
}
#endif

/**
# Results
~~~gnuplot testpoint
plot 'testpoint' u 1:3 w l t'DEBIT' ,\
'testpoint' u 1:4 w l t'QQQ' ,\
0.50265482457*cos(x) t 'Q theory = 0.5'
~~~

~~~gnuplot pressure at y = L0/2
plot 'pressure.dat' u 3:4 w l t'pressure' ,\
'pressure.dat' u 3:5 w l t'pressure out'
~~~
*/
