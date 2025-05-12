/**
# Triple point study test case

The domain is divided into three phases and they form a triple point (where three phases meet) at the center of the domain, we'll see the behavoir of the triple point. Specifically,

* Position (positioning of the triple point)
* Size (diameter of the triple point)
* Fractions error
$$
\epsilon_c=\frac{\int_V \left| \sum_i^3 c_i-1\right| dV}{V_0}.
$$

![Evolution of the three phases and the triple point (where three phases meet)](triple_point_study/movie.mp4)

* Infinity norm of the velocity
 */

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "three-phase.h"
#include "tension.h"
#include "norm-3p.h"
#include "view.h"

#define theta0 0
#define theta1 60
#define theta2 60
#define dtheta 10

#define MAXLEVEL 4

#define OH	0.25		

#define fluid1(x,y) max(x*tan((theta0+dtheta)*pi/180),x*tan((theta1+dtheta)*pi/180)) - y
#define fluid2(x,y) min(x*tan((theta0+dtheta)*pi/180),-x*tan((theta2-dtheta)*pi/180)) - y

/**

   We use the capillary-viscous velocity $U_\mu=\sigma_{12}/\mu_1$ because we don't have the carateristic length in this test case. The non-dimensional NS equation can be expressed as

$$   
\mathrm{Oh}^{-2}\rho\frac{\mathrm{D}\mathbf{u}}{\mathrm{D}t}=-\nabla p + \mu \nabla^2 \mathbf{u} + \sigma\kappa\delta\mathbf{n}
$$

 */
#define mu_1 1. 
#define mu_2 1.
#define mu_3 1.

#define rho_1 1. * sq(1./OH)
#define rho_2 1. * sq(1./OH)
#define rho_3 1. * sq(1./OH)

#define SIGMA12 1. // eau-air
#define SIGMA13 1. // huile-air
#define SIGMA32 1. // eau-huile

vector h1[], h2[], h3[];
FILE * fp, * fp1;

scalar f1e[], f2e[], f3e[];


h1.t[left]   = contact_angle (pi/2-dtheta*pi/180);
h1.t[top]    = contact_angle ((theta1+dtheta)*pi/180);
h2.t[left]   = contact_angle (pi/2+dtheta*pi/180);
h2.t[bottom] = contact_angle ((theta2-dtheta)*pi/180);
h3.t[top] = contact_angle (pi-(theta1+dtheta)*pi/180.);
h3.t[bottom]    = contact_angle (pi-(theta2-dtheta)*pi/180.);

scalar phi[]; 				//distance function used to define the fractions

int main(int argc, char * argv[]){
  origin(-0.5,-0.5);
  size(1);
  f1.height = h1;
  f2.height = h2;
  f3.height = h3;
  f1.sigma = (SIGMA13 + SIGMA12 - SIGMA32)/2.;
  f2.sigma = (SIGMA12 + SIGMA32 - SIGMA13)/2.;
  f3.sigma = (SIGMA32 - SIGMA12 + SIGMA13)/2.;
  rho1     =	rho_1;
  rho2     =	rho_2;
  rho3     =	rho_3;
  mu1      =	mu_1;
  mu2      =	mu_1;
  mu3      =	mu_1;

  N = 1 << MAXLEVEL;
  init_grid (N);
  char name[80];
  sprintf(name, "numerical_error");
  fp1 = fopen(name, "w");
  run();
  fclose(fp1);
  return 0;
}

event init(i=0) {
  if(!restore("dumpfilekk")){
    scalar f4[],f6[];
    fraction (f4,fluid1(x,y));
    fraction (f6,fluid2(x,y));

    foreach(){
      f1[] = 1-f4[];
      f2[] = f6[];
      f3[] = clamp(1.-f1[]-f2[],0.,1.);
      f1e[] = f1[];
      f2e[] = f2[];
      f3e[] = f3[];
    }
    boundary ((scalar *){f1, f2, f3, f1e, f2e, f3e, rhov});
    boundary ((scalar *){h1, h2, h3});
  }
  else{
    ;
  }
}

event movie(t+=10.){
  clear();
  box();
  cells();
  draw_vof("f1");
  draw_vof("f2", filled = 1, fc = {0.9, 0.4, 0.2});  
  draw_vof("f2");
  draw_vof("f3", filled = 1, fc = {0.2, 0.4, 0.9});
  draw_vof("f3");
  save("movie.mp4");
  coord p1;
  fprintf(fp1, "%g %g %g %g\n", t, locate_triple_point(&p1, f1, f2, f3, -L0/2., L0/2., -L0/2., L0/2., 1.2), fractions_error(f1, f2, f3, sq(L0)), shape_error(f1, f2, f3, f1e, f2e, f3e, sq(L0)));
}

event logfile (t=end){
  coord p1;
  locate_triple_point(&p1, f1, f2, f3, -L0/2., L0/2., -L0/2., L0/2., 1.2);
  output_cells(stdout);
  output_facets(f1, stderr);
  output_facets(f2, stderr);
  output_facets(f3, stderr);
  char filename[80];
  sprintf(filename, "triple_point");
  fp = fopen(filename, "w");
  fprintf(fp, "%g %g\n", p1.x, p1.y);
}


event end(t = 5000.){}

/**
## Results

### Diameter of the triple point

We can plot the evolution of the triple point zone "diameter"

~~~gnuplot Evolution of the triple point zone "diameter"
set xlabel 't'
set ylabel 'diameter'
plot 'numerical_error' u 1:2 w l notitle
~~~

The location of the center of the triple point is not good at the beginning of the simulation.

### Fractions error

and the fractions error

~~~gnuplot Evolution of the fractions error
set xlabel 't'
set ylabel 'fractions error'
plot 'numerical_error' u 1:3 w l notitle
~~~

The fractions error increases slowly along the simulation.

### Shape error evolution

~~~gnuplot Time evolution of the shape error
set xlabel 't'
set ylabel 'shape error'
plot 'numerical_error' u 1:4 w l notitle
~~~


### Triple point center position

~~~gnuplot Position of the triple point center
set size ratio -1
plot 'log' u 1:2 w l notitle,\
'out' u 1:2 w l notitle,\
'triple_point' u 1:2 w p lc 7 t 'triple point center'
~~~


 */
