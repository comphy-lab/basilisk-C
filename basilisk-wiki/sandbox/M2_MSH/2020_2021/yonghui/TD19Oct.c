/**
## TD1 Gradiant des volumique
EX1 
x
EX2 
x
EX3 
x

*/
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
//#include "view.h"
# define LEVEL 8

/** 
## Boundary conditions */

u.n[left] = dirichlet(0.);
u.n[right] = dirichlet(0.);
u.n[top] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);

u.t[left] = neumann(0.);
u.t[right] = neumann(0.);
u.t[top] = neumann(0.);
u.t[bottom] = neumann(0.);
/**
p[left] = neumann(0.);
p[right] = neumann(0.);
p[top] = neumann(0.);
p[bottom] = neumann(0.);
*/
 
/** 
## Main */
FILE *fp; 
int main() {
  init_grid (pow(2,LEVEL));
	origin(-L0/2.,-L0/2.);

  rho1 = 1.;
  rho2 = 1.e-3;

  mu1 = 1.e-5;
  mu2 = 1.e-5;
  f.sigma = 0.2;

	//DT =0.001;

	fp =fopen ("ecine", "w");

 // TOLERANCE = 1e-6;
  run();
}
// compute the local normal of norm 2

coord normal (Point point, scalar c) {
  coord n = mycs (point, c);
  double nn = 0.;
  foreach_dimension()
    nn += sq(n.x);
  nn = sqrt(nn);
  foreach_dimension()
    n.x /= nn;
  return n;
}

// compute the local normal of norm 2 all over the domain

void compute_normal (scalar f, vector normal_vector) {
  foreach() {
    coord n = normal (point, f);
    foreach_dimension() 
      normal_vector.x[] = n.x;
  }
  boundary((scalar*){normal_vector});
}

	scalar kappa[];

event init (t = 0) {
	// generation of initial circular bubble
  fraction (f, -1.*(sq(x) + sq(y) - sq(0.2)) )  ;
	// output the shape of bubble 
  output_ppm (f, linear = false, max=1., min=0., n=512, file = "fold.png");	
  /**
  ## Curvature computation */
	boundary({f});
  curvature (f, kappa);
  vector n[];
  compute_normal (f, n);

  output_ppm (n.x, linear = false,  n=512, file = "nx.png");	
  output_ppm (n.y, linear = false,  n=512, file = "ny.png");	


  fraction (f, sq(0.2) - sq(x*2.) - sq(y/2.) )  ;
	// output the shape of elliptic bubble
  output_ppm (f, linear = false, max=1., min=0., n=512, file = "fnew.png");	
	foreach()
		foreach_dimension()
			u.x[] = 0.;
	
}

/**
# Mesures
*/
double ecine, ecap,etot;
event myossilation ( i += 10 ; t < 0.3 ){

  output_ppm (p, linear = false,  n=512, file = "pression.mp4");	
  output_ppm (f, linear = true,  n=512, file = "shape.mp4");	
  output_ppm (u.x, linear = false,  n=512, file = "ux.mp4");	
  output_ppm (u.y, linear = false,  n=512, file = "uy.mp4");	

	ecine = 0.;
	ecap = 0.;
	etot = 0.;
	foreach()
			ecine += 0.5*f[]*dv()*( sq(u.x[]) + sq(u.y[]) ) ;	

	double dss;
	dss = interface_area(f);
  //curvature (f, kappa);
	ecap = f.sigma*dss ;
	etot = ecap + ecine;

	fprintf (fp, "%g %g %g %g %g\n", t, dss, ecine,ecap,etot);	
	fprintf (stderr, "%d %g %g\n", i, t, dt);	

  output_ppm (f, linear = false, max=1., min=0., n=512, file = "fnow.png");	

}


/**
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){1e-3,1e-1,1e-1}, LEVEL,4);
}

*/
