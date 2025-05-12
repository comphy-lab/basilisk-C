#include "adapt_wavelet_leave_interface.h"
#include "navier-stokes/centered.h"
#include "threephase.h"
#include "conserving-3f.h"
#include "view.h"
#include "tension.h"

#define JETTIME 2.5
#define VENTHALFWIDTH 1
#define VENTHEIGHT 2
#define BOXWIDTH 15.0
#define VEL myconst1

u.n[bottom] = dirichlet (6.5*sin(0.24)*(x > 1. && x < 2.));
u.t[bottom] = dirichlet (-6.5*cos(0.34)*(x > 1. && x < 2.)); //Are these right?! It doesn't look right in the simulation....

#define TRIANGLE (x < (-2./L0)*y)


#define RHOG (0.01) // Density 10 x greater than air. Approxmately the density of Radon gas (use this for viscoity)
#define RHOL (1.0)
#define RHOO (1.05) //Approximate 3rd phase as salt water. 
#define MUG (0.0008) // Used viscosity of Radon gas @ 600 degrees
#define MUL (0.02)
#define MUO (0.02)

#define SIGMAGL (0) // Assume at high temps surface tension is zero
#define SIGMALO (0)
#define SIGMAGO (0)

#define WIDTH 0.2
#define Hw 0 
#define MAXLEVEL 8


int main ()
{

  size (15.);
  origin (-L0/15.);
 
  rho1 = RHOG;
  rho2 = RHOL;
  rho3 = RHOO;
  mu1 = MUG;
  mu2 = MUL;
  mu3 = MUO;
  f1.sigma = (SIGMAGO+ SIGMAGL- SIGMALO)/2.;
  f2.sigma = (-SIGMAGO+ SIGMAGL+ SIGMALO)/2.;
  f3.sigma = (SIGMAGO- SIGMAGL+ SIGMALO)/2.;

  run();

}

event init (t = 0) {


   foreach() {
      if (x <= Hw ) {
	f1[] = 0;
	f2[] = 1.0;
	f3[] = 0;
	}
      if (x >= Hw) {
	f1[] = 1.;
	f2[] = 0.;
	f3[] = 0.;
	}
      if (x > 0. && x < 0.5 && y < 0.5 && y > 0.) {
	f1[] = 0.;
	f2[] = 0.;
	f3[] = 1.;
	}
	
    }
    boundary ({f1,f2,f3});
}



event fcorrect ( i++){
  double s=0.0;
  foreach(){
    s =(f1[]+f2[]+f3[]);
    f1[] =clamp((f1[]/s),0.,1.);
    f2[] =clamp((f2[]/s),0.,1.);
    f3[] =clamp((f3[]/s),0.,1.);
  }
  boundary ({f1,f2,f3});

}

scalar tri[];
event Stephanes_trick (i++) { // Sorry Stephane I like the name Stephane's Trick :P 
    fraction (tri, TRIANGLE);
    foreach(){
      foreach_dimension()
	u.x[] -= u.x[]*tri[];
    }
    //boundary (f1,f2,f3,(scalar *){u});
    boundary ((scalar *){u});
}



event logfile (i++) {
 
  double ke = 0., pe = 0., c = 0., d = 0.; //ignore this

  foreach(reduction(+:ke) reduction(+:pe)) {
   if (y > 6) {
   double r = rho(f1[], f2[], f3[]);
    ke += r*sq(norm(u))/2.*dv();
    pe += -1.*r*x*dv();
}
}
fprintf (ferr, "%g %g %g\n", t, ke, pe);
}


event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    //av.x[] -= 9.81;
      av.x[] -= 1.;
}

event outputs (i += 5) {
  dump();
}


event movie (t =0.; t += 0.03) {

view (fov = 24, quat = {0.00768262,0.0193538,-0.709094,0.704809}, tx = 0.525673, ty = -0.393655, bg = {1,1,1}, width = 600, height = 600, samples = 1);
  
  box();
  draw_vof ("f1", lw = 2, lc = {1,0,0});
  draw_vof ("f2", lw = 2, lc = {0,1,0});
  draw_vof ("f3", lw = 2, lc = {0,1,1});
  draw_vof ("tri");
  squares ("u.x", linear = true);
   save ("movie.mp4");

  clear();
  box();
  cells();
  save ("cells.mp4");
}


event adapt (i++) {
  adapt_wavelet({f1,f2,f3,u},(double[]){0.01,0.01,0.01,0.001,0.001}, MAXLEVEL);
}


event end (t = 10) {
  printf ("i = %d t = %g\n", i, t);
}

/**

## Results

![Representation of the different fluid fraction](injection_3f_test/movie.mp4)(width="800" height="600")


![Refinement of the mesh during](injection_3f_test/cells.mp4)(width="800" height="600")

*/