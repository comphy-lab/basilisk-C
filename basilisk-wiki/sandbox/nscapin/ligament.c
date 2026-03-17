
#include "axi.h"
#include "navier-stokes/centered.h"
#if CLSVOF
#include "two-phase-clsvof.h"
#include "integral.h"
#include "curvature.h"
#else
#include "two-phase.h"
#include "tension.h"
#endif
#include "tag.h"
#include "view.h"

//Physical ratios and box size
#define RHOR 850.0 // density ratio
#define MUR 25.0   // kinematic viscosity ratio (real value 55; to be consistent with turbubble simulations)
#define HALFBOX 5.

/** 
   Input parameters of the simulations.  */

double Oh  = 0.015;   // Ohnesorge number
double thk = 0.1;     // thickness of the film
double mu1_in = 0.01; // input viscosity
int MAXLEVEL = 9;     // max level of the simulation
int MINLEVEL = 6;     // min level of the simulation
double MAXTIME = 1;   // maximum simulation time
int istart = 1;       // start

/** 
   We define some auxiliary variables. */

double SIGMA = 1.; // we change later
double P1x = 1., P1y = 1.; // we change later
double P2x = 1., P2y = 1.; // we change later
double Lx = 1.;

/**
   For the mesh refinement. */

double uemax = 1.e-3;
double femax = 1.e-4;

/** 
   Boundary conditions are chosen to match a zero-pressure outflow. */

u.n[left]   = dirichlet(0.);
u.n[right]  = neumann(0);
u.n[top]    = neumann(0);
u.n[bottom] = dirichlet(0.);

p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);
p[bottom] = neumann(0.);
pf[bottom] = neumann(0.);

p[left] = neumann(0.);
pf[left] = neumann(0.);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);

int main(int argc, char * argv[]){
  
  fprintf(stderr, "************************\n"), fflush (stderr);
  fprintf(stderr, "Check of input parameters only\n"), fflush (stderr);
  fprintf(stderr, " Ohnesorge number = %.10e\n ", Oh), fflush (stderr);
  fprintf(stderr, " Thickness = %.10e\n ", thk), fflush (stderr);
  fprintf(stderr, " Input viscosity = %.10e\n ", mu1_in), fflush (stderr);
  fprintf(stderr, " MAX LEVEL = %d\n ", MAXLEVEL), fflush (stderr);
  fprintf(stderr, " MIN LEVEL = %d\n ", MINLEVEL), fflush (stderr);
  fprintf(stderr, " MAXTIME = %.10e\n ", MAXTIME), fflush (stderr);
  fprintf(stderr, " istart = %d\n ", istart), fflush (stderr);
  fprintf(stderr, "************************\n"), fflush (stderr);
  
  //Boxsize and initial grid
  size(2.*HALFBOX);
  init_grid(1<<MINLEVEL);

  //Physical properties
  rho1 = 1.0;//liquid is the reference
  rho2 = rho1/RHOR;

  mu1 = mu1_in;
  mu2 = mu1/MUR;

  //Filament length
  Lx = 9.;

  //Surface tension
  SIGMA = sq(mu1)/(rho1*sq(Oh)*Lx);
#if CLSVOF
  (const) scalar sigmaf[] = SIGMA;
  d.sigmaf = sigmaf;
#else
  f.sigma = SIGMA;
#endif

  //Run!
  run();

}

double profile (double x, double y, 
		double P1x, double P2x, double P1y) {

  double dist = nodata;
  if(     x <= P1x) {
    double a = (P2y-P1y)/(sq(P1x)-sq(P2x));
    double b = -2.*a*P1x;
    double c = P1y-a*sq(P1x)-b*P1x;
    dist = y - (a*sq(x)+b*x+c);
  }
  else if(x >= P1x && x <= P1x+P1y) {
    dist = sq(x-P1x)+sq(y-0.)-sq(P1y);
  }

  return dist;

}

# define POPEN(name, mode) fopen (name ".ppm", mode)

event init(i = 0) {
  
  /**
     A-posteriori check */

  fprintf(stderr, "************************\n"), fflush (stderr);
  fprintf(stderr, "Check of input parameters only\n"), fflush (stderr);
  fprintf(stderr, " Ohnesorge number = %.10e\n ", mu1/sqrt(rho1*SIGMA*Lx)), fflush (stderr);
  fprintf(stderr, " Thickness = %.10e\n ", thk), fflush (stderr);
  fprintf(stderr, " Input viscosity = %.10e\n ", mu1_in), fflush (stderr);
  fprintf(stderr, " MAX LEVEL = %d\n ", MAXLEVEL), fflush (stderr);
  fprintf(stderr, " MIN LEVEL = %d\n ", MINLEVEL), fflush (stderr);
  fprintf(stderr, " MAXTIME = %.10e\n ", MAXTIME), fflush (stderr);
  fprintf(stderr, " istart = %d\n ", istart), fflush (stderr);
  fprintf(stderr, " Surface tension = %.10e\n ", SIGMA), fflush (stderr);
  fprintf(stderr, "************************\n"), fflush (stderr);
  
  /**
     Set parameters for the parabolic profile */

  P1x = Lx, P1y = thk;
  P2x = 0., P2y = P1y+1.;
  
  /**
     Initial condition or restart */

  if(istart == 1) {
    refine( y <= 0.5*L0 && level < MAXLEVEL);
    
    // First, we define the profile
    vertex scalar ls_t[]; 
#if CLSVOF
    fprintf(stderr, "Start a new simulation with CLSVOF\n"), fflush (stderr);
    foreach() {
      d[] = profile(x, y, P1x, P2x, P1y);
    }
    foreach_vertex() {
      ls_t[] = (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.; 
    }
    fractions(ls_t, f);
#else
    fprintf(stderr, "Start a new simulation with VOF\n"), fflush (stderr);
    foreach_vertex() {
      ls_t[] = profile(x, y, P1x, P2x, P1y);
    }
    fraction(f, ls_t[]);
#endif

    // Finally, set the velocity to zero
    foreach() {
      foreach_dimension() {
        u.x[] = 0.;
      }
    }

    // Print an image for checking
    view(fov=20.,tx=-0.5);
    clear();
    draw_vof("f", filled = 1);
    box(notics=true);
    mirror({0, 1}) {
      draw_vof("f", filled = 1);
      box(notics=true);
    }
    {
      static FILE * fp = POPEN ("f_0", "a");
      save (fp = fp);
      fflush (fp);
    }
  
  }
  else {
    fprintf(stderr, "Restart from a prexisting field\n"), fflush (stderr);
    FILE * fp = fopen("restart.bin", "r");
    if( fp == NULL ) {
      fprintf(stderr, "Restarting file absent. Please provide it!\n"), fflush (stderr);
      return 1; 
    }
    else {
      restore (file = "restart.bin");
    }
  }
  
}

/** 
 We adapt on both the velocity field and the position of the interface. */

event adapt(i++) { //grid adaptation
  adapt_wavelet((scalar *){f, u}, (double[]){femax, uemax, uemax}, MAXLEVEL, MINLEVEL);
}

event movie(t += 0.1) {

  // Movie with the volume fraction
  view(fov=20.,tx=-0.5);
  clear();
  draw_vof("f", filled = 1);
  box(notics=true);
  mirror({0, 1}) {
    draw_vof("f", filled = 1);
    box(notics=true);
  }
  {
    static FILE * fp = POPEN ("f", "a");
    save (fp = fp);
    fflush (fp);
  }
 
}

/** 
   We output a dump file at the end of the simulation. */

event end(t = 0.5) {

  printf("%g", t);
  dump("end.bin");

} 