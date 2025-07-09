/**
### Asymmetric Couette flow
# Variable viscosity + embed
*/
#define EPS 1e-2
double HEIGHT = 1.0+EPS; // channel height
double uFlow = 1.0;
double nu     = 1.; // fluid viscosity;
double deltaT = 1.0;
double nu_alpha = 1.;
#define VISCOSITY_LAW exp(nu_alpha*y)
#define ANALYTIC_PROFILE (1.-exp(-nu_alpha*y))/(1.-exp(-nu_alpha*HEIGHT))

//#include "embed.h"
#include "ghigo/src/myembed.h"
//#include "navier-stokes/centered.h"
#include "fpicella/src/centered_embed_variable_viscosity.h"

#include "view.h"

#define FLOW (1*y*uFlow)

#define channel difference(-y+HEIGHT,-y+EPS)

#define SOLID solid(cs,fs,channel)

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We define the mesh adaptation parameters. */

int lmin = 4;  // Min mesh refinement level (l=6 is 2pt/d)
int lmax = 7; // Max mesh refinement level (l=10 is 36pt/d, l=13 is 292pt/d)
#define cmax (1.e-6) // Absolute refinement criteria for the velocity field

FILE * output_file; // global file pointer

int main()
{
	display_control(lmin,4,10);
	display_control(lmax,4,10);
	display_control(deltaT,0.001,1000);
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (4.0 [0]);
  DT = 1.0 [0];

  origin (-L0/2.,0.); // you could add whathever value here to shift the domain in y
											// you should get the same results
											// tweaking this value, you could have one (top) or two
											// (top and bottom) embedded boundaries.

  stokes = true;

  periodic(left);
  //periodic(top);

  N = 1 << (lmax);//initialize at maximum refinement, so to have _child_

  init_grid (N);

  output_file = fopen ("profiles.dat", "w");

	run();	

	fclose(output_file);
}

event init(i=0)
{
#if TREE
  /**
  When using *TREE*, we refine the mesh around the embedded
  boundary. */

  astats ss;
  int ic = 0;
  do {
    ic++;
		SOLID;
    ss = adapt_wavelet ({cs}, (double[]) {1.e-30},
                        maxlevel = (lmax), minlevel = (lmin));
  } while ((ss.nf || ss.nc) && ic < 100);
#endif // TREE
  /**
  Set boundary conditions */
	u.t[top]     = dirichlet(FLOW);
	u.n[top]     = dirichlet(0.);
	u.t[bottom]  = dirichlet(FLOW);
	u.n[bottom]  = dirichlet(0.);

  for (scalar s in {u, p, pf})
    s.third = true;
}



event properties (i++) // refresh particle's BC on embed at each iteration...!
{
  /**
  Set boundary conditions, on embed, must be updated with mesh refinement... */
  u.n[embed]  = dirichlet(FLOW);
  u.t[embed]  = dirichlet(0.0 );
  uf.n[embed] = dirichlet(FLOW);
  uf.t[embed] = dirichlet(0.0 ); 
 
//// Compute viscosity...
  foreach_face()
    muv.x[] = (nu)*fs.x[]*VISCOSITY_LAW;
  /**
  We set the viscosity field in the event *properties*. */
  mu = muv;
  boundary ((scalar *) {muv});
}


event adapt (i+=1) // and not classic adapt, so to keep mesh identical during subiterations...
{
  SOLID;
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(cmax),(cmax)},
                 maxlevel = (lmax), minlevel = (lmin));
}


scalar un[]; // field that will contain previous solution...
event logfile (t += deltaT; i <= 1000) {
/**
We look for a stationary solution. */
  double du = change (u.y, un);
  if (i > 10 && du < 1e-5){ // since I'm looking for a steady solution, this must be quite small...
  	for (double y = 0.; y <= HEIGHT; y += HEIGHT/100.)
    	fprintf (output_file, "%+6.5e %+6.5e %+6.5e %+6.5e\n"
					, y, interpolate (u.x, 0., y), 
					ANALYTIC_PROFILE,interpolate(muv.x,0,y));
  	fflush(output_file);
    return 1; /* stop */
  }
}

/**
### Velocity profile in a Couette flow, variable viscosity + embed,
~~~gnuplot Couette flow, eta = exp(alpha y)
set grid
set xlabel "u(y)"
set ylabel "y"
set title "Couette flow, eta = exp(nu_{alpha} y))"

# Plot data from file, only lines with column 6 equal to 0 (i.e. yShift = 0)
plot 'profiles.dat' using ($1):($2) with points pt 7 ps 2 lc rgb "green" title "numerical",\
     'profiles.dat' using ($1):($3) with line lw 3 ps 2 lc rgb "blue" title "analytical",
~~~
*/
