/**
In this test case, we provide a simple Fluid Particle implementation
for the simulation of a single cylinder sitting inside a plane couette flow.
## Fluid Particle Method
[Tanaka Araki PRL 2000](https://doi.org/10.1103/PhysRevLett.85.1338)
### Embed + variable viscosity!
function is inspired from ghigo/src/myviscosity.h
and consist of a match of src/viscosity.h and src/viscosity-embed.h

## Free cylinder inside a plane couette
Dvinsky Popel 1987, Stokes Flow
[Dvinsky Popel 1987](https://doi.org/10.1016/0045-7930(87)90031-4).

*/
double RADIUS = 0.25; // cylinder radius
double xShift = 0.0;
double yShift = 0.0;
double HEIGHT = 1.0+1e-2; // channel height
double uFlow = 1.0;
double nu     = 1.; // fluid viscosity;
double nuRatio= 100.0;// particle-fluid viscosity ratio
double deltaT = 0.01;

//#include "embed.h"
#include "ghigo/src/myembed.h"
//#include "navier-stokes/centered.h"
#include "fpicella/src/centered_embed_variable_viscosity.h"
#include "fpicella/src/periodic-shift-treatment.h"


#include "view.h"

#define FLOW (2*y*uFlow)

#define channel difference(-y+HEIGHT/2.,-y-HEIGHT/2.)

#define SOLID solid(cs,fs,channel)

#define circle (sq (PS(x,xShift)) + sq (PS(y,yShift)) - sq (RADIUS))

#define rCenter sqrt(sq (PS(x,xShift)) + sq (PS(y,yShift)))

scalar csFP[];
face vector fsFP[];

scalar omega[];

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We define the mesh adaptation parameters. */

int lmin = 4;  // Min mesh refinement level (l=6 is 2pt/d)
int lmax = 6; // Max mesh refinement level (l=10 is 36pt/d, l=13 is 292pt/d)
#define cmax (1.e-2) // Absolute refinement criteria for the velocity field

FILE * output_file; // global file pointer

int main()
{
	display_control(xShift,-10,10);
	display_control(yShift,-1,1);
	display_control(nuRatio,0.,1000);
	display_control(RADIUS,0.1,0.45);
	display_control(lmin,4,10);
	display_control(lmax,4,10);
	display_control(deltaT,0.001,1000);
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (4.0 [0]);
  DT = 1.0 [0];

  origin (-L0/2., -L0/2.);

  stokes = true;

  periodic(left);
  //periodic(top);

  N = 1 << (lmax);//initialize at maximum refinement, so to have _child_

  init_grid (N);

  output_file = fopen ("Free_Particle.dat", "w");

	for(lmax = 8; lmax>=7; lmax -=1){
		for(yShift = 0.0; yShift<= 0.2; yShift += 0.05)
		{
			RADIUS = 0.25;
  		run();
		}

		for(RADIUS = 0.1; RADIUS <= 0.4; RADIUS += 0.05)
		{
			yShift = 0.;
  		run();
		}
	}
//	run();	

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
//  u.n[left]  = dirichlet(FLOW);
//  u.t[left]  = dirichlet(0.);
//    p[left]  = neumann(0.);
//  u.n[right] = dirichlet(FLOW);
//  u.t[right] = dirichlet(0.);
//    p[right] = neumann(0.);
	u.t[top]     = dirichlet(FLOW);
	u.n[top]     = dirichlet(0.);
	u.t[bottom]  = dirichlet(FLOW);
	u.n[bottom]  = dirichlet(0.);

  foreach_face()
    muv.x[] = (nu)*fs.x[];
  boundary ((scalar *) {muv});

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
 
  /**
	Compute fraction associated to a cylinder */
	fraction(csFP,circle);
/**
  use of PWP's macros from periodic-shift-treatment
  to aboid clumsy periodicity definitions...
*/
//// Compute viscosity...
	foreach_face()
		fsFP.x[] = face_value(csFP,0); // "naive" face interpolation
	// naive interpolation works better (for this purpose)
  // than solid(csFP,fsFP,circle).
  // WHY? To verify.

  foreach_face()
    muv.x[] = (nu)*fs.x[]+(1.-fsFP.x[])*nuRatio*nu - fs.x[]*(1.-fsFP.x[])*nu;
  /**
  We set the viscosity field in the event *properties*. */
  mu = muv;
  boundary ((scalar *) {muv});
}


event adapt (i+=10) // and not classic adapt, so to keep mesh identical during subiterations...
{
  SOLID;
  adapt_wavelet ({cs,csFP,u}, (double[]) {1.e-2,1e-2,(cmax),(cmax)},
                 maxlevel = (lmax), minlevel = (lmin));
}


scalar un[]; // field that will contain previous solution...
event logfile (t += deltaT; i <= 1000) {

	double Ux=0., Uy=0., OMEGA=0.;
	double counter=0;
  vorticity (u, omega);
	foreach(){
		if(csFP[]<1.0){
			Ux     += u.x[]  *sq(Delta)*(1.-csFP[]);
			Uy     += u.y[]  *sq(Delta)*(1.-csFP[]);
			OMEGA  += omega[]*sq(Delta)*(1.-csFP[]);
			counter+=1.0     *sq(Delta)*(1.-csFP[]);
		}
	}
	Ux    /= counter;
	Uy    /= counter;
	OMEGA /= counter*2.;
	double target_area = M_PI*sq(RADIUS);
	fprintf(stderr,"TOTO %+6.5e %+6.5e %+6.5e %+6.5e \n",Ux,Uy,OMEGA,(counter-target_area)/target_area*100);
	
/**
We look for a stationary solution. */
  double du = change (u.y, un);
  if (i > 0 && du < 1e-4){ // since I'm looking for a steady solution, this must be quite small...
                           // making the simulation quite expensive
  	fprintf(output_file,"%+3.2e %02d %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e \n",
                L0, lmax, RADIUS, yShift, Ux, Uy, OMEGA);
  	fflush(output_file);
    return 1; /* stop */
  }
}

/**
### Fixed cylinder location (center), varying radius, 
~~~gnuplot Free cylinder in plane Couette, y=0
set grid
set xlabel "Radius"
set ylabel "Omega"
set title "Free cylinder in plane Couette, y=0"

# Plot data from file, only lines with column 6 equal to 0 (i.e. yShift = 0)
plot '../Dvinsky_Popel_1987_fig_11.dat' using ( $1):( $2) \
     with lines lw 5 lc rgb "red" title "Dvinsky Popel 1987, fig 11",\
     'Free_Particle.dat' using ( $4==0 && $2==7 ? $3 : 1/0 ):( $4==0  && $2==7 ? +$7 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "blue" title "lmax 7",\
     'Free_Particle.dat' using ( $4==0 && $2==8 ? $3 : 1/0 ):( $4==0  && $2==8 ? +$7 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "green" title "lmax 8"
~~~
Good match!

### Fixed cylinder radius, changing it's location. Translational velocity.

~~~gnuplot Free cylinder in plane Couette, R=0.25
set grid
set xlabel "y"
set ylabel "U_x"
set title "Free cylinder in plane Couette, R=0.25"

# Plot data from file, only lines with column 6 equal to 0 (i.e. yShift = 0)
plot '../Dvinsky_Popel_1987_fig_10a_R_025.dat' using ($1-0.5):( $2) \
     with lines lw 5 lc rgb "red" title "Dvinsky Popel 1987, fig 10a",\
     'Free_Particle.dat' using ( $3==0.25 && $2==7 ? $4 : 1/0 ):( $3==0.25  && $2==7 ? +$5 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "blue" title "lmax 7",\
     'Free_Particle.dat' using ( $3==0.25 && $2==8 ? $4 : 1/0 ):( $3==0.25  && $2==8 ? +$5 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "green" title "lmax 8"
~~~
Overestimating translational velocity...and it does not seem to be dependent on mesh refinement nor dt...

### Fixed cylinder radius, changing it's location. Angular Velocity.

~~~gnuplot Free cylinder in plane Couette, R=0.25
set grid
set xlabel "y"
set ylabel "Omega"
set title "Free cylinder in plane Couette, R=0.25"

# Plot data from file, only lines with column 6 equal to 0 (i.e. yShift = 0)
plot '../Dvinsky_Popel_1987_fig_10b_R_025.dat' using ($1-0.5):( $2) \
     with lines lw 5 lc rgb "red" title "Dvinsky Popel 1987, fig 10b",\
     'Free_Particle.dat' using ( $3==0.25 && $2==7 ? $4 : 1/0 ):( $3==0.25  && $2==7 ? +$7 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "blue" title "lmax 7",\
     'Free_Particle.dat' using ( $3==0.25 && $2==8 ? $4 : 1/0 ):( $3==0.25  && $2==8 ? +$7 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "green" title "lmax 8"
~~~
Not bad at all for all cases...

...only for _low resolution_ (i.e. lmax =7) the gap between channel and cylinder is discretized with only 2 cells...

...then that's ok if accuracy is lost!
*/
