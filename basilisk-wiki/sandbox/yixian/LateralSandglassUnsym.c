/**
# Flow in a Sandglass / Hourglass with a lateral orifice in a silo with height H = 3.9L 
We propose an implementation of the Jop Pouliquen Forterre µ(I) rheology for the flow in a hourglass.
Like the case with an orifice at the bottom. We find that the flow through the lateral orifice alse follows the Beverloo-Hagen discharge law in a pure 2D case.
And finally we add the influence of a moderate friction of the front and back wall to simulate a 3D case
 
 # Code
 Includes and definitions
 */

#include "grid/quadtree.h"
#include "navier-stokes/centered.h"
#include "vof.h"
// Domain extent
#define LDOMAIN 4.
// heap definition
double  H0,R0,D,W,tmax,Q,Wmin,muwall,WDOMAIN,Hmask,Lz;

/**
   passive fluid small density to preserve 0 pressure
   and small viscocity
*/
#define RHOF 1e-4
#define mug  1e-5
// Maximum refinement
#define LEVEL 8
//#define MINLEVEL 7
char s[80];
FILE * fpf,*fwq;
scalar f[];
scalar * interfaces = {f}; 
face vector alphav[];
face vector muv[];
/**
   Boundary conditions for granular flow, pressure must be zero at the surface.
   The pressure is zero in the hole $x=1$ and $4\deltax < y< W+4\deltax$, but the lithostatic gradient is given elswhere 
   on the right wall.
   No slip boundary conditions on the other walls.
*/
p[top] = dirichlet(0);
u.n[top] = neumann(0);
u.t[bottom] =  dirichlet(0);
u.n[bottom] =  dirichlet(0);
u.n[left] = dirichlet(0);
u.t[left] = dirichlet(0);
f[left]= neumann(0);
u.t[right] = (fabs(y)>= 0.1 &&
	      fabs(y)<= (W+0.1)) ? neumann(0):  dirichlet(0);
u.n[right] = (fabs(y)>= 0.1&&
	      fabs(y)<= (W+0.1)) ? neumann(0):  dirichlet(0);
p[right]   = (fabs(y)>= 0.1 &&
	      fabs(y)<= (W+0.1)) ? dirichlet(0): neumann(0);


/**
   the three cases in the main */

int main(int argc, char ** argv) {
  L0 = LDOMAIN;
  // number of grid points
  N = 1 << LEVEL;
  // maximum timestep
  DT = 0.001;
  // coefficient of friction of wall
  muwall=0.1;
  TOLERANCE = 1e-3;
  //silo height
  H0 = 3.9;
  R0=20.000;
  // Grain size
  D=1./90.;
  // size of the hole
  W=0.5;
  fwq = fopen ("outWQ", "w");
  fclose(fwq);
  //film size
  Lz = LDOMAIN;
  //Silo thickness
  WDOMAIN = 0.5;
  const face vector g[] = {0.,-1.,0};
  a = g;
  alpha = alphav;
  mu = muv;
  Q = 0; 
  tmax = 1.;
  fpf = fopen ("interface.txt", "w");
  run();
  fclose (fpf);
  fprintf (stdout,"\n");
  fwq = fopen ("outWQ", "a");
  fprintf(fwq," %lf %lf \n", W, Q);
  fclose (fwq);
}


/**
   initial heap, a rectangle
*/
// note the v
event init (t = 0) {
  mask (x < 3.*L0/4. ? left : none);
  scalar phi[];
  foreach_vertex()
    phi[] = min(H0 - y, R0 - x);
  fractions (phi, f);
  /**
     lithostatic pressure, with a zero pressure near the hole
     to help 
  */
  foreach()
    p[] = (fabs(y-(W/2.+ 0.1 ))<= W/2. &&
	   fabs(x-LDOMAIN)<= .1) ?  0 : max(H0 - y,0) ;
}
/**
   total density 
*/

#define rho(f) ((f) + RHOF*(1. - (f)))
/**
   Viscosity computing $D_2=D_{ij}D_{ji}$; 

   In the pure shear flow 
   $D_{11}=D_{22}=0$ et $D_{12}=D_{21}=(1/2)(\partial u/\partial y)$,
   so that  
   $D_2=\sqrt{D_{ij}D_{ij}} =\sqrt{ 2 D_{12}^2} = \frac{\partial u}{ \sqrt{2}\partial y}$.
   In a pure shear flow, $\partial u/\partial y= \sqrt{2} D_2$.
   The inertial number $I$ is $D \sqrt{2} D_2/\sqrt(p)$ 
   and $\mu = \mu_s+ \frac{\Delta \mu}{1+I/I_0}$ 
   the viscosity is $\eta = \mu(I)p/D_2$:

   note that if $\eta$ is too small an artificial small viscosity $\rho D \sqrt{gD}$
   is taken see Lagrée et al. 11 § 2.3.1
*/
event properties (i++) {
  trash ({alphav});
  scalar eta[];
  foreach() {
    eta[] = mug;
    if (p[] > 0.) {
      double D2 = 0.;
      foreach_dimension() {
	double dxx = u.x[1,0] - u.x[-1,0];
	double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/2.;
	D2 += sq(dxx) + sq(dxy);
      }
      if (D2 > 0.) {
	D2 = sqrt(2.*D2)/(2.*Delta);
	double In = D2*D/sqrt(p[]);
	double muI = .4 + .28*In/(.4 + In);
	double etamin = sqrt(D*D*D);
	eta[] = max((muI*p[])/D2, etamin);
	eta[] = min(eta[],100);
      }
    }
  }
  boundary ({eta});
  scalar fa[];
  foreach()
    fa[] = (4.*f[] + 
	    2.*(f[-1,0] + f[1,0] + f[0,-1] + f[0,1]) +
	    f[1,1] + f[-1,1] + f[1,-1] + f[-1,-1])/16.;
  boundary ({fa});
  foreach_face() {
    double fm = (fa[] + fa[-1])/2.;
    muv.x[] = (fm*(eta[] + eta[-1])/2. + (1. - fm)*mug);
    // mu.x[] = 1./(2.*fm/(eta[] + eta[-1,0]) + (1. - fm)/mug);
    alphav.x[] = 1./rho(fm);
  }
  // boundary_normal ({mu,alpha});
  boundary ((scalar *){muv,alphav});
}
/**
   convergence outputs
*/
void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,  
	     exp (log (mg.resb/mg.resa)/mg.i));  
}
/**
   convergence stats
*/
event logfile (i += 1) {
  stats s = statsf (f);
  fprintf (stderr, "%g %d %g %g %g %g\n", t, i, dt, s.sum, s.min, s.max - 1.);
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
  fflush (stderr);
}


/**
   wall friction
   $$\frac{du}{dt} = \frac{-2 \mu_w p u}{W ||u||}$$
   equation with discretization:
   $$\frac{u^{n+1}-u^n}{\delta t} = \frac{-2 \mu_w p u^{n+1}}{W ||u^n||} $$
*/

#if dimension == 2
event friction (i++) {	
  foreach() { 
    double a = norm(u)==0 ? HUGE : 1. + 2.*muwall*dt*p[]/(norm(u)*WDOMAIN);
    foreach_dimension()
      u.x[] /= a;
  }
}
#endif



/**
   save some interfaces
*/
event interface (t = 0 ; t+=1. ; t<= tmax) {
#if dimension == 2
  output_facets (f, fpf);
#endif
  char s[80];
  sprintf (s, "field-%g.txt", t);
  FILE * fp = fopen (s, "w");
  output_field ({f,p,u,uf,pf}, fp, linear = true);
  fclose (fp);
}



/**
   Rate of flowing materials across the hole
*/
event debit (t += 0.1 ) {
  static double V=1; 
  V=0;
  foreach()
    V = V + f[]* Delta * Delta; 
  if(t>=0.) fprintf (stdout,"%lf %lf \n",t,V);
  fflush (stdout);
}  
/**
   film output
*/
#if 1
event movie (t += 0.05) {
  static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, fp1, min = 0, max = LEVEL, 
	      n = 256, box = {{0,0},{Lz,Lz}});

  foreach()
    l[] = f[]*(1+sqrt(sq(u.x[]) + sq(u.y[])));
  boundary ({l});
  static FILE * fp2 = popen ("ppm2mpeg > velo.mpg", "w");
  output_ppm (l, fp2, min = 0, max = 2., linear = true, 
	      n = 256, box = {{0,0},{Lz,Lz}});

  static FILE * fp3 = popen ("ppm2mpeg > f.mpg", "w");
  foreach()
    l[] = f[]*p[];
  output_ppm (l, fp3, min = 0, linear = true,
	      n = 256, box = {{0,0},{Lz,Lz}});
}
event pictures (t==3) {
  output_ppm (f, file = "f.png", min=0, max = 2,  spread = 2, n = 256, linear = true,
	      box = {{0,0},{2,2}});
}
#endif


/**
# Run

to run

~~~bash
qcc -g -O2 -Wall -o LateralUnsym LateralSandglassUnsym.c -lm
./LateralUnsym > out
~~~
 
Plot of pressure at time 1
~~~gnuplot pressure
reset
set pm3d; set palette rgbformulae 22,13,-31;unset surface;
set ticslevel 0;
unset border;
unset xtics;
unset ytics;
unset ztics;
unset colorbox;
#set xrange[-3:3];set yrange[-3:3];
set view 0,0
sp'./field-1.txt' u 1:2:($1>3. ? $4*$3 :NaN) t 'p'
~~~
 
 
Plot of velocity at time 1
 
~~~gnuplot velocity
reset
set pm3d; set palette rgbformulae 22,13,-31;unset surface;
set ticslevel 0;
unset border;
unset xtics;
unset ytics;
unset ztics;
unset colorbox;
#set xrange[-3:3];set yrange[-3:3];
set view 0,0
sp'field-1.txt' u 1:2:($1>3. ? sqrt($7*$7+$6*$6)*$3 :NaN)  t '|v|'
~~~
*/


/**
#Bibliography
L. Staron, P.-Y. Lagrée, & S. Popinet (2014)
"Continuum simulation of the discharge of the granular silo, A validation test for the μ(I) visco-plastic flow law"
Eur. Phys. J. E (2014) 37: 5 DOI 10.1140/epje/i2014-14005-6
 
L. Staron, P.-Y. Lagrée & S. Popinet (2012)
"The granular silo as a continuum plastic flow: the hour-glass vs the clepsydra"
Phys. Fluids 24, 103301 (2012); doi: 10.1063/1.4757390
*/
