/**
# Flow in a Sandglass / Hourglass with a bottom orifice

We propose an implementation of the Jop Pouliquen Forterre $\mu(I)$
rheology for the flow in a hourglass.  We find that the flow through
the orifice follows the Beverloo-Hagen discharge law in a pure 2D
case.  Here is an exemple of sandglass with an injection of air trough the orifice.  We
solve Poisson equation of the Darcy air flo in teh granular media and add the force of pressure to simulate the
injection of air

# Code

Includes and definitions
*/

#include "navier-stokes/centered.h"
#include "vof.h"
// Domain extent
#define LDOMAIN 4.
// heap definition
double  H0,R0,D,W,tmax,Q,Qair,Wmin,muwall,WDOMAIN,Hmask,maskr;
//film size
double  Lz;

/**
passive fluid small density to preserve 0 pressure
and small viscocity */

#define RHOF 1e-4
#define mug  1e-5
// Maximum refinement
#define LEVEL 8
// ratio of mask section
#define maskr 0.75
//#define MINLEVEL 7
char s[80];
FILE * fpf,*fwq;
scalar f[];
scalar * interfaces = {f};
face vector alphav[];
face vector muv[];
face vector av[];

//velocity of granular materials
double ugranular;

/**
Darcy fluid */

double betam = 100;
scalar pDarcy[],sourceDarcy[];
face vector beta[];
face vector gradpDarcy[];
mgstats mgpDarcy;

/**
Boundary conditions for air flow, we suppose that the regime is laminar
the flow rate of air is $Q_{air}$ at the top, pressure is zero at the ouput
from Darcy law:
$\nabla P = \frac{Q_{air}/S - U_p}{\beta} $
we set $\beta$ to a large value so that le pressure drop is very small (penalisation method)
*/

pDarcy[left]   = neumann(0);
pDarcy[right]  = neumann(0);
pDarcy[top]    = neumann((Qair / (LDOMAIN * (1 - maskr) ) - ugranular) / betam );
pDarcy[bottom] = fabs( x - (maskr + 1) * 0.5 * LDOMAIN ) <=  W / 2. ?
  dirichlet(0) : neumann(0);

/**
Boundary conditions for granular flow, pressure must be zero at the
surface.  The pressure is zero in the hole $|x|<=W/2$, but the
lithostatic gradient is given elsewhere on the bottom wall.  No slip
boundary conditions on the other walls. */

p[top]      = dirichlet(0);
u.n[top]    = neumann(0);
u.t[bottom] = fabs(x - (maskr + 1) * 0.5 * LDOMAIN) <= W / 2. ?
  neumann(0) : dirichlet(0);
u.n[bottom] = fabs(x - (maskr + 1) * 0.5 * LDOMAIN) <= W / 2. ?
  neumann(0) : dirichlet(0);
p[bottom]   = fabs(x - (maskr + 1) * 0.5 * LDOMAIN) <= W / 2. ?
  dirichlet(0) : neumann(0);
u.t[right]  = dirichlet(0);
u.t[left]   = dirichlet(0);

int main(int argc, char ** argv) {
  L0 = LDOMAIN;
  // number of grid points
  N = 1 << LEVEL;
  // maximum timestep
  DT = 0.001;
  // coefficient of friction of wall
  muwall = 0.1;
  TOLERANCE = 1e-3;
  H0 = 3.9;
  R0 = 20.000;
  // Grain size
  D = 1. / 90.;
  fwq = fopen("outWQ", "w"); // ????
  fclose(fwq);               // ????
  Lz = LDOMAIN;
  // size of the hole
  W = 0.125;
  Qair = 0.;
  WDOMAIN = 2.;
  a = av;
  alpha = alphav;
  mu = muv;
  Q = 0;
  tmax = 80.;
  fpf = fopen("interface.txt", "w");
  run();
  fclose(fpf);
  fprintf(stdout, "\n");
  fwq = fopen("outWQ", "a");
  fprintf(fwq, " %lf %lf \n", W, Q);
  fclose(fwq);
}

#if 0
int adapt() {
#if TREE
  astats s = adapt_wavelet ({f}, (double[]){5e-3}, LEVEL, MINLEVEL);
  return s.nf;
#else // Cartesian
  return 0;
#endif
}
#endif

/**
initial heap, a rectangle */

// note the v
event init (t = 0) {
  mask (x < maskr * L0 ? left : none);
  fraction (f, min(H0 - y, R0 - x));

  /**
  lithostatic pressure, with a zero pressure near the hole
  to help */

  foreach()
    p[] = (fabs(x-(maskr + 1)*0.5*LDOMAIN) <= W/2. && fabs(y) <= .1) ?
    0 : max(H0 - y,0);
  // the boundary conditions for the pressure need to be handled by
  // the Navier--Stokes solver
  
  /**
  initial for pressure darcy */

  foreach(){
    pDarcy[] = 0;
    sourceDarcy[] = 0.;
  }
  boundary ({pDarcy, sourceDarcy});
}

/**
total density */

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
is taken see Lagrée et al. 11 § 2.3.1. */

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
    alphav.x[] = 1./rho(fm);
  }
  boundary ((scalar *){muv,alphav});
}

/**
convergence outputs */

void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
             exp (log (mg.resb/mg.resa)/mg.i));
}

/**
convergence stats */

event logfile (i += 1) {
  stats s = statsf (f);
  fprintf (stderr, "%g %d %g %g %g %g\n", t, i, dt, s.sum, s.min, s.max - 1.);
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
  fflush (stderr);
}

event pressureDarcy (i++)
{
    
  /**
  'betam' is a penalisation parameter, so that pressure is constant ouside
  the grains and as
  $(\beta dp/dn)$ is constant across the interface
  and $\beta_m dp/dy = Qair$ at the top */

  foreach_face()
    beta.x[] = (f[] + betam*(1. - f[]));
  boundary ((scalar *){beta});
  mgpDarcy = poisson (pDarcy, sourceDarcy, beta);
}

event acceleration (i++)
{    
  double eps = 1.;

  /**
  Coupling gravity : Grad(p) added to gravity */
    
  foreach_face(y)
    av.y[] += - 1. - eps*(pDarcy[0,0] - pDarcy[0,-1])/Delta;
  foreach_face(x)
    av.x[] += - eps*(pDarcy[0,0] - pDarcy[-1,0])/Delta ;
  boundary ((scalar *){av});
}

/**
save some interfaces */

event interface (t = 0 ; t += 1. ; t <= tmax) {
#if dimension == 2
  output_facets (f, fpf);
#endif
  char s[80];
  sprintf (s, "field-%g.txt", t);

  foreach_face()
    gradpDarcy.x[] = (pDarcy[0] - pDarcy[-1])/Delta ;
  boundary ((scalar *){gradpDarcy});
  
  FILE * fp = fopen (s, "w");
  output_field ({f,p,u,uf,pf,pDarcy,gradpDarcy}, fp, linear = true);
  fclose (fp);
}

/**
Rate of flowing materials across the hole */

event debit (t = 0 ; t += 0.1 ; t <= tmax) {
  static double V = 1;
  V = 0;
  foreach()
    V = V + f[]* Delta * Delta;
  if (t >= 0.) fprintf (stdout,"%lf %lf %lf\n",t,V,ugranular);
  fflush (stdout);
}


/**
Calculate velocity of flowing materials in the bulk and add it to Qair at the top
to simulate the counter flow */

event velocitybulk (i++) {
  static double VVold, VV = 3.9,Qinst = 0;
  VVold = VV;
  VV = 0;
  foreach()
    VV = VV + f[]* Delta * Delta;
  Qinst = -(VV-VVold)/dt;
  ugranular = Qinst/(LDOMAIN*(1-maskr));
}

/**
film output */

#if 1
event movie (t += 0.05) {    
  static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
  scalar l[];
  foreach()
    l[] = level;
  boundary ({l});
  output_ppm (l, fp1, min = 0, max = LEVEL,
	      n = 256, box = {{0,0},{Lz,Lz}});
    
  foreach()
    l[] = f[]*(1 + sqrt(sq(u.x[]) + sq(u.y[])));
  boundary ({l});
    
  static FILE * fp2 = popen ("ppm2mpeg > velo.mpg", "w");
  output_ppm (l, fp2, min = 0, max = 2., linear = true,
	      n = 256, box = {{0,0},{Lz,Lz}});
    
  static FILE * fp3 = popen ("ppm2mpeg > pDarcy.mpg", "w");
  foreach()
    l[] = pDarcy[];
  boundary ({l});
  output_ppm (l, fp3, min = 0, linear = true,
	      n = 256, box = {{0,0},{Lz,Lz}});
}

event pictures (t==3) {
  output_ppm (f, file = "f.png",
	      min = 0, max = 2,  spread = 2, n = 256, linear = true,
	      box = {{0,0},{2,2}});
}
#endif


#if 0
event gfsview (i++)
{
  static FILE * fp = popen ("gfsview2D -s", "w");
  output_gfs (fp, t = t);
}
#endif
 
/**
# Run

to run

~~~bash
qcc -g -O2 -Wall -o Bottom BottomSandglassDarcy.c -lm
./Bottom > out
~~~

# Bibliography 

L. Staron, P.-Y. Lagrée, & S. Popinet (2014) "Continuum
simulation of the discharge of the granular silo, A validation test
for the μ(I) visco-plastic flow law" Eur. Phys. J. E (2014) 37: 5 DOI
10.1140/epje/i2014-14005-6.

L. Staron, P.-Y. Lagrée & S. Popinet (2012) "The granular silo as a
continuum plastic flow: the hour-glass vs the clepsydra" Phys. Fluids
24, 103301 (2012); doi: 10.1063/1.4757390.
*/
