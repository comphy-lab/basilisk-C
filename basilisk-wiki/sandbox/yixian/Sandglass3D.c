/**
# Flow in a Sandglass / Hourglass with an orifice at the bottom in 3D case
   We propose an implementation of the Jop Pouliquen Forterre µ(I) rheology for the flow in a hourglass.
   
   
## Code 
   Includes and definitions
*/
// #include "grid/multigrid.h"
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "vof.h"
// Domain extent
#define LDOMAIN 1.
// heap definition
double  H0,R0,D,W,tmax,Q,Wmin,DW;
//film size
double  Lz;
/**
   passive fluid small density to preserve 0 pressure
   and small viscocity
*/
#define RHOF 1e-4
#define mug  1e-5
// Maximum refinement
#define LEVEL 5 //7
char s[80];
FILE * fpf,*fwq;
scalar f[];
scalar * interfaces = {f}; 
face vector alphav[];
face vector muv[];
scalar rhov[];

/**
   Boundary conditions for granular flow, pressure must be zero at the surface.
   The pressure is zero in the hole $x=1$ and $4\deltax < y< W+4\deltax$, but the lithostatic gradient is given elswhere 
   on the right wall.
   No slip boundary conditions on the other walls.
*/

p[top] = dirichlet(0);
u.n[top] = neumann(0);
u.t[bottom] =  (x-LDOMAIN/2)*(x-LDOMAIN/2)+ (z-LDOMAIN/2)*(z-LDOMAIN/2) <= W/2.*W/2. ? neumann(0):  dirichlet(0);
u.n[bottom] =  (x-LDOMAIN/2)*(x-LDOMAIN/2)+ (z-LDOMAIN/2)*(z-LDOMAIN/2) <= W/2.*W/2. ? neumann(0):  dirichlet(0);
//u.n[bottom] = dirichlet(0);
p[bottom]   =  (x-LDOMAIN/2)*(x-LDOMAIN/2)+ (z-LDOMAIN/2)*(z-LDOMAIN/2) <= W/2.*W/2. ? dirichlet(0): neumann(0); 
u.n[right] = dirichlet(0);
u.n[left] = dirichlet(0);
u.t[right] = dirichlet(0);
u.t[left] = dirichlet(0);
u.n[front] = dirichlet(0);
u.t[front] = dirichlet(0);
u.n[back] = dirichlet(0);
u.t[back] = dirichlet(0);
f[left]= neumann(0);



/**
   the three cases in the main */

int main(int argc, char ** argv) {
  L0 = LDOMAIN;
  // number of grid points
  N = 1 << LEVEL;
  // maximum timestep
  DT = 0.001;
 
  TOLERANCE = 1e-3; 
  // Initial conditions a=.5
  H0=0.96;
  R0=20.000;
  // Grain size
  D=1./90.;
  // size of the hole
  W = 0.25;
  //Wmin = 0.15625;
  DW = LDOMAIN/N;
  fwq = fopen ("outWQ", "w");
  fclose(fwq);
  Lz = LDOMAIN;
  const face vector g[] = {0.,-1.,0};
  a = g;
  alpha = alphav;
  mu = muv;
  rho = rhov;

  
  Q = 0; 
  tmax = 3;
  fpf = fopen ("interface.txt", "w");
  run();
  fclose (fpf);
  fprintf (stdout,"\n");
  fwq = fopen ("outWQ", "a");
  fprintf(fwq," %lf %lf \n", W, Q);
  fclose (fwq);
  
}

/**
   initial heap, a paralepipede
*/
// note the v
event init (t = 0) {
 // mask (x < 3.*L0/4. ? left : none);
  
  scalar phi[];
  foreach_vertex()
    phi[] = min(H0 - y, R0 - x);
  fractions (phi, f);
  /**
     lithostatic pressure, with a zero pressure near the hole
     to help 
  */
  foreach()
    p[] = ((x-LDOMAIN/2)*(x-LDOMAIN/2)+ (z-LDOMAIN/2)*(z-LDOMAIN/2) <= W/2.*W/2. && fabs(y)<= .1) ?  0 : max(H0 - y,0) ;
  boundary ({p});
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

/**
 foreach_dimension() {
  double dxx = u.x[1,0] - u.x[-1,0];
  double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/2.;
  D2 += sq(dxx) + sq(dxy);
      }
*/        
  foreach() {
    eta[] = mug;
    if (p[] > 0.) {
      double D2 = 0.;
      double dxx = u.x[1,0,0] - u.x[-1,0,0];
      double dyy =  u.y[0,1,0] - u.y[0,-1,0];
      double dzz =  u.z[0,0,1] - u.z[0,0,-1];
      double dxy = (u.x[0,1,0] - u.x[0,-1,0] + u.y[1,0,0] - u.y[-1,0,0])/2.;
      double dxz = (u.x[0,0,1] - u.x[0,0,-1] + u.z[1,0,0] - u.z[-1,0,0])/2.;
      double dyz = (u.y[0,0,1] - u.y[0,0,-1] + u.z[0,1,0] - u.z[0,-1,0])/2.;
      D2 = sq(dxx) + sq(dyy) +sq(dzz) +2.*sq(dxy) +2.*sq(dxz) +2.*sq(dyz);
      if (D2 > 0.) {
	D2 = sqrt(2.*D2)/(2.*Delta);
	double In = D2*D/sqrt(p[]);
	double muI = .32 + .28*In/(.4 + In);
	double etamin = sqrt(D*D*D);
	eta[] = max((muI*p[])/D2, etamin);
	eta[] = min(eta[],100);
      }
    }
  }
  boundary ({eta});
  scalar fa[];
  foreach()
    fa[] = (8.*f[] + 
	    4.*(f[-1,0,0] + f[1,0,0] + f[0,-1,0] + f[0,1,0]+ f[0,0,1]+ f[0,0,-1]) +
	        2*(f[1,1,0] + f[-1,1,0] + f[1,-1,0] + f[-1,-1,0]+f[0,1,1] + f[0,-1,1] +
                   f[0,1,-1] + f[0,-1,-1]+f[1,0,1] + f[-1,0,1] + f[1,0,-1] + f[-1,0,-1])+ 
          f[1,1,1]+ f[1,1,-1]+f[-1,1,-1]+f[-1,1,1]+f[1,-1,1]+ f[1,-1,-1]+f[-1,-1,-1]+f[-1,-1,1])/64.;
  boundary ({fa});
  foreach_face() {
    double fm = (fa[] + fa[-1])/2.;
    muv.x[] = (fm*(eta[] + eta[-1])/2. + (1. - fm)*mug);
    // mu.x[] = 1./(2.*fm/(eta[] + eta[-1,0]) + (1. - fm)/mug);
    alphav.x[] = 1./rho(fm);
  }
   foreach()
    rhov[] = rho(fa[]); 
    boundary ({muv,alphav,rhov});

  //boundary ((scalar *){muv,alphav});
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
   Rate of flowing materials across the hole
*/
event debit (t = 0 ;t += 0.1 ; t<= tmax  ) {
  static double V=1; 
  V=0;
  foreach()
    V = V + f[]* Delta * Delta* Delta; 
  if(t>=0.) fprintf (stdout,"%lf %lf \n",t,V);
  fflush (stdout);
}  
/**
   film output for 2D
*/
#if 0
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
gfsview...
*/
#if 0
event gfsview (i += 1) {
  static FILE * fp = popen (dimension == 2 ?
			    "gfsview2D gfsview.gfv" :
			    "gfsview3D gfsview3D.gfv",
			    "w");
  output_gfs (fp, t = t);
}
#endif



/**
   if "grid/multigrid.h" is not included, then this is the quadtree with possibility of adaptation
*/
#if 0 // QUADTREE
event adapt (i++) {
adapt_wavelet ({f,u.x,u.y}, (double[]){5e-3,0.001,0.001}, LEVEL, LEVEL-2,
	       list = {p,u,pf,uf,g,f});
}
#endif
/**
   # Run

   to run

~~~bash
   qcc -g -O2 -Wall -o Sandglass3D Sandglass3D.c -lm
   ./Sandglass3D > out
~~~

## Result

We plot the volume as function of time. it is linear after a transition, this is Beverloo law:
~~~gnuplot V(t)
p'out' w l
~~~

*/ 


/**
#Bibliography

* L. Staron, P.-Y. Lagrée, & S. Popinet (2014)
   "Continuum simulation of the discharge of the granular silo, A validation test for the μ(I) visco-plastic flow law" 
   Eur. Phys. J. E (2014) 37: 5 DOI 10.1140/epje/i2014-14005-6

*   L. Staron, P.-Y. Lagrée & S. Popinet (2012)
   "The granular silo as a continuum plastic flow: the hour-glass vs the clepsydra" 
   Phys. Fluids 24, 103301 (2012); doi: 10.1063/1.4757390  
   
*  Y. Zhou PhD
*/
