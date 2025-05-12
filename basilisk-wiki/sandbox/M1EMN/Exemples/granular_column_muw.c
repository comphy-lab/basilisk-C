/**
# Collapse of granular columns with basal friction
We propose an implementation of the Jop Pouliquen Forterre µ(I) rheology with friction at the base. This is the counterpart  of the [no slip](http://basilisk.fr/sandbox/M1EMN/Exemples/granular_column.c) case.
 
# Code 
Includes and definitions
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vof.h" 
// Domain extent
#define LDOMAIN 5.
// heap definition
double  H0,R0,D;
//film size
double  Lz;
//passive fluid
#define RHOF 1e-3
#define mug  1e-4
#define SIGMA 0.01
// Maximum refinement
#define LEVEL 8
char s[120];
FILE * fpf;
scalar f[];
scalar * interfaces = {f}; 
face vector alphav[];
face vector muv[];
scalar rhov[];
scalar eta[];
scalar etadudy[];
scalar muIp[];
vector u1[];
double tauw;
// note the v


// Boundary conditions
p[top] = dirichlet(-RHOF*LDOMAIN);
pf[top] = dirichlet(-RHOF*LDOMAIN);
u.n[top] = neumann(0);
u1.t[bottom] = neumann(0);
/**
We may impose no slip 
*/
//u.t[bottom] = dirichlet(0);
/** 
but here we impose a friction at the wall
$$\tau = \mu_w p$$
so that we impose a Neumann BC: 
$$\frac{\partial u}{\partial y}|_0 = \frac{\mu_w p}{\eta}$$
where $\eta$ is computed by the $\mu(I)$ law in the bulk. Note the minus sign (the normal os $-y$)), and note that we do some regularisation as we impose the stress in the reverse direction of the velocity
which may be zero (two possible regularizations). The solution is dependant on the regularisation parameter which should be not too small 'eps'.
*/
double eps=3e-3,mu_w=0.48;
//u.t[bottom] = neumann(mu.y[] ? -(mu_c*p[]/mu.y[])*u1.x[]/sqrt( sq(u1.x[]) + sq(eps)) : 0.);
u.t[bottom] = neumann(mu.y[] ? -(mu_w*p[]/mu.y[] * u1.x[]/(fabs(u1.x[] + eps)))  : 0.);
/**
see [couette with stress](http://basilisk.fr/sandbox/M1EMN/Exemples/couette_muw.c) case.*/

int main() {
  L0 = LDOMAIN;
  // number of grid points
  N = 1 << LEVEL;
  // maximum timestep
  DT = 2e-2;
  TOLERANCE = 1e-3;
// Initial conditions for the heap
  H0=1;
  R0=1.4286;
  // Grain size
  D=1./200;
  const face vector g[] = {0.,-1.};
  a = g;
  alpha = alphav;
  mu = muv;
  rho = rhov;
  Lz = LDOMAIN;
  run(); 
}

event init (t = 0) {
  fraction (f, min(H0 - y, R0 - x));
/*
initialisation of hydrostatic pressure for granular phase  
*/
  foreach()
    p[] = f[]*(H0-y);
}

#define rho(f) ((f) + RHOF*(1. - (f)))

event properties (i++) {
  trash ({alphav}); 
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
	D2 = sqrt(2.*D2)/(2.*Delta); // this D2 is sqrt(2) D2
	double In = D2*D/sqrt(p[]);
	double muI = .48 + .25*In/(.279 + In);
	double etamin = sqrt(D*D*D);
	eta[] = max((muI*p[])/D2, etamin);
	eta[] = min(eta[],100);      
	muIp[] = muI*p[]; }
    }
  }
  boundary ({eta,muIp});
  scalar fa[];
  foreach()
    fa[] = (4.*f[] + 
	    2.*(f[-1,0] + f[1,0] + f[0,-1] + f[0,1]) +
	    f[1,1] + f[-1,1] + f[1,-1] + f[-1,-1])/16.;
  boundary ({fa});
  foreach_face() {
    double fm = (fa[] + fa[-1,0])/2.;
    muv.x[] = (fm*(eta[] + eta[-1,0])/2. + (1. - fm)*mug);
    // mu.x[] = 1./(2.*fm/(eta[] + eta[-1,0]) + (1. - fm)/mug);
    alphav.x[] = 1./rho(fm);
  }
  foreach()  
  rhov[] = rho(fa[]);
 // boundary_normal ({mu,alpha});
 boundary ({muv,alphav,rhov});
}
/**
save stress and save old field
*/
event sylvain (i++) { 
  foreach()
    etadudy[] = muv.y[]*(u.x[0,0] - u.x[0,-1])/(Delta);
  boundary({etadudy});

  foreach()  
     u1.x[] = u.x[];
 boundary ({u1.x});
}

void mg_print (mgstats mg)
{
   if (mg.i > 0 && mg.resa > 0.)
        fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
                 mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0);
}

event logfile (i++) {
  stats s = statsf (f);
  fprintf (stderr, "%g %d %g %g %g %g\n", t, i, dt, s.sum, s.min, s.max - 1.);
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
  fflush (stderr);
}

event interface (t = {0,  0.8371, 1.5068, 1.6742, 2.5113 }) {
  //output_facets (f, fpf); 
  char s[120];
  sprintf (s, "field-%g.txt", t);
  FILE * fp = fopen (s, "w");
  output_field ({f,p,u,uf,pf,etadudy,muIp}, fp, linear = true);
  fclose (fp);
}

#if 0
event movie (t += 0.05) {
  static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, fp1, min = 0, max = LEVEL, 
	      n = 2048, box = {{0,0},{Lz,Lz}});

  foreach()
    l[] = f[]*(1+sqrt(sq(u.x[]) + sq(u.y[])));
  boundary ({l});
  static FILE * fp2 = popen ("ppm2mpeg > velo.mpg", "w");
  output_ppm (l, fp2, min = 0, max = 1.5, linear = true, 
	      n = 1024, box = {{0,-1},{3,1.5}});

  static FILE * fp3 = popen ("ppm2mpeg > f.mpg", "w");
  foreach()
    l[] = f[]*p[];
  output_ppm (l, fp3, min = 0, linear = true,
	      n = 2048, box = {{0,0},{Lz,Lz}});
}
 event pictures (t==3) {
  output_ppm (f, file = "f.png", min=0, max = 2,  spread = 2, n = 2048, linear = true,
  box = {{0,0},{2,2}});
}
#endif
/**

#Plots
 

Check that final $ p$ at the base multiplied by $\mu_s=0.48$ (chosen value) is the computed $\mu(I)p$

~~~gnuplot 
set xlabel "x"
  p[:3]'field-2.5113.txt'u 1:($3>.95? $2:NaN) t 'heap','' u 1:($2==0.00976562? 0.48*$4:NaN) t' mu p', ''u 1:($2==0.00976562? $11:NaN) t'mu(I)p'
~~~

 
Check that   $\mu(I)p$ is $\eta \partial u/\partial y$. this is not  when there is regularisation (which appears here when the flow is very slow)


~~~gnuplot 
set xlabel "x"
  p[:3]'field-2.5113.txt' u 1:($3>.95? $2:NaN)t 'heap',''u 1:($2==0.00976562? $10:NaN) t'eta du/dy' ,''u 1:($2==0.00976562? $11:NaN) t'mu(I)p'
~~~
 

*/
