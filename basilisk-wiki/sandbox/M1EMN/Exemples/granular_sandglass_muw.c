/**
# Flow in a Sandglass / Hourglass / Silo with a lateral orifice in a Hele-Shaw configuration
We propose an implementation of the Jop Pouliquen Forterre µ(I) rheology for the flow in a silo with a lateral orifice of hight $D_w$.
Like the case with an orifice at the bottom, we find that the flow follows the Beverloo-Hagen discharge law in a pure 2D case. 
For the 3D case, we use the ideas of Hele-Shaw cell, integrate across the width $W_{domain}$. 
This adds the influence of a moderate friction at the front and back wall to simulate a shallow 3D case by a 2D set of equations with an extra source term.

Equations are incompressibility:
$\nabla \cdot u =0$

and momentum: 
$\frac{du}{dt} = - \nabla p + \nabla \cdot (2 \eta D) + \rho g  - \frac{2 \mu_w p\; u}{W_{domain} |u|}$ with the viscosity $\eta = \frac{\mu(I)p}{\sqrt{2}D_2}$ ($D$, rate of strain tensor, $D_2$ his second invariant)

If we use $D_w$ the vertical size of the hole as length scale (do not confuse $D_w$, with $D_{grain}$ the grain diameter and $D$ the rate of strain tensor), the time scale is $\sqrt{D_w/g}$, the velocity scale is $\sqrt{g D_w}$, the  pressure scale is $\rho gD_w$.


without dimension momentum is:

$$\frac{du}{dt} = - \nabla p + \nabla \cdot (\sqrt{2} \mu(I) \frac{D}{D_2}) - e_z  - \frac{2 \mu_w D_w p\; u}{W_{domain} |u|}$$

Grain diameter is in the $\mu(I)$ function and 
$\frac{2 \mu_w D_w  }{W_{domain} }$ is the Hele-Shaw friction coefficient which contains the geometry and the wall friction.
This is the parameter that we will change.



# Code 
Includes and definitions
*/
#include "grid/multigrid.h"
//#include "grid/quadtree.h"
////#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "vof.h"
// Domain extent
#define LDOMAIN 4.
// heap definition
double  H0,R0,Dgrain,DW,tmax,Q,muwall,sintheta,WDOMAIN;
/**
passive fluid small density to preserve 0 pressure
and small viscocity
*/
#define RHOF 1e-4
#define mug  1e-5
// Maximum refinement
#define LEVEL 7
char s[80];
FILE * fpf,*fwq;
scalar f[];
scalar * interfaces = {f}; 
face vector alphav[];
face vector muv[];
scalar rhov[];
/**
Boundary conditions for granular flow, pressure must be zero at the surface.
The pressure is zero in the hole $x=1$ and $h_w < y< D_w+h_w$, but the lithostatic gradient is given elswhere 
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
double hw=0.09375;  // a given hight
u.t[right] = (fabs(y)>= hw && fabs(y)<= (DW+hw)) ? neumann(0):  dirichlet(0);
u.n[right] = (fabs(y)>= hw && fabs(y)<= (DW+hw)) ? neumann(0):  dirichlet(0);
p[right]   = (fabs(y)>= hw && fabs(y)<= (DW+hw)) ? dirichlet(0): neumann(0);
/** 
main */
int main(int argc, char ** argv) {

  disable_fpe (FE_DIVBYZERO|FE_INVALID);

  L0 = LDOMAIN;
  // number of grid points
  N = 1 << LEVEL;
  // coefficient of friction of wall
  muwall=.1;
  TOLERANCE = 1e-3; 
  // Initial conditions a=.5
  H0=3.9;
  R0=20.000;
  // Grain size
  Dgrain=1./90.;
  
  
  const face vector g[] = {0.,-1.,0};
  a = g;
  alpha = alphav;
  mu = muv;
  rho=rhov;

  
  fwq = fopen ("outWQW", "w");
  fclose(fwq); 
    
  int iWW[11]={ 12 , 14 , 18, 20, 28, 76, 140 , 268 , 1024, 2048, 5000 }; 
  
  for (int iW=0;iW<11;iW++){ 
   Q = 0; 
   // maximum timestep DT = 0.002 for L7  DT =0.001 for L 8
   DT = 0.004;    
   tmax = 4.;  
   DW=0.4375; // size of the hole  
   DW=1; 
   WDOMAIN = (iWW[iW])*Dgrain; // width of the cell in grain size
   run(); 
   fprintf (stdout,"\n");
   fwq = fopen ("outWQW", "a");
   fprintf(fwq," %lf %lf %lf %lf %d  \n", DW, Q , WDOMAIN, sintheta, iW);
   fclose (fwq);
   }
}
/**
initial heap, a rectangle
*/ 
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
     p[] = (fabs(y-(DW/2.+ 4.*LDOMAIN/pow(2,LEVEL) ))<= DW/2. && fabs(x-LDOMAIN)<= .1) ?  0 : max(H0 - y,0) ;   
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
the viscosity is $\eta = \frac{\mu(I)p}{\sqrt{2}D_2}$:

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
  D2 = sqrt(D2)/(2.*Delta);  // this is D2     
  double sD2 = sqrt(2.)*D2; // this sD2 is (sqrt(2) D2)
  double In = sD2*Dgrain/sqrt(p[]);
  double muI = .4 + .28*In/(.4 + In);
  double etamin = sqrt(Dgrain*Dgrain*Dgrain);
  eta[] = max((muI*p[])/sD2, etamin);
  eta[] = min(eta[],100);      }
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
    double fm = (fa[] + fa[-1,0])/2.;
    muv.x[] = (fm*(eta[] + eta[-1,0])/2. + (1. - fm)*mug);
    // mu.x[] = 1./(2.*fm/(eta[] + eta[-1,0]) + (1. - fm)/mug);
    alphav.x[] = 1./rho(fm);
  }
foreach()
    rhov[] = rho(fa[]); 
 boundary ({muv,alphav,rhov});
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
event logfile (i++) {
  stats s = statsf (f);
  fprintf (stderr, "%g %d %g %g %g %g\n", t, i, dt, s.sum, s.min, s.max - 1.);
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
  fflush (stderr);
}


/**
wall friction is introduced as a source term in the 2D equations (Hele-Shaw).
We split the influence of the wall friction:
$$\frac{du}{dt} = \frac{-2 \mu_w p\; u}{W_{domain} |u|}$$
If we write velocity as  $u = |u| T$, where the tangential vector is $T=\frac{u}{|u|}$. 
We introduce the [Frenet base](https://en.wikipedia.org/wiki/Jean_Frédéric_Frenet) with $T$ and $N$.


So 
$\frac{du}{dt} = \frac{-2 \mu_w p\; }{W_{domain} } T$, but $\frac{du}{dt} = \frac{d|u|}{dt} T + |u|^2  \frac{N}{R}$, 
 as $\frac{dT}{ds}=\frac{N}{R}$, hence there is no curvature, and the velocity remains in the same direction
$$\frac{d|u|}{dt} = \frac{-2 \mu_w p }{W_{domain} }$$
therefore the exact integration gives:
$$u_i^{n+1}= u_i^n  - \Delta t \frac{2 \mu_w p  u_i^{n}}{W_{domain}|u^n|} $$
as velocity decreases, it can not be less than 0, hence $|u|^{n+1}=max(|u|^{n+1},0)$.
*/
event friction (i++) {
    foreach()
    {
        double m = 2.*muwall*dt*p[]/WDOMAIN;
        double U = norm(u);
        foreach_dimension()
        u.x[]= U > 0 ? max(U-m,0)*u.x[]/U : 0 ; 
    }
}
/**
note that for $\mu_{wall}=0.1$, and $W_{domain}$ is  measured in grain diameter ($D=1/90$) for comparison with experiments,
we have 
$2 \mu_{wall}/W_{domain}= 18/n_{grains}$, this coefficient is one for 18 grains. This explains that the computation breaks for about 13 grains.
*/


/**
Rate of flowing materials across the hole
*/
event debit (t += 0.05; t <= tmax ) {
  static double Vold,V=1,Qinst=0;
  Vold=V; 
  V=0;
  foreach()
    V = V + f[]* Delta * Delta;
  Qinst = -(V-Vold)/0.05;  
  if(Qinst > Q) Q = Qinst;
   double ux=interpolate (u.x, LDOMAIN-0.01, hw+DW/2);
   double uy=interpolate (u.y, LDOMAIN-0.01, hw+DW/2);; 
   double U;
   U=sqrt(sq(ux)+sq(uy));
   sintheta = (U>0 ? fabs(ux/U) : 0);
  if(t>=.1) fprintf (stdout,"%lf %lf %lf %lf %lf \n",t,V/L0/H0,DW,Q,sintheta); 
  fflush (stdout);
}  
 
/**
# Run

to run

~~~bash
qcc  -g -O2 -o granular_sandglass_muw granular_sandglass_muw.c -lm 
./granular_sandglass_muw > out
~~~

Plot of the curve showing the influence of $D_w/W_{domain}$
 

~~~gnuplot extended Beverloo 
set logscale
set xlabel "D/W"
set ylabel "Q/W^{5/2}"
set key top left  
p [:10][:100] 'outWQW' u ($1/$3):($2/($3**(1.5))) t'calc'   w lp, .76*x**1.5*sqrt(1/(1+0.26*x)) t'fit article .76*x**1.5*sqrt(1/(1+0.26*x))', 1.49*x linec -1 , .76*x**1.5 linec -1 

~~~



*/ 


/**
# Related examples
 
* [The granular column](http://basilisk.fr/_showraw/sandbox/M1EMN/Exemples/granular_column.c)

* [Silo with aperture at the bottom](http://basilisk.fr/sandbox/M1EMN/Exemples/granular_sandglass.c)
 
#Bibliography

* [Y. Zhou](http://www.theses.fr/2016AIXM4731) Ejection de gaz et de grains suite à la rupture d'un crayon de combustible nucléaire : modélisation de la dynamique, thèse soutenue le 2 Novembre 2016

* L. Staron, [P.-Y. Lagrée](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/epje130141.pdf), & S. Popinet (2014)
"Continuum simulation of the discharge of the granular silo, A validation test for the μ(I) visco-plastic flow law" 
Eur. Phys. J. E (2014) 37: 5 DOI 10.1140/epje/i2014-14005-6

* L. Staron, [P.-Y. Lagrée](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/PhysFluids_24_103301.pdf) & S. Popinet (2012)
"The granular silo as a continuum plastic flow: the hour-glass vs the clepsydra" 
Phys. Fluids 24, 103301 (2012); doi: 10.1063/1.4757390  

* Y. Zhou, [P.-Y. Lagrée](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/ZhouJFM2017.pdf), S. Popinet, P. Ruyer, P. Aussillous (2017)
"Experiments on, and discrete and continuum simulations of, the discharge of granular media from silos with a lateral orifice", 
Journal of Fluid Mechanic   https://doi.org/10.1017/jfm.2017.543
*/
