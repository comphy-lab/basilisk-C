/** 
# A self-similar solver for the viscous collapse of a heap 

This is an adaptation of the 
[code developed by *P.-Y. Lagrée*](http://basilisk.fr/sandbox/M1EMN/Exemples/column_viscous.c) 
for exhibiting the self-similar behaviour of the collapse of 
a viscous heap for a newtonian fluid.

Here, we want to know if a [*self-similar solver*](http://basilisk.fr/sandbox/cailler/self_sim_DNS/README) 
can be adaptated to this problem, as self-similar coordinates'scalings differ in 
the two directions (see the problem transformed in self-similar coordinates 
in the below `selfsim_centered_huppert.h` file). 
*/

#include "selfsim_centered_huppert.h"
#include "contact.h"
#include "selfsim_two-phase_huppert.h"



// Domain extent
#define LDOMAIN 4.000001
//passive fluid
#define RHOF 1e-3
#define mug  1e-3
double mub;
double Bi=0;
double tmax;
double etamx;
// Maximum refinement 4./2**8 = 0.015625
// Minimum refinement 2./2**4 = 0.125
#define LEVEL 7 //7 // 9 OK
#define LEVELmin 2

scalar eta[];
scalar omega[]; // vorticity field
/**
## Boundary conditions
 */
 u.t[bottom] = dirichlet(0);

/**
## Main Parameters
*/
int main() {
  L0 = LDOMAIN;
  // number of grid points
  N = 1 << LEVEL;
  // maximum timestep
  DT = 0.025/4 ;
  TOLERANCE = 1e-3;
  
  /**
   **Numerical values**
    
  The problem is a heap of viscoplastic material of height = length
    released on a flat plate.
    We may define a cone for Abrams slup test (concrete).
    */
  // Initial conditions
#define Ltas 0.9999999999
#define Ls 0.9999999999
#define Htas 0.99999999999

/**
     
The fluid has a density $\rho_0$, time is $T=\sqrt{L/g}$,
We define a quantity, $\mu_b$ which is in fact the inverse of
$Re_1=\frac{\rho_0 U_0L}{\mu_N }$     and wich is here 1.
and we define the *Bingham* number $Bi=(\frac{\tau_y L}{\mu_N U_0})$, 
we change it.
 
*/
    
  mub= 1;
  etamx=1000.;
  //tmax = 199.99;
  tmax = 100;
  DT = 0.025;

  rho1 = 1., rho2 = 1.e-3 ;
  mu1 = 1., mu2 = 1e-3; 

  run();   
}

/**
## Main Events
*/
event init (t = 0) {
// Defining an initial trapezoidal shape:
  fraction (f, min(Htas - y, Ltas - (fabs(x)+ y * (Ltas-Ls))));
}

event vorti (i++){
  vorticity (u, omega);
}

#if QUADTREE
// if no #include "grid/multigrid.h" then adapt
event adapt (i++) {
  scalar K1[],K2[];
  foreach()
    K1[]=(f[0,1]+f[0,-1]+f[1,0]+f[-1,0])/4*noise();
  boundary({K1});

  for(int k=1;k<8;k++)
  {
    foreach()
    K2[]=(K1[0,1]+K1[0,-1]+K1[1,0]+K1[-1,0])/4;
    boundary({K2});
    foreach()
    K1[]=K2[]*noise();
  }
  adapt_wavelet({K1,f},(double[]){0.001,0.01}, maxlevel = LEVEL, minlevel = LEVELmin);
  //adapt_wavelet ({f,u,muv.x}, (double[]){5e-3,0.02,0.02,0.01}, LEVEL, LEVELmin,list = {p,u,pf,uf,g,f});
}
#endif

event stop (t = tmax);


event monitoring (i++){
  printf("time = %g ; iteration #%i\n", t, i);
  fflush(stdout); // empty buffer
}


/**
## Outputs

**Convergence outputs:**
*/
void mg_print (mgstats mg)
{
#if 0
    if (mg.i > 0 && mg.resa > 0.)
        fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
                 exp (log (mg.resb/mg.resa)/mg.i));
#else
  if (mg.i > 0 && mg.resa > 0.)
    fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
              mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0); 
#endif
}

event logfile (i++) {
  stats s = statsf (f);
  fprintf (stderr, "%g %d %g %g %g %g\n", t, i, dt, s.sum, s.min, s.max - 1.);
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
  fflush (stderr);
}

/** 
 **Some interface profiles for comparison:**
*/
event pij0 (t = 0.) {
  char s[80];
  sprintf (s, "out0-%g",Bi);
  FILE * fp = fopen (s, "w");
  output_facets (f, fp);
  fclose (fp);
}
event pij2 (t = 2.5001234) {
  char s[80];
  sprintf (s, "out1-%g",Bi);
  FILE * fp = fopen (s, "w");
  output_facets (f, fp);
  fclose (fp);
}
event pij5 (t = 4.9876) {
  char s[80];
  sprintf (s, "out2-%g",Bi);
  FILE * fp = fopen (s, "w");
    output_facets (f, fp);
  fclose (fp);
}
event pij10 (t = 9.9867212){
  char s[80];
  sprintf (s, "out3-%g",Bi);
  FILE * fp = fopen (s, "w");
  output_facets (f, fp);
  fclose (fp);
}
event pij40 (t = 39.9867212){
  char s[80];
  sprintf (s, "out4-%g",Bi);
  FILE * fp = fopen (s, "w");
  output_facets (f, fp);
  fclose (fp);
}
event pij50 (t = 49.9867212){
  char s[80];
  sprintf (s, "out5-%g",Bi);
  FILE * fp = fopen (s, "w");
  output_facets (f, fp);
  fclose (fp);
}
event pijtmax (t = tmax){
  char s[80];
  sprintf (s, "out6-%g",Bi);
  FILE * fp = fopen (s, "w");
  output_facets (f, fp);
  fclose (fp); 
}


/**
 **Saving the top position and runout position for slump measurements:**
*/
vector h[];
event timeseries (t += 0.1 ) {
  heights (f, h);
  double maxy = - HUGE,maxx = - HUGE;;
  foreach()
  if ((h.y[] != nodata) && (h.x[] != nodata)) {
    double yi = y + height(h.y[])*Delta;
    double xi = x + height(h.x[])*Delta;
    if (yi > maxy)
      maxy = yi;
    if (xi > maxx)
      maxx = xi;
  }
  char s[80];
  sprintf (s, "hmax-%g",Bi);
  static FILE * fp0 = fopen (s, "w");
  fprintf (fp0, "%g %g %g\n", t, maxx, maxy);
  fflush (fp0);
}


// #ifdef gnuX
// event output (t += .1 ) {
//     fprintf (stdout, "p[0:4][0:2]  '-' u 1:2    not  w   l \n");
//     output_facets (f, stdout);
//     fprintf (stdout, "e\n\n");
// }
// #endif




/**

## Links
 
 * 1D viscous collapse
 * Multilayer viscous collapse
 * 1D Bingham collapse
 * Multilayer Bingham collapse
 * granular collapse
 * http://basilisk.fr/sandbox/M1EMN/Exemples/column_SCC.c
 
 
## Bibliography
 

* [Huppert H.](http://www.itg.cam.ac.uk/people/heh/Paper47.pdf)
 ”The propagation of two-dimensional and axisymmetric viscous gravity currents over a rigid horizontal surface” J . Fluid Mech. (1982), vol. 121, p p . 43-58
 * Lagrée  [M1EMN
Master 1 Ecoulements en Milieu Naturel](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
 * Lagrée  [M2EMN
Master 2 Ecoulements en Milieu Naturel](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/mainM2EMN.pdf)


 
*/