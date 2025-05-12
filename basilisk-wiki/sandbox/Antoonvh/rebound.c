/**
#Perservere

If you do not succeed at first, try again. We target a dipolar vortex
towards an opening that is not big enough for the vortex to pass
through. It turns out this is due to the vorticity that is generted by
"back flow" as the vortex 'tries' to pass. The opening geometry is
implemented via the embedded boundary formulation. The set-up is
inspired by the results presented in J.A. van Hooft's and C.H. Wong's
Master Thesises (see references).
*/
#include "embed.h" 
#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"
#include "view.h"

#define RAD (pow(pow((x - xo), 2)+(pow((y - yo), 2)), 0.5))
#define ST (-(x - xo)/RAD)

double xse = 1., yse = 0.1, xo = 0, yo = -5;
int maxlevel = 11;

/**
   A function to initialize the dipolar vortex is given below.
 */
void init_lamb (double sig, double xo, double yo){
  double k = 3.83170597;
  refine (sq(x - xo) + sq(y - yo) < 2 && level< maxlevel);
  scalar psi[], omg[];
  foreach() 
    omg[] = -sig*((RAD<1)*((-2*j1(k*RAD)*ST/(k*j0(k)))))*sq(k);
  boundary ({omg});
  poisson (psi, omg);
  boundary ({psi});
  foreach(){
    u.x[] += ((psi[0,1] - psi[0,-1])/(2*Delta));
    u.y[] += -(psi[1,0] - psi[-1,0])/(2*Delta);
  }
}

face vector muc[];
double Re = 3000;

u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.); 

int main(){ 
  L0 = 7*pi;
  X0 = Y0 = -L0/2;
  mu = muc;
  foreach_dimension() 
    periodic (left); 
  init_grid (512);
  run();
}

event properties (i++){
  foreach_face()
    muc.x[] = fm.x[]/Re;
  boundary ((scalar *){muc});
}

event init(t = 0){
  init_lamb(1, xo, yo);
  vertex scalar vs[];
  foreach_vertex()
    vs[] = sq(sq(2*(y - yse))) - (fabs(x) - xse) ;
  foreach_vertex() //Close the top boundary
    vs[] -= 100000*(y - 6)*(y > 6) ;
  fractions(vs, cs, fs);
  boundary({cs});
  foreach(){
    u.x[] *= cm[];
    u.y[] *= cm[];
  }
  TOLERANCE = 1e-4;
  NITERMIN = 2;
  CFL = 0.6;
}

event adapt(i++){
  adapt_wavelet({cs, u.x, u.y}, (double[]){0.005, 0.04, 0.04}, maxlevel, 5);
}
/**
We generate a movie that displays the vortex dynamics.
*/

event mov (t += 0.1){
  scalar omega[];
  view(fov = 10);
  vorticity(u, omega);
  draw_vof ("cs", filled = -1, fc = {0.6, 0.4, 0.6});
  draw_vof ("cs");
  squares("omega", min = -7, max = 7, map = cool_warm);
  save("nja.mp4");
}

event end (t = 25);

/**
##Result

After the initial vortex rebound, the original dipole halves recombine
into a small vortex pair. This pair is able to pass through the
opening. 

![Such exotic dynamics](rebound/nja.mp4)
 
This can be compared against the Numerical simulations, a Point-Vortex model and experimental investigations,

![Numerical, point vortex (van Hooft) and experimental (Wong) results](http://www.basilisk.fr/sandbox/Antoonvh/rebound.reboundimg.png)
 
## References
* [C.H. Wong, *Electromagnetically Generated Dipolar Vortices in Shallow Stratified Fluids*](https://research.tue.nl/en/studentTheses/electromagnetically-generated-dipolar-vortices-in-shallow-stratif) 

* [J.A. van Hooft, *The Dynamical Behaviour of a Dipolar Vortex near Sharp-Edged Boundaries*](https://research.tue.nl/en/studentTheses/the-dynamical-behavior-of-a-dipolar-vortex-near-sharp-edged-bound)
 */
