/**
# The flow towards a cone
It's a simple test of a uniform inflow from a channel towards a triangle, we want to check the flow dynamics with different shape of exit shape and triangle position.
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "view.h"
/**
We defind some MACRO here.

The origin (0,0) is set to the lower conner of the exit, with channel height $diameter$, its length (from left border to exit) $aa1$ and the minimum exit wall thickness at bottom $bb1$ (used for $L_0 = 1$).

The distance of traingle's left conner towards exit $conepx$, 
with a possible height differece $eta$ between the triangle's lower interface and the one of channel. */
#define L_p 1 // domain size L0 = L_p^2
#define LVM 8 + L_p 
#define aa1 0.3  
#define bb1 0.4  
#define diameter 0.15  
#define eta 0.02  
#define conepx 0.1  
/**
the NORMALIZED (with $L_0 = 1$) constante viscosity mu = 1/Re in the simulation.
*/
#define mucst 1./2000   
#define tend 10.0  
#define ADAPT 1
// DO NOT SET 2 SWITCH OPEN AT SAME TIME,round will override chamber
#define round_conner 0
#if !round_conner
#define chamfer_conner 0  
#endif

/**
we set a tracer for fun.
*/
scalar f[];
scalar * tracers = {f};
face vector muv[];

/**
## main
*/
int main() {
  size(pow(2,L_p));
  init_grid (64);
  origin( -aa1, -bb1-(L0-1.)/2. );
  mu = muv;
  TOLERANCE = 1e-3;
  run(); 
}

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*mucst;
}
/**
The fluid is injected on the left boundary with a unit velocity. 
*/

u.n[left]  = dirichlet(1.*cs[]);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet( y < eta || (y < eta+diameter/2. && y > diameter/2.) );

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);
/**
Free outflow condition at other three boundaries.
Comment the follow for using free slip conditions.
*/
u.n[top] = neumann(0.);
p[top]   = dirichlet(0.);
pf[top]  = dirichlet(0.);

u.n[bottom] = neumann(0.);
p[bottom]   = dirichlet(0.);
pf[bottom]  = dirichlet(0.);


/**
ALL the solid surfaces are no-slip. */

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

/**
## init event 
We use the function proposed in Basilisk

   #define intersection(a,b)   min(a,b)

   #define union(a,b)          max(a,b)

    #define difference(a,b)     min(a,-(b))

    phi[] > 0 means the fluid domain.

with fraction field $cs$, fluids : cs[]=1 and solid : cs[]=0. 
*/

event init (t = 0)
{
  refine (( x<1.-aa1 && x > -aa1 && y < 1.-bb1 && y > -bb1 ) && level < LVM ) ;
  vertex scalar phi[];
  vertex scalar rec1[];
  vertex scalar rec2[];
  vertex scalar cone[];

  /**
A rounded conner shape exit.

We use "simply" two rect and one circle!!
*/
  #if round_conner
  double rr1 = 0.1; //the bottom
  vertex scalar rec11[];
  vertex scalar rec12[];
  vertex scalar rec1t[];
  double rr2 = 0.1; //the top
  vertex scalar rec21[];
  vertex scalar rec22[];
  vertex scalar rec2t[];
  vertex scalar cir2[];

  foreach_vertex() {
    rec11[] = union (x, y+rr1);
    rec12[] = union (x+rr1, y);
    rec1t[] = intersection (rec11[], rec12[]);
    //left lower rectangle <0, rest >0
    rec1[] = intersection (rec1t[],  sq(x+rr1) + sq(y+rr1)-sq(rr1));

    rec21[] = union (x+rr2, diameter - y);
    rec22[] = union (x, diameter + rr2 - y);
    rec2t[] = intersection (rec21[], rec22[]);
    //left upper rectangle <0, rest >0
    rec2[] = intersection (rec2t[], sq(x+rr2) + sq(y-rr2-diameter)-sq(rr2));
  }
  /** 
A chamfer conner shape exit.

We use "simply" two trapezoids at angle 45 degree!!
*/
  #elif chamfer_conner
  double tt1 = 0.1; //the bottom
  vertex scalar trap11[];
  vertex scalar trap12[];
  double tt2 = 0.1; //the top
  vertex scalar trap21[];
  vertex scalar trap22[];

  foreach_vertex() {
    trap11[] = union (x, y+x+tt1);
    trap12[] = union (y, y+tt1+x);
    //left lower rectangle <0, rest >0
    rec1[] = union (trap11[], trap12[]);

    trap21[] = union (x, diameter + tt2 + x- y);
    trap22[] = union (diameter-y, diameter + x + tt2 - y );
    //left upper rectangle <0, rest >0
    rec2[] = union (trap21[], trap22[]);
  }
  /** 
A normal rectangle shape exit.
*/
  #else
  foreach_vertex() {
    //left lower rectangle <0, rest >0
    rec1[] = union (x, y);
    //left upper rectangle <0, rest >0
    rec2[] = union (x, diameter - y);
  }
  #endif 

  /**
With the choosen exit shape, now we define the traingle shape and remove the 2 solid exit domain form the domain.
*/
  foreach_vertex() {
    //right cone <0, rest >0,assume slope = 0.3
    cone[] = union ( (y - 0.3*( x - conepx )) ,eta-y);
    // remove rec1 from cone
    phi[] = intersection (cone[],rec1[]);
    // remove rec1 from phi
    phi[] = intersection (phi[], rec2[]);
  }

  boundary ({phi});
  fractions (phi, cs, fs);

  //foreach()
  //	u.x[] = cs[] ? 1. : 0.;

  // ============ whole domain view=================
  char legend[1000];
  sprintf(legend, "      mu=%g, D=%g, eta=%g, t=%0.2g", mucst, diameter, eta, t);
  view (fov = 24., tx = aa1/L0 - 0.5, ty =  (bb1-0.5)/L0, bg = {1,1,1}, width = 1024, height = 1024, samples = 1);
  clear();
  box();
  cells();
  draw_string(legend, 1, size = 50.,lw = 3.);
  draw_vof ("cs", filled = -1, fc = {1,1,1});
  squares("cs",max=1.,min=0.);
  #if round_conner
  save ("round.png");
  #elif chamfer_conner
  save ("chamber.png");
  #else
  save ("rectangle.png");
  #endif

  // ============ local view=================
  view (fov = 24./L0, tx = (aa1-0.5)/L0, ty = (bb1-0.5)/L0, bg = {1,1,1}, width = 1024, height = 1024, samples = 1);
  clear();
  cells();
  draw_string(legend, 1, size = 50.,lw = 3.);
  draw_vof ("cs", filled = -1, fc = {1,1,1});
  squares("cs",max=1.,min=0.);
  save ("test.png");

  //exit(0);
}

/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++){
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
}
/**
We output some videos. we neglect the "converged regime" before t = 1
*/

event movies (t += tend/300. ; t <= tend)
{
  if (t > 1.) {
    scalar omega[];
    vorticity (u, omega);
    char legend[2000]; //time
    sprintf(legend, "      mu=%g, D=%g, eta=%g, t=%0.2g", mucst, diameter, eta, t);


    view (fov = 24./L0, tx = (aa1-0.5)/L0, ty = (bb1-0.5)/L0, bg = {1,1,1}, width = 1024, height = 1024, samples = 1);
    clear();
    box();
    draw_vof ("cs", filled = -1, fc = {1,1,1});
    draw_string(legend, 1, size = 50.,lw = 3.);
    squares("omega");
    #if round_conner
    save ("omega_round.mp4");
    #elif chamfer_conner
    save ("omega_chamber.mp4");
    #else
    save ("omega_rectangle.mp4");
    #endif

    view (fov = 24., tx = aa1/L0 - 0.5, ty =  (bb1-0.5)/L0, bg = {1,1,1}, width = 1024, height = 1024, samples = 1);
    clear();
    cells();
    draw_vof ("cs", filled = -1, fc = {1,1,1});
    draw_string(legend, 1, size = 50.,lw = 3.);
    squares("f",linear = false, min = 0, max = 1,);
    #if round_conner
    save ("f_round.mp4");
    #elif chamfer_conner
    save ("f_chamber.mp4");
    #else
    save ("f_rectangle.mp4");
    #endif
  }
}

/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */
#if ADAPT
event adapt (i++) {
  adapt_wavelet ({cs,u,f}, (double[]){1e-2,3e-2,3e-2,3e-2}, LVM, 4);
}
#endif
