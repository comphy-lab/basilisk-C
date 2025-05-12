/**
# Roll-wave Simulation VOF version (no quadtree adaptivity)

A trial to reproduce Brock's experimental results (1967). Code modified from [this example case](/src/example/gaussian-ns.c)
Not using the Momentum-conserving advection of velocity <span style="color:red"> **necessary???** </span>*/

//#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "reduced.h"
#include "view.h"
#include "navier-stokes/perfs.h"
//#include "navier-stokes/conserving.h"

// #define NU 1e-2
#define NU 0.0e-3
#define MAXLEVEL 14

// problem-sepcific parameters
double So = 0.05011;
double normalDepth = 0.00798;
double normalVelocity = 1.0377;
double gravityCoeff = 9.81;
double disMag = 0.15; // large amplitude
double disPeriod = 0.933;
double simTime = 20.0;
double Lx = 16.34304;
double Ly = 1.0;

scalar f0[];

int main()
{
  size (Lx);
  init_grid(1 << 11);
  rho1 = 0.001;
  rho2 = 1.;
  mu2 = NU;
  mu1 = mu2*rho1/rho2/10.;
  G.y = - 9.81;
  G.x = 9.81*So;
  Z.y = normalDepth;
  
  for (scalar s in {f0})
    s.refine = s.prolongation = fraction_refine;
  
  run();
}

/** inlet: specify h and u?? Follow tangaroa example for inlet depth. */
u.n[left] = dirichlet(normalVelocity*(1.0+disMag*(y>normalDepth ? 1.0 : 0.0)*sin(2. * pi * t / disPeriod)));
p[left] = neumann(0);
pf[left] = neumann(0);
f[left] = f0[];

// outlet
u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

// u.n[embed] = dirichlet(0.);
// u.t[embed] = dirichlet(0.);

/** slip BC or no-slip BC--to comment this line or not */
u.t[bottom] = dirichlet(0.);

event init (i = 0)
{
  
//   refine ((y < Ly) && level < 14);
  mask(y > 0.016 ? top : none);
  refine ((y < 0.01596) && level < 14);
  
  fraction (f0, y - normalDepth);
  
  foreach() {
      f[] = f0[];
      u.x[] = normalVelocity;
    }
  boundary ({f, u.x});
  
//   vertex scalar phi[];
//   foreach_vertex()
//     phi[] = y - BA*exp(- sq(x - 10.)/5.) - 1e-3;
//   boundary ({phi});
//   fractions (phi, cs, fs);
}

// this part related to subsonic inlet?
// event update_hl (i++)
// {
//   hl = 0.;
//   foreach_boundary (left, reduction(+:hl))
//     hl += Delta*f[];
//   hl = L0 - hl;
//   printf ("%g %g\n", t, hl);
// }

// event snapshot (i += 10; t <= 70) {
//   p.nodump = false;
//   dump();
// }

// to restrict dt
event maxdt (t <= simTime; t += 0.020);

event outputFile (t += 0.20)
{
  char name1[25];
  char name2[25];

  sprintf(name1, "out-prof-%.2f.txt", t);
  sprintf(name2, "out-field-%.2f.txt", t);

  FILE *fp1 = fopen(name1, "w");
  FILE *fp2 = fopen(name2, "w");

  output_facets (f, fp1);
  foreach ()
  {
    fprintf(fp2, "%g %g %g %g %g\n", x, y, f[], u.x[], u.y[]);
  }
}

// event pictures (t = end)
// {
//   view (fov = 4.04484, tx = -0.498476, ty = -0.0923365, sy = 5,
// 	bg = {1,1,1},
// 	width = 1869, height = 390);
//   draw_vof ("cs", filled = -1, fc = {1,1,1});
//   draw_vof ("f", filled = 1, fc = {1,1,1});
//   squares ("u.x", min = -0.5, max = 4, linear = true);
//   isoline ("u.x", 0., lc = {1,1,1}, lw = 2);
//   save ("u.x.png");
// 
//   draw_vof ("cs", filled = -1, fc = {1,1,1});
//   draw_vof ("f", filled = 1, fc = {1,1,1});
//   squares ("u.y", min = -0.8, max = 0.8, linear = true);
//   save ("u.y.png");  
// }
