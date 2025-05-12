/**
# Bubble rising in a elasto-viscoplastic fluid

Using the Mendes-Thompson model in a header file

## Problem nondimensionalization

## Code
*/


#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"

#include "mendes-thompson.h"

#define Pl 0.05                           // Plastic number
#define Ar 2.                             // Reynold number
#define Bo 5.                             // Weber number
#define Deb 0.5                           // Deborah number
#define MUr 0.01                          // Ratio of outer(matrix) to inner(drop) viscosity
#define RHOr 0.01                         // Ratio of outer to inner density
#define R0 0.5                            // Bubble radius
#define A0 0.793701
#define B0 0.39685
#define HEIGHT 20*R0                      // Domain size
#define X0 5*R0                           // Bubble initial position 
#define U 1.                              // Characteristic velocity
#define D 1.                              // The domain height
#define BETA 0.5                          // Ratio of the solvent viscosity to the total viscosity
#define ALPHA 1e2                         // The ratio of the totally unstructured to the totally structured material viscosity
#define IOTA 1.                           // The ratio of the characteristic elastic modulus to the fully structed material elastic modulus
#define density 1.                        // Density
#define MUC 1./Ar                         // Total viscoelastic viscosity 
#define MUS ((1. - BETA)*MUC)             // Structural viscosity
#define MUI (BETA*MUC)                    // Infinity viscosity 
#define MU0 (MUI*ALPHA)                   // Totally structured material viscosity
#define Lambda_c log(BETA)/log(1/ALPHA)   // Characteristic structural parameter
#define GC (U/D)*(MUS/Deb)                // Characteristic Elastic Modulus
#define G0 GC/IOTA                        // Fully structed material elastic modulus
#define MM (log(IOTA)+1)*Lambda_c         // Constant of the elastic modulus equation
#define t_eq 1e-2                         // Equilimbrium time

int LEVEL = 9;
int MIDLEVEL = 5;
int MINLEVEL = 4;

double gamma_c;
double nnn = 0.7;                         // Flow index
double Nreg = 1.e0;                       // Dimensionless regularization parameter.

double tEnd = 0.3;

double xdd = 0.;
double vdd = 0.;
double uddx = 0.;

u.t[top] = dirichlet (0.);
u.t[left] = dirichlet (0.);
u.t[right] = dirichlet (0.);


int main() 
{
  size (HEIGHT);
  init_grid (1 << MIDLEVEL);
  //  DT = 1e-4;
  
  gamma_c = U/D;                                 // Characteristic strain rate
  mu1 = MUI;                                     // Fully unstructured material viscosity
  rho1 = 1.;                                     // Outer density

  tau_y = Pl*MUS*gamma_c;                        // Yield stress
  K = (MUS-tau_y/gamma_c)/pow(gamma_c, nnn-1);   // Consistency index
  fi = nnn;                                      // Flow index
  mui = MUI;                                     // Fully unstructured material viscosity
  mu0 = MU0;                                     // Fully structured material viscosity
  g0 = G0;                                       // Fully structured material elastic modulus
  mm = MM;                                       // Positive constant of the structure elastic modulus model
  t_eqq = t_eq;                                  // Thixotropic equilibrium time

  mu2 = MUr*MUC;                                 // Inner phase Newtonian phase viscosity
  rho2 = RHOr;                                   // Inner phase density

  f.sigma = 1./Bo;                               // Surface tension

  TOLERANCE = 1e-4;
  NITERMAX = 100;
  run();
}


event acceleration (i++)
{
  face vector av = a;
  foreach_face(x)
    av.x[] = -1.;
}


event init (i = 0) 
{
  if (!restore (file = "restart"))
  {
//    refine (sq(x - X0) + sq(y) - sq(1.2*R0) < 0 && level < LEVEL);
//    fraction (f, - sq(R0) + (sq(x - X0) + sq(y)));

    refine (sq((x - X0)/(1.2*A0)) + sq(y/(1.2*B0)) - 1 < 0 && level < LEVEL);
    fraction (f, sq((x - X0)/A0) + sq(y/B0) - 1);
  }

  foreach()
  {
    Lambda[] = 1.;
    Lambda_eq[] = 1.;
    mus[] = (MU0-MUI)*f[];
    G[] = G0;
    lamb[] = f[]*mus[] / G[];
  }
  boundary ((scalar *){Lambda, Lambda_eq, mus, G, lamb});
}




scalar le[];
event logfile(i+=1)
{
  int j = 0;
  foreach()
  {
    le[] = level;
    j++;
  }
  stats s = statsf(le);

  printf ("i = %d t = %g, # = %d, l.mi = %g, l.ma = %g tau_y = %g K = %g mu1 = %g mu2 = %g MUI = %g MUS = %g MU0 = %g G0 = %g mm = %g\n", i, t, j, s.min, s.max, tau_y, K, mu1, mu2, MUI, MUS, MU0, G0, mm);
  fflush(stdout);
}


event adapt (i++)
{
  adapt_wavelet ({f, u, yielded, mus}, (double[]){1e-3, 1e-1, 1e-1, 1e-1, 1e-1}, LEVEL, MINLEVEL);
}


/*
event snapshot (t = 0; t <= tEnd; t += 0.5)
{
  char named[80];
  sprintf (named, "dump-%.2f", t);
  dump (file = named);
}*/


double xdd;
double vdd;
double uddx;
event saveDatas (i++)
{
  double xd1 = 0.;
  double vd1 = 0.;
  double udx1 = 0.;

  foreach(reduction(+:vd1) reduction(+:xd1) reduction(+:udx1))
  {
     vd1 += dv()*(1.-f[]);
     xd1 += dv()*(1.-f[])*x;
     udx1 += u.x[]*dv()*(1.-f[]);
  } 

  xdd = xd1;
  vdd = vd1;
  uddx = udx1;

  char nameu2[50];
  sprintf (nameu2, "pos.txt");
  static FILE * fname2 = fopen (nameu2, "w");
  fprintf (fname2, "%.8f %.8f %.8f %.8f\n", t, vd1, xd1/vd1, udx1/vd1);
  fflush (fname2);
}


/*
event uprofile (t += 1.0) 
{
  char namei1[80];
  sprintf (namei1, "interface-%.2f.txt", t);
  FILE* fp1 = fopen (namei1, "w");
  output_facets (f, fp1);
  fclose (fp1);

  char field2[80];
  sprintf (field2, "fields-output2-%.2f.txt", t);
  FILE* fld2 = fopen (field2, "w");
  output_field ((scalar *){u, f, mus, lam, tau_p, yielded}, fld2, n = 500, linear = true,  box = {{0., 0.},{10.,10.}});
  fclose (fld2);
}*/




event movie_mesh (i += 10)
{
  clear(); 
  view(fov = 20., tx = -0.5, ty = 0.001);

  scalar vel[];
  foreach()
    vel[] = pow(sq(u.x[])+sq(u.y[]), 0.5);
  boundary({vel});
  draw_vof ("f");
  squares("vel", linear = true);
  mirror (n = {0,-1})
        {
          draw_vof ("f");
          squares("vel", linear = true);
        }
  save("movie_vel.mp4");

  draw_vof ("f");
  squares("mus", linear = true);
    mirror (n = {0,-1})
        {
          draw_vof ("f");
          squares("mus", linear = true);
        }
  save("movie_mus.mp4");

  draw_vof ("f");
  squares("Lambda", linear = true);
    mirror (n = {0,-1})
        {
          draw_vof ("f");
          squares("Lambda", linear = true);
        }
  save("movie_Lambda.mp4");

  draw_vof ("f");
  squares("Lambda_eq", linear = true);
    mirror (n = {0,-1})
        {
          draw_vof ("f");
          squares("Lambda_eq", linear = true);
        }
  save("movie_Lambda_eq.mp4");

  draw_vof ("f");
  squares("stress_norm_dev", linear = true);
  mirror (n = {0,-1})
        {
          draw_vof ("f");
          squares("stress_norm_dev", linear = true);
        }
  save("movie_stress.mp4");

  draw_vof ("f");
  squares("yielded", linear = true);
  mirror (n = {0,-1})
        {
          draw_vof ("f");
          squares("yielded", linear = true);
        }
  save("movie_yield.mp4");
}


event end (t = tEnd) {}



/**
## Results

![Velocity field](bubble_rising_mendes-thompson/movie_vel.mp4)

![Structural viscosity field](bubble_rising_mendes-thompson/movie_mus.mp4)

![Total Stress field](bubble_rising_mendes-thompson/movie_stress.mp4)

![Yielded and unyielded regions](bubble_rising_mendes-thompson/movie_yield.mp4)

![Structural parameter](bubble_rising_mendes-thompson/movie_Lambda.mp4)

![Equilibrium structural parameter](bubble_rising_mendes-thompson/movie_Lambda_eq.mp4)

~~~gnuplot
set xlabel 'Time'
set ylabel 'Bubble Velocity'
set grid
plot 'pos.txt' u 1:4 w l lw 3
~~~
*/