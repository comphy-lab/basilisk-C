/**
# Bubble rising in a elasto-viscoplastic fluid

This files presents the case of a bubble rising in a elasto-viscoplastic material (Saramito material).

## Problem nondimensionalization


## Code
*/


#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#include "saramito.h"

#define Pl 0.04         // Plastic number
#define Ar 2.           // Archimedes number
#define Bo 10.          // Bond number
#define Deb 10.         // Deborah number
#define MUr 0.01        // Ratio of outer(matrix) to inner(bubble) viscosities
#define RHOr 0.01       // Ratio of outer to inner densities
#define Beta 0.5        // Ratio of the solvent visc. to the total viscoelastic visc.
#define R0 0.5          // Bubble radius
#define A0 0.793701     // For an elliptical bubble 
#define B0 0.39685
#define HEIGHT 20*R0    // Domain size
#define X0 5*R0         // Bubble initial position 


// Mesh levels
int LEVEL = 9;
int MIDLEVEL = 5;
int MINLEVEL = 4;

double gamma_c, muc1, mupc1;
double D = 2*R0;         // Droplet diameter
double U = 1.;           // Characteristic velocity U = V * D / H
double nnn = 0.7;        // Flow index
double Nreg = 1.e0;      // Dimensionless regularization parameter.
double tEnd = 10.;
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

  /**
  We set the viscosities, densities, surface tension and visco-elastic
  parameters. */
  
  gamma_c = U/D;                                   // Characteristic strain rate
  muc1 = 1/Ar;                                     // Matrix characteristic viscosity
  mupc1 = (1-Beta)*muc1;                           // Polymeric characteristic viscosity of the matrix
  tau_y = Pl*mupc1*gamma_c;                        // Yield stress
  K = (mupc1-tau_y/gamma_c)/pow(gamma_c, nnn-1);   // Powe law index
  epsilon = tau_y/(Nreg * muc1);                   // Regularization parameter

  mu1 = Beta*muc1;                                 // Solvent viscosity of the outer phase
  mu2 = MUr/Ar;                                    // Inner phase (bubble) Newtonian viscosity
  rho1 = 1.;                                       // Outer phase density
  rho2 = RHOr;                                     // Inner phase density
  f.sigma = 1./Bo;                                 // Surface tension
  
  lamb_c = Deb;
  mupp = mupc1;
  fi = 0.7;


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
// Spherical bubble
//    refine (sq(x - X0) + sq(y) - sq(1.2*R0) < 0 && level < LEVEL);
//    fraction (f, - sq(R0) + (sq(x - X0) + sq(y)));

// Elliptical bubble    
    refine (sq((x - X0)/(1.2*A0)) + sq(y/(1.2*B0)) - 1 < 0 && level < LEVEL);
    fraction (f, sq((x - X0)/A0) + sq(y/B0) - 1);
  }
}

/*
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
  printf ("i = %d t = %g, # = %d, l.mi = %g, l.ma = %g ep = %g tau_y = %g K = %g mu1 = %g mu2 = %g\n", i,t,j, s.min, s.max, epsilon, tau_y, K, mu1, mu2);
  fflush(stdout);
}*/



event adapt (i++)
{
  adapt_wavelet ({f, u, yielded, mupc}, (double[]){1e-3, 1e-1, 1e-1, 1e-1, 1e-1}, LEVEL, MINLEVEL);
}


/*
event snapshot (t = 0; t <= tEnd; t += 1.)
{
  char named[80];
  sprintf (named, "dump-%.2f", t);
  dump (file = named);
}*/


/**
Bubble volume, position and velocity.
*/
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
  output_field ((scalar *){u, f, mupd, lam, tau_p, yielded}, fld2, n = 500, linear = true,  box = {{0., 0.},{10.,10.}});
  fclose (fld2);
}*/




event movies (i += 10)
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
  squares("mupc", linear = true);
    mirror (n = {0,-1})
        {
          draw_vof ("f");
          squares("mupc", linear = true);
        }
  save("movie_mup.mp4");

  scalar tau_pp[];
  foreach()
  {
      tau_pp[] = pow(0.5*(sq(tau_p.x.x[])+sq(tau_p.y.y[])+sq(tau_p.x.y[])+sq(tau_p.y.x[])),0.5);
  }
  boundary({tau_pp});
  draw_vof ("f");
  squares("tau_pp", linear = true);
  mirror (n = {0,-1})
        {
          draw_vof ("f");
          squares("tau_pp", linear = true);
        }
  save("movie_tau_p.mp4");

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

### Velocity field
![Velocity field](bubble_rising_saramito/movie_vel.mp4)

![Polymeric viscosity field](bubble_rising_saramito/movie_mup.mp4)

![Polymeric Stress field](bubble_rising_saramito/movie_tau_p.mp4)

![Yielded and unyielded regions](bubble_rising_saramito/movie_yield.mp4)

~~~gnuplot
set xlabel 'Time'
set ylabel 'Bubble Velocity'
set grid
plot 'pos.txt' u 1:4 w l lw 3
~~~
*/