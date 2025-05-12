#include "shakti.h"

#define MAXLEVEL 8
#define MINLEVEL 4

double rhow = 1000., rhoi = 917., g = 9.81, nu = 1.787e-6, n = 3., A = 2.5e-25, G = 0.05, om = 0.001, L = 334000., drag = 20.;
double ieb = 1.e-11, dzb = 0.02, hi = 500.;

scalar f[];
scalar ct[];
scalar * tracers = {f,ct};

/**outflow pressure boundary*/
p[left] = dirichlet(0);

/**here we set the grid size */
int main()
{
  size(64); 
  init_grid (1 << MAXLEVEL);
  run();
}

event init (i = 0) {
  foreach() {
    f[] = 0.01;
    ct[] = 0;
    }
}

event logfile (i++)
{
  stats s = statsf (f);
  fprintf (stderr, "%d %g %d %g %g %g\n", 
	   i, t/(24.*3600.), mgp.i, s.sum, s.min, s.max);
}

event movies (t += 2.592e5)
{
  output_ppm (f, min = 0.0007, max = 0.001, file = "f.mp4");
  
  output_ppm (p, min = 0, max = 120, file = "p.mp4");

  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, min = 0, max = MAXLEVEL, linear = false, file = "level.mp4");
}

event output (t += 2.592e5; t <= 12.96e6)
{
  static int nf = 0;
  printf ("file: field-%d\n", nf++);
  output_field ({p,f}, stdout, N);
}

event coefficients (i++)
{
  foreach_face() {
    double ff = (f[] + f[-1])/2.;
    double Re = om*norm(u)/nu;
    beta.x[] = clamp(-ff*ff*ff/(1.+Re)*g/12./nu,-1000.,-1.e-10);
  }
  foreach(){
    double uxx = (u.x[]+u.x[1,0])/2;
    double uyy = (u.y[]+u.y[0,1])/2;
    double q = sqrt(sq(uxx)+sq(uyy));
    double N = rhoi*g*hi + rhow*g*(dzb*x - p[]);
    double col = A*N*N*N*f[];
    double mdiff = ((m[1,0]*f[1,0]+m[]*f[])*(f[1,0]-f[])-((m[-1,0]*f[-1,0]+m[]*f[])*(f[]-f[-1,0])))/sq(Delta)/2 + ((m[0,1]*f[0,1]+m[]*f[])*(f[0,1]-f[])-((m[0,-1]*f[0,-1]+m[]*f[])*(f[]-f[0,-1])))/sq(Delta)/2;
    m[] = (G + 12.*nu*rhow*sq(q)/(f[]*f[]*f[])*(1.+om*q/nu))/L+mdiff;
    double moulin = clamp(ct[]/3600./24./30,0.,1.)*4.77464829276*exp(-(sq(x-16.015625)+sq(y-32.015625))/2.0);
    zeta[] = (1./rhow-1./rhoi)*m[] + col + ieb + moulin;
    opening[] = m[]/rhoi - col;
    ones[] = 1.;
    }
}

#if TREE
event adapt (i++) {
  double tolerance = 0.000005;
  adapt_wavelet ({f}, &tolerance, MAXLEVEL,MINLEVEL);
  event ("coefficients");
}
#endif

