/* tangential interpolation on face vector fields  */

#include "utils.h"

scalar h[];
face vector u[];

static double solution (double x, double y, double z)
{
#if 0
  double R0 = 0.1;
  return exp(-(x*x+y*y)/sq(R0));
#else
  return dimension == 2 ? x - y : 2.*x - y - z;
#endif
}

double tolerance = 1e-4;

static int error()
{
  double max = 0, maxv = 0., maxw = 0.;

  scalar eu[];
  foreach(reduction(max:max) reduction(max:maxv) reduction(max:maxw))
    for (int i = 0; i <= 1; i++) {
      double xu = x + (i - 0.5)*Delta, yu = y, zu = z ;
      eu[] = fabs (solution(xu,yu,zu) - u.x[i,0]);
      if (eu[] > max)
	max = eu[];

      double xv = x, yv = y + (i - 0.5)*Delta, zv = z;
      double e = solution(xv,yv,zv) - u.y[0,i];
      if (fabs(e) > maxv)
	maxv = fabs(e);
#if DEBUG
      printf ("%g %g %g %d %d %g %g\n", 
	      xu, yu, zu, level, cell.neighbors, u.x[i,0], eu[]);
#endif

#if dimension == 2
      maxw = maxv;
#else // dimension == 3
      double xw = x, yw = y, zw = z + (i - 0.5)*Delta;
      e = solution(xw,yw,zw) - u.z[0,0,i];
      if (fabs(e) > maxw)
	maxw = fabs(e);
#endif      
    }

  fprintf (stderr, "maximum error: %g %g %g\n", max, maxv, maxw);
  stats s = statsf (eu);
  fprintf (stderr, "eu: avg: %g stddev: %g max: %g\n", 
	   s.sum/s.volume, s.stddev, s.max);

  return (max != maxv || max != maxw);
}

int main (int argc, char ** argv)
{
  int maxlevel = dimension == 2 ? 10 : 7;
  int n = 1 << maxlevel;
  init_grid (n);

  u.x.refine = refine_face_solenoidal;
  
  foreach_dimension() {
    u.n[right] = solution(x,y,z);
    u.n[left]  = solution(x,y,z);
  }

  foreach_dimension() {
    u.t[right] = dirichlet(solution(x,y,z));
    u.t[left]  = dirichlet(solution(x,y,z));
#if dimension > 2    
    u.r[right] = dirichlet(solution(x,y,z));
    u.r[left]  = dirichlet(solution(x,y,z));
#endif
  }

  origin (-0.5, -0.5, -0.5);
  double R0 = 0.1;
  foreach()
    h[] = exp(-(x*x + y*y + z*z)/sq(R0));

  foreach_face()
    u.x[] = solution(x,y,z);

  astats s = adapt_wavelet ({h}, &tolerance, maxlevel);
  while (s.nc) {
    fprintf (stderr, "refined: %d coarsened: %d\n", s.nf, s.nc);
    s = adapt_wavelet ({h}, &tolerance, maxlevel);
  }
  
#if DEBUG
  FILE * fp = fopen ("cells", "w");
  output_cells (fp);
  fclose (fp);
#endif

  if (error())
    return 1;

  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += u.x[] - u.x[1];
  }
  stats sdiv = statsf(div);
  fprintf (stderr, "div before: %g\n", sdiv.max);

  tolerance = 1e-5;
  s = adapt_wavelet ({h}, &tolerance, maxlevel, list = {h,u});
  fprintf (stderr, "refined: %d coarsened: %d\n", s.nf, s.nc);

  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += u.x[] - u.x[1];
  }
  sdiv = statsf(div);
  fprintf (stderr, "div after: %g\n", sdiv.max);

#if DEBUG
  fp = fopen ("cells1", "w");
  output_cells (fp);
  fclose (fp);
#endif
  
  if (error())
    return 1;
}
