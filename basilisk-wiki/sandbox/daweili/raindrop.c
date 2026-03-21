/**
# Films do not break (on GPUs either)(Copy from fuster)

![Inside distance function](raindrop/movie.mp4)(width="80%") 

~~~gnuplot
set ylabel 'Minimal distance'
set xlabel 'Time'
set logscale y
set grid
plot 'log' u 1:($2*2.) w l t ''
~~~
*/

#include "grid/gpu/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase-clsvof.h"
#include "integral.h"

const double tend = 15.;
const double We = 5.;
const double Re = 1000;
const int maxlevel = 11;

double dmin0;

FILE * fp;

int main()
{
  periodic(right);

  dmin0 = 1.;
  L0 = 10.;
  X0 = - L0/2.;
  Y0 = - L0/2. - 3.;
  const scalar sigma[] = 1./We;
  d.sigmaf = sigma;

  fp = fopen("dmin.dat","w");
  rho1 = 1, rho2 = 0.1;
  mu1 = mu2 = 1./Re;

  N = 1 << maxlevel;
  run();
  
  fclose(fp);
}

event acceleration (i++)
{
  foreach_face (y)
    a.y[] -= 1.;
}

event init (t = 0) {
  foreach()
    d[] = - sq(y + 0.01*sin(2*pi*x/L0*2)) + sq(0.5);
}

event movies (i += 30; t <= tend)
{
  scalar df[];
  foreach()
    df[] = f[] > 0.5 ? d[] : nodata;
  output_ppm (df, spread = -1, file = "movie.mp4", n = 1024);
  //  output_ppm (df, spread = -1, file = "zoom.mp4", box = {{-1.25,-1.75},{1.25,2}});
}

scalar dmin[];

event logfile (i = 100; i += 5)
{
  foreach() {
    dmin[] = nodata;
    if (d[] > 0) {
      double dx = (d[1] - d[-1])/2.;
      double dy = (d[0,1] - d[0,-1])/2.;
      double dn = sqrt(1. + sq(dx) + sq(dy));
      double dxx = d[1] + d[-1] - 2.*d[];
      double dyy = d[0,1] + d[0,-1] - 2.*d[];
      double dxy = (d[1,1] - d[-1,1] - d[1,-1] + d[-1,-1])/4.;
      // First fundamental form
      double E = 1. + sq(dx), F = dx*dy, G = 1. + sq(dy);
      // Second fundamental form (also noted L, M, N)
      double e = dxx/dn, f = dxy/dn, g = dyy/dn;
      // Gaussian curvature
      double K = (e*g - f*f)/(E*G - F*F);
      // Mean curvature
      double H = (e*G - 2.*f*F + g*E)/(2.*(E*G - F*F));
      // Principal curvatures
      double /* k1 = H + sqrt(sq(H) - K),*/ k2 = H - sqrt(sq(H) - K);
      if (k2 < - Delta) {
	// Principal direction i.e. dy/dx
	double lambda = - (f - k2*F)/(g - k2*G);
	// the derivative cancels at (s,lambda*s)
	double s = - (dx + lambda*dy)/(dxx + sq(lambda)*dyy + 2.*lambda*dxy);
	// distance where the derivative cancels
        if (y > 0.)
	dmin[] = d[] + s*(dx + lambda*dy) + sq(s)/2.*(dxx + sq(lambda)*dyy + 2.*lambda*dxy);
        else
          dmin[]=L0;
      }
    }
  }
  
  double dmintmp = statsf(dmin).min;
  
  if (dmintmp/dmin0 > 0.1) {
    dmin0 = dmintmp;
    fprintf (stderr, "%g %g\n", t, dmintmp );
    fprintf (fp, "%g %g\n", t, dmintmp );
  }
  output_ppm (dmin, spread = -1, file = "dmin.mp4", n = 1024);
}