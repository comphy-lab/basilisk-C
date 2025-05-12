/**
# Gradient sharpening by a dipolar vortex

A dipolar vortex propages itself though a modest scalar
gradient. hereby the scalar gradients sharply increase at the edge of
the dipole's atmosphere.

![The process of grdient sharpening](gsc/movv.mp4)(width=90%)

![Zoom on the gradient and grid](gsc/movvc.mp4)(width=70%)
*/

#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "view.h"

scalar s[], *tracers = {s};

double Re;
double Pr = 1.;
int maxlevel = 99; //Unlimited
FILE * fp;
double se = 0.02, ue = 0.02;

/**
The equations of motion are solved in to frame of reference that is
co-moving with the initialized dipolar vortex. As such, in and outflow
conditions are set at the left and right wall, respectively.
*/
u.n[left] = neumann (0.);
p[left]   = neumann (0.);
pf[left]  = neumann (0.);

u.n[right] = dirichlet (-1.);
p[right]   = dirichlet (0.);
pf[right]  = dirichlet (0.);

s[right] = dirichlet (tanh((x - 5. + t)/2.));

int main () {
  L0 = 25;
  X0 = Y0 = -L0/2;
  N  = 512;
  TOLERANCE = 1e-5;
  NITERMIN = 2;
  foreach_dimension()
    u.x.refine = refine_linear;
  for (Re = 50; Re <= 800; Re *= 2)
    //for (ue = 0.03; ue >= 0.01; ue -= 0.01)
    //for (se = 0.03; se >= 0.01; se /= 0.01)
    run();
}

double xo = -0.01, yo = 0.1;
#define RAD (sqrt(sq(x - xo) + sq(y - yo)))
#define ST ((yo - y)/RAD)
event init (t = 0) {
  char fname[99];
  sprintf (fname, "datac%g,%g%g", Re, se, ue);
  fp = fopen (fname, "w");
  scalar psi[];
  psi[left]  = dirichlet ((1/RAD - RAD)*ST);
  psi[right] = dirichlet ((1/RAD - RAD)*ST);
  psi[top]  = dirichlet ((1/RAD - RAD)*ST);
  psi[bottom] = dirichlet ((1/RAD - RAD)*ST);
  const face vector nu[] = {1./Re, 1./Re};
  mu = nu;
  double k = 3.83170597;
  refine (RAD < 2. && level < 10);
  refine (RAD < 1.1 && level < 11);
  foreach() 
    psi[] = ((RAD > 1)*(1/RAD - RAD)*ST +
	     (RAD <= 1)*(-2*j1(k*RAD)*ST/(k*j0(k))));
  boundary ({psi});
  foreach() {
    u.x[] = -((psi[0, 1] - psi[0, -1])/(2*Delta));
    u.y[] = (psi[1, 0] - psi[-1, 0])/(2*Delta);
    s[] = tanh((x - 5)/2);
  }
  boundary ({s, u});
}

event tracer_diffusion (i++) {
  const face vector kappa[] = {1./(Pr*Re), 1./(Pr*Re)};
  diffusion (s, dt, kappa);
}

event adapt (i++) {
  int minlevel = 6;
  adapt_wavelet ({s, u}, (double[]){se, ue, ue}, 99, minlevel);
  unrefine (x > X0 + 14*L0/15 && level > minlevel - 1);
}

event output (t += 0.1) {
  scalar lev[], gr[], omega[];
  face vector grf[];
  vertex scalar sv[];
  double maxgr = 0;
  boundary({s});
  foreach_face() {
    if (fabs(x - X0) < Delta/2. || fabs(x - X0 - L0) < Delta/2. ||
	fabs(y - Y0) < Delta/2. || fabs(y - Y0 - L0) < Delta/2.)
      grf.x[] = 0;
    else
      grf.x[] = (15.*(s[] - s[-1]) - (s[1] - s[-2]))/(12.*Delta);
  }
  boundary ((scalar*){grf});
  foreach_vertex() 
    sv[] = (s[] + s[-1] + s[0,-1] + s[-1,-1])/4.;
  foreach() {
    lev[] = level;
    gr[] = 0;
    foreach_dimension()
      // gr[] += sq((s[1] - s[-1])/(2.*Delta)); //sconder order diagnosis
      gr[] += sq((grf.x[] + grf.x[1])/2.); 
    if (gr[] > sq(maxgr))
      maxgr = sqrt(gr[]);
  }
  vorticity (u, omega);
  boundary ({gr, omega, sv});
  /**
For $Re = 800$ a movie is generated.
   */
  if (Re == 800) {
    view (fov = 25, width = 1200, height = 700);
    translate (x = t - 5) {
      squares ("omega",  linear = true, map = cool_warm, min = -8, max = 8);
      for (double sval = -0.8 ; sval <= 0.8; sval += 0.4)
	isoline ("sv", sval, lw = 2);
      //cells();
      save ("movv.mp4");
    }
    
    view (fov = 4, width = 1200, height = 700);
    squares ("gr", linear = true, map = cool_warm);
    cells();
    save ("movvc.mp4");
  }
    fprintf (fp, "%g %d %g %g %d %ld %g \n",
	     t, i, maxgr, statsf(gr).sum, grid->maxdepth, grid->n, perf.t);
    fflush (fp);
    printf ("%g %d %d %d %g %g %d\n",
	    t, i, mgp.i, mgpf.i, maxgr, statsf(gr).sum, depth());
    
}

event stop (t = 15) {
  char dname[99];
  sprintf (dname, "dump%g", Re);
  dump (dname);
  fclose (fp);
}
