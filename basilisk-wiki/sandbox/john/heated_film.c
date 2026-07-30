/**
# Evaporating film on an inclined, non-uniformly heated wall

Liquid film ($f=1$) flowing on a heated wall at the bottom, with
passive gas above. $x$ along the incline and $y$ wall-normal.Momentum equation is:
$$
\rho\,(\partial_t\mathbf u + \mathbf u\!\cdot\!\nabla\mathbf u) =
-\nabla p + \nabla\!\cdot\!(2\mu\mathbf D)
+ \sigma_0\kappa\,\mathbf n\,\delta_s + \nabla_s\sigma\,\delta_s
+ \rho_1 f\,g\sin\beta\,\mathbf e_x - \rho\,g\cos\beta\,\mathbf e_y
$$


Energy with latent heat, and the one-sided kinetic law:
$$
\rho c_p \frac{DT}{Dt} = \nabla\!\cdot\!(k\nabla T) + J h_{lv}\delta_s,
\qquad J = z\,(T_{sat} - T_i)
$$

*/
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"     // constant-sigma normal (Laplace) force
#include "tracer.h"      // advection of T
#include "diffusion.h"   // implicit diffusion of T
#include "view.h"


/**
## Parameters
*/
double h0 = 1.;                    // flat-film thickness
double beta0 = 45.*pi/180.;        // inclination angle
double k1 = 9.52e-3, k2 = 5e-4;    // conductivities
double cp1 = 1., cp2 = 1.;
double Tsat = 0.;                  // saturation temperature
double DT0 = 1., DT1 = 0.8;        // mean superheat + lateral modulation for **TWALL**
double sigma0 = 3.;                // constant surface tension
double sigmaT = 0.03;              // -d(sigma)/dT
double zc  = 1e-2;                 // kinetic coefficient of J = zc*(Tsat-T)
double hlv = 20.;                  // latent heat
double qin = 1.;                   // run time computable inlet liquid flux
double V0 = 0., tauV = 8.;         // reference volume
double gval;

int LEVEL = 8, MINLEVEL = 4;
double TEND = 5.;


/* wall temperature distribution */
#define TWALL(x) (Tsat + DT0 + DT1*cos(2.*pi*(5.*x)/L0))

/* evaporation flux */
#define JFLUX(TT) (zc*(Tsat - (TT)))

/* inlet velocity to optain correct qin*/
#define UPROFILE(y) ((y) < h0 ? h0*(y) - sq(y)/2. : sq(h0)/2.)

/* inlet temperature at t=0*/
#define TI_0 (Tsat + (TWALL(0.) - Tsat)/(1. + zc*hlv*h0/k1))

/* inlet temperature BC */
#define TPROFILE(y) ((y) < h0 ? TWALL(0.) + (TI_0 - TWALL(0.))*(y)/h0 : TI_0)

/* volume fraction tolerance */
#define F_TOL 1e-6


scalar T[];
scalar * tracers = {T};

u.t[bottom] = dirichlet (0);         // no-slip substrate
T[bottom]   = dirichlet (TWALL(x));  // heated/patterned wall
T[top]      = neumann (0);           // passive gas far field

int main()
{
#if _OPENMP
  fprintf (stderr, "# OpenMP: running with %d threads\n",
	   omp_get_max_threads());
#endif

  L0 = 5.;
  rho1 = 1.,  rho2 = 1e-2;
  mu1 = 1./10., mu2 = mu1/50.;     // nu1 = 1/10 ->  Re = 15
  f.sigma = sigma0;
  gval = 3.*(mu1/rho1)/(sq(h0)*sin(beta0));  // Ubar === 1
  N = 1 << LEVEL;
  CFL = 0.4;

  /* inflow: h(0) = h0, flux qin, conduction-kinetic T profile */
  u.n[left] = dirichlet (qin*UPROFILE(y)/(cube(h0)/3.));
  u.t[left] = dirichlet (0);
  T[left]   = dirichlet (TPROFILE(y));
  p[left]   = neumann (0);
  pf[left]  = neumann (0);
  f[left]   = (y < h0);

  /* outflow */
  u.n[right] = neumann (0);
  T[right]   = neumann (0);
  p[right]   = dirichlet (0);
  pf[right]  = dirichlet (0);
  f[right]   = neumann (0);

  run();
}

/**
## Initial condition
Flat film with a small pertrubation */
event init (t = 0)
{
  fraction (f, h0*(1. + 0.05*cos(8.*pi*x/L0)) - y);

  double nu1 = mu1/rho1, gs = gval*sin(beta0);
  foreach() {
    double yl = min (y, h0);
    u.x[] = gs/nu1*(h0*yl - sq(yl)/2.);
    T[] = Tsat + (TWALL(x) - Tsat)*clamp(1. - y/h0, 0., 1.);
  }
}

/**
## Phase change
$\mathbf v_{pc} = (J/\rho_1)\,\mathbf n$ with
$\mathbf n = -\nabla f/|\nabla f|$ (liquid $\to$ gas) */
face vector uf_save[];

event stability (i++)
{
  vector pcv[];
  foreach() {
    pcv.x[] = pcv.y[] = 0.;
    if (f[] > F_TOL && f[] < 1. - F_TOL) {
      double gx = (f[1] - f[-1])/(2.*Delta);
      double gy = (f[0,1] - f[0,-1])/(2.*Delta);
      double gm = sqrt (sq(gx) + sq(gy));
      if (gm > 1e-10) {
		double J = JFLUX(T[]);
		pcv.x[] = - J/rho1*gx/gm;
		pcv.y[] = - J/rho1*gy/gm;
      }
    }
  }
  foreach_face() {
    uf_save.x[] = uf.x[];
    uf.x[] += fm.x[]*(pcv.x[] + pcv.x[-1])/2.;
  }
}

event tracer_advection (i++)
{
  foreach_face()
    uf.x[] = uf_save.x[];
}

/**
## Inlet flux
We fill the liquid that escapes and evaporates. Not exact but 
via interpolation.
*/
event inflow_control (i++)
{
  /* liquid outflow just inside the right boundary */
  double qo = 0., dy = L0/(1 << LEVEL), xq = L0 - 2.*dy;
  for (double y = dy/2.; y < 4.*h0; y += dy)
    qo += interpolate (f, xq, y)*interpolate (u.x, xq, y)*dy;

  /* net phase-change rate */
  double Jtot = 0., V = 0.;
  foreach (reduction(+:Jtot) reduction(+:V)) {
    V += f[]*dv();
    if (f[] > F_TOL && f[] < 1. - F_TOL) {
      double gx = (f[1] - f[-1])/(2.*Delta);
      double gy = (f[0,1] - f[0,-1])/(2.*Delta);
      Jtot += JFLUX(T[])/rho1*sqrt(sq(gx) + sq(gy))*dv();
    }
  }
  if (V0 == 0.) V0 = V;

  qin = qo - Jtot + (V0 - V)/tauV;
}

/**
## Temperature
$\rho c_p\,\partial_t T = \nabla\!\cdot\!(k\nabla T) + J h_{lv}\delta_s$, with $\delta_s \simeq |\nabla f|$ and phase-weighted
$\rho c_p$, $k$. */
event tracer_diffusion (i++)
{
  scalar thetav[], rr[];
  face vector kappav[];
  foreach() {
    thetav[] = cm[]*(f[]*rho1*cp1 + (1. - f[])*rho2*cp2);
    rr[] = 0.;
    if (f[] > F_TOL && f[] < 1. - F_TOL) {
      double gx = (f[1] - f[-1])/(2.*Delta);
      double gy = (f[0,1] - f[0,-1])/(2.*Delta);
      rr[] = cm[]*JFLUX(T[])*hlv*sqrt (sq(gx) + sq(gy));
    }
  }
  foreach_face() {
    double ff = (f[] + f[-1])/2.;
    kappav.x[] = fm.x[]*(ff*k1 + (1. - ff)*k2);
  }
  diffusion (T, dt, D = kappav, r = rr, theta = thetav);
}

/**
## Gravity and Marangoni forces
$$\mathbf g = g\sin\beta\,\mathbf e_x - g\cos\beta\,\mathbf e_y, 
\qquad \mathbf F_{Ma} =
\left(\nabla\sigma - (\mathbf n\!\cdot\!\nabla\sigma)\,\mathbf n\right)
|\nabla f|,
\qquad \sigma = \sigma_0 - \sigma_T (T - T_{sat}). $$
*/
event acceleration (i++)
{
  scalar sig[];
  foreach()
    sig[] = sigma0 - sigmaT*(T[] - Tsat);

  face vector av = a;
  double gx = gval*sin(beta0);
  double gy = gval*cos(beta0);

  foreach_face(x) {
    double ff = (f[] + f[-1])/2.;
    av.x[] += gx*rho1*ff*alpha.x[]/fm.x[];    // liquid-weighted g_x
  }
  foreach_face(y)
    av.y[] -= gy;                             // uniform g_y

  foreach_face() {
    double fx = (f[] - f[-1])/Delta;
    double fy = (f[0,1] + f[-1,1] - f[0,-1] - f[-1,-1])/(4.*Delta);
    double ds = sqrt (sq(fx) + sq(fy));
    if (ds > 1e-6) {
      double nx = fx/ds, ny = fy/ds;
      double sx = (sig[] - sig[-1])/Delta;
      double sy = (sig[0,1] + sig[-1,1] - sig[0,-1] - sig[-1,-1])/(4.*Delta);
      double sn = sx*nx + sy*ny;
      av.x[] += alpha.x[]/fm.x[]*(sx - sn*nx)*ds;
    }
  }
}

/**
## Logging
$t$, $V$, $\int_\Gamma J/\rho_1\,ds$, $h_{max}$, $T_s$, $\Delta t$. */
event logfile (i += 10)
{
  if (i == 0)
    fprintf (stderr, "# t\tV\tJtot\thmax\tTs\tdt\n");

  double vol = 0., Jtot = 0., hmax = 0., TsN = 0., TsD = 0.;
  foreach (reduction(+:vol) reduction(+:Jtot)
	   reduction(max:hmax) reduction(+:TsN) reduction(+:TsD)) {
    vol += f[]*dv();
    if (f[] > 0.5)
      hmax = max (hmax, y);
    if (f[] > F_TOL && f[] < 1. - F_TOL) {
      double gx = (f[1] - f[-1])/(2.*Delta);
      double gy = (f[0,1] - f[0,-1])/(2.*Delta);
      double gm = sqrt (sq(gx) + sq(gy));
      Jtot += JFLUX(T[])/rho1*gm*dv();
      TsN += T[]*gm*dv();
		 	TsD += gm*dv();
    }
  }
  fprintf (stderr, "%g\t%g\t%g\t%g\t%g\t%g\n", t, vol, Jtot, hmax,
	   TsD > 0. ? TsN/TsD : 0., dt);
}

/** 
In case film becomes too small, we terminate pre-emptively.
*/
event dryout (i += 20)
{
  double fw = 1., band = 3.*L0/(1 << LEVEL);
  foreach (reduction(min:fw))
    if (y < band && f[] < fw)
      fw = f[];
  	  if (fw < 0.5) {
        fprintf (stderr, "# near-dryout at t = %g — stopping\n", t);
    return 1;
  }
}

/** 
AMR at the zones of interest.
*/
#if TREE
event adapt (i++) {
  adapt_wavelet ({f, T, u},
		 (double[]){1e-3, 1e-2, 1e-2, 1e-2}, LEVEL, MINLEVEL);
}
#endif

/**
## Outputs
Interface faces.
*/
event interface_out (t += 0.5)
{
  char name[80];
  sprintf (name, "facets-%.0f", t);
  FILE * fp = fopen (name, "w");
  output_facets (f, fp);
  fclose (fp);
}

/**
### Movies
output_ppm strips of $T$, $f$, $u_x$ and refinement level, and a
view.h composite (temperature + mesh + interface) in `bview.mp4`. */
event movie (t += 0.10)
{
  scalar m[];
  foreach()
    m[] = f[] - 0.5;

  output_ppm (T, file = "T.mp4", n = 1024, box = {{0.,0.},{L0,3.}},
	      min = Tsat, max = Tsat + DT0 + DT1,
	      linear = true, map = cool_warm, mask = m);
  output_ppm (f, file = "f.mp4", n = 1024, box = {{0.,0.},{L0,3.}},
	      min = 0, max = 1, linear = true);
  output_ppm (u.x, file = "ux.mp4", n = 1024, box = {{0.,0.},{L0,3.}},
	      min = -0.5, max = 2., linear = true, map = cool_warm);
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, file = "level.mp4", n = 1024, box = {{0.,0.},{L0,3.}},
	      min = MINLEVEL, max = LEVEL);
}

event viewer (t += 0.1)
{
  view (fov = 17, tx = 0.02, ty = 0.03,   // composite ~[-L0,L0] x [-3.5,3]
	width = 1500, height = 500);      // tune fov/tx/ty by eye once
  clear();

  /* top right: temperature */
  squares ("T", min = Tsat, max = Tsat + DT0 + DT1,
	   linear = true, map = cool_warm);
  draw_vof ("f", lw = 1.5);

  /* bottom right: streamwise velocity */
  translate (y = -3.8) {
    squares ("u.x", min = -0.5, max = 2., linear = true, map = cool_warm);
    draw_vof ("f", lw = 1.5);
  }

  /* top left: volume fraction */
  translate (x = -L0 - 0.3) {
    squares ("f", min = 0, max = 1, linear = true);
    draw_vof ("f", lw = 1.5);
  }

  /* bottom left: mesh */
  translate (x = -L0 - 0.3, y = -3.8) {
    cells (lw = 0.5);
    draw_vof ("f", lw = 1.5);
  }

  /* screen-anchored corner labels + time stamp */
  char s[99];
  sprintf (s, "T    t = %.1f", t);

  draw_string ( s   , 0, size=10, lw=2);
  draw_string ( "f" , 1, size=10, lw=2);
  draw_string ( "T" , 2, size=10, lw=2);
  draw_string ("u.x", 3, size=10, lw=2);

  save ("bview.mp4");
}


event end (t = TEND);

/**
## Plots

~~~gnuplot Total liquid volume
set xlabel 't'
set ylabel 'V'
plot 'log' u 1:2 w l lw 2 t 'V(t)'
~~~

~~~gnuplot Net Phase change
set xlabel 't'
set ylabel 'integral J/rho_1'
plot 'log' u 1:3 w l lw 2 t 'integral J/rho_1'
~~~

~~~gnuplot Maximum film height
set xlabel 't'
set ylabel 'h_{max}'
plot 'log' u 1:4 w l lw 2 t 'h_{max}(t)'
~~~

~~~gnuplot Mean interface temperature
set xlabel 't'
set ylabel 'T_s'
plot 'log' u 1:5 w l lw 2 t 'T_s(t)'
~~~

![Film simulation](heated_film/bview.mp4)
*/
