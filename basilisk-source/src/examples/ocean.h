/**
# A generic ocean model

This generic header is used to setup both the [Gulf
Stream](gulf-stream.c) and [Global Oceanic Circulation](global.c)
examples.

We use the hydrostatic solver with time-implicit integration of the
barotropic mode, on a regular Cartesian grid and in spherical
coordinates. */

#include "grid/multigrid.h"
#include "spherical.h"
#include "layered/hydro.h"
#include "layered/implicit.h"
#include "profiling.h"

/**
## Coriolis acceleration and bottom friction

The Coriolis acceleration takes its standard definition. */

const double Omega = 7.292205e-5;
#define F0() (2.*Omega*sin(y*pi/180.))

/**
The quadratic bottom friction coefficient is set to 2 x 10^-3^ (see
Table 2 in [Hurlburt & Hogan, 2000](#hurlburt2000)). The friction is
only applied in the deepest wet layer.

We have three options for bottom friction. The first one, which we
use, weighs the friction coefficient with the cube of `zb/zbs`
i.e. the "compressed" bathymetry over the smoothed bathymetry. The
main effect of this weighing is an increase of the friction
coefficient (by a factor up to approx. seven) in shallow seas. The
rationale for this weighing is discussed in [Metzger & Hurlburt,
1996](#metzger1996) Appendix A. This has a non-negligeable impact on
Gulf Stream separation (but a small impact on the overall
circulation). */

const double Cb = 2e-3;
scalar zbs[];
#if 1
#define K0() (point.l > 0 && h[0,0,-1] > dry ? 0. : h[] < dry ? HUGE :	\
	      Cb*cube(zb[]/zbs[])*norm(u)/h[])
#elif 1
#define K0() (y > 50 ? (y - 50)/3600. :					\
	      y < 9.5 ? 1./3600. :					\
	      point.l > 0 && h[0,0,-1] > 10. ? 0. : h[] < dry ? HUGE :	\
	      Cb*norm(u)/h[])
#else
#define K0() (point.l > 0 && h[0,0,-1] > 10. ? 0. : h[] < dry ? HUGE : Cb*norm(u)/h[])
#endif
#include "layered/coriolis.h"

/**
## Isopycnal layers

The isopycnal layers thicknesses and relativity density differences
need to be defined by the calling program. */

#include "layered/isopycnal.h"
#define rho0 1000.
extern double * dh, * drho;

/**
## Diapycnal entrainment

See [entrainment.h](/src/layered/entrainment.h) for explanations. The
mininum and maximum layer thicknesses, the average entrainment
velocity and coefficient of additional interfacial friction associated
with entrainment, must be specified by the calling program. */

#include "layered/entrainment.h"
extern double * hmin, * hmax, omr, Cm;

/**
The volume of each layer is conserved, using (vertical) fluxes between
layers, if necessary (for example when the AMOC is included through
[boundary fluxes](bflux.h) for the [Gulf Stream](gulf-stream.c)
model). */

#include "layered/conservation.h"

/**
## Various utilities */

#include "terrain.h" // for bathymetry
#include "input.h"   // for wind inputs
#include "layered/perfs.h"

/**
Time averages. */

const double hour = 3600., day = 86400., month = 30.*day, year = 365.*day;

vector ua, ud;
scalar Ha, etam[], eta2[];
vector uga[], ugd[];

/**
The latitude is bounded my $\pm$ maxlat. */

double maxlat = 90;

/**
The default starting time for averaging (i.e. spinup time) is set to
5 years. This is a minimum. */

double tspinup = 5.*year;

/**
## Initial conditions 

The function below applies some Laplacian smoothing to the real
bathymetry, as done in [Hulburt & Hogan, 2000 &
2008](#hurlburt2000). Note that the runs are robust even without this
smoothing, which may or may not be necessary to obtain realistic
results. We keep it for consistency with [Hulburt & Hogan,
2000](#hurlburt2000). */

void laplacian_smoothing()
{
  for (int i = 0; i < 2; i++) {
    foreach() {
      if (zb[] < 0.)
	zbs[] = (zb[1] + zb[-1] + zb[0,1] + zb[0,-1] +
		 zb[1,1] + zb[-1,-1] + zb[-1,1] + zb[1,-1])/8.;
      else
	zbs[] = zb[];
    }
    foreach()
      zb[] = zbs[];
  }  
}

/**
We have the option to restart (from a previous "dump" file, see below)
or start from initial conditions (i.e. a "flat" ocean at rest). */

event init (i = 0)
{
  if (restore ("restart"))
    event ("metric");
  else {

    /**
    The terrain uses the ETOPO2 bathymetric KDT database, which needs
    to be generated first. See the [*xyz2kdt*
    manual](http://gerris.dalembert.upmc.fr/xyz2kdt.html) for
    instructions. */
    
    terrain (zb, "~/terrain/etopo2", NULL);
    laplacian_smoothing();

    /**
    We have the option to use the real bathymetry (with a "coastline"
    at - 10 meters) or the "compressed bathymetry" described in Note c
    for Table 1 in [Hulburt & Hogan, 2000](#hurlburt2000). Note that
    the model used in H&H, 2000 cannot deal with isopycnals
    intersecting the bathymetry, which is the main reason for the
    "topography compression". This is not the case with Basilisk (see
    [/src/test/bleck.c]()) and the present code runs fine with the
    real bathymetry, however the results are less realistic than with
    the compressed bathymetry, probably due to a tuning of the
    isopycnal layers, boundary fluxes etc. which is specific to this
    bathymetry. This should be investigated further. */
    
#if 0 // !COMPRESSED
    foreach()
      if (zb[] > - 10)
	zb[] = 100.;
#else // COMPRESSED
    foreach() {
      double zbmin = - 6500.;
      if (zb[] > - 200 || fabs(y) > maxlat)
	zb[] = 1000.;
      else
	zb[] = zbmin + 0.82*(zb[] - zbmin);
    }
#endif // COMPRESSED

    /**
    This initializes the isopycnal layers, based on their nominal thicknesses `dh`. */
    
    foreach() {
      double z = 0.;
      for (point.l = nl - 1; point.l >= 0; point.l--) {
	if (point.l > 0 && z - dh[point.l] > zb[])
	  h[] = dh[point.l];
	else
	  h[] = max(z - zb[], 0.);
	z -= h[];
      }
    }

    /**
    We reset the fields used to store various averages/diagnostics. */
    
    reset ({etam, eta2, ua, ud, Ha, uga, ugd}, 0.);
  }

  /**
  ## Boundary conditions 

  We set a dry, high terrain on all the "wet" domain
  boundaries. Without these the Coriolis acceleration seems to be
  "less balanced" on these boundaries. */
  
  foreach_dimension()
    if (!Period.x) {
      u.t[right] = dirichlet(0);
      u.t[left] = dirichlet(0);
      zb[right] = 1000;
      zb[left] = 1000;
      h[right] = 0;
      h[left] = 0;
    }
}

/**
## Wind stress

The surface wind stress is modelled as in [H & H,
2000](#hurlburt2000), page 293. Note that `Cw` and `rho_air` are only
used with the "COADS" wind climatology since the Hellerman &
Rosenstein climatology directly gives the wind stress.*/

double Cw = 1.5e-3, rho_air = 1.2;

/**
This function loads the [Hellerman & Rosenstein, 1983](#hellerman1983)
(default) or the COADS wind climatology. */

void load_wind (vector wind, int index)
{
  char name[80];
#if COADS  
  sprintf (name, "coads-%d_5.asc", index + 1);
  input_grd (wind.x, file = name, linear = true, periodic = {true, false}, nodatavalue = 0.);
  sprintf (name, "coads-%d_6.asc", index + 1);
  input_grd (wind.y, file = name, linear = true, periodic = {true, false}, nodatavalue = 0.);
#else // HR
  sprintf (name, "wind/hr-%d-x.asc", index + 1);
  input_grd (wind.x, file = name, linear = true, periodic = {true, false}, nodatavalue = 0.,
	     smooth = 1);
  sprintf (name, "wind/hr-%d-y.asc", index + 1);
  input_grd (wind.y, file = name, linear = true, periodic = {true, false}, nodatavalue = 0.,
	     smooth = 1);
#endif // HR
}

/**
At initialisation, we check whether the wind climatology files already
exist, if they don't we get them from their websites and convert them
to the [GRD ASCII raster format](/src/input.h#input_grd). This
requires the [GDAL](https://gdal.org) conversion utilities, which can
easily be installed on Debian systems using

~~~bash
sudo apt install gdal-bin
~~~

Alternatively, you can directly retrieve the preprocessed files for
the HR climatology with something like

~~~bash
wget http://basilisk.fr/src/examples/gulf-stream/wind.tgz
tar xzvf wind.tgz
~~~
*/

event init (i = 0)
{
#if COADS
  system ("if ! test -f coads-1_1.asc; then "
	  "  wget https://github.com/NOAA-PMEL/FerretDatasets/raw/master/data/coads_climatology.cdf"
	  "   -O coads_climatology.cdf; "
	  "  for i in `seq 1 1 12`; do "
	  "    gdal_translate -of AAIGrid -ot float32 -b $i -sds -q "
	  "    coads_climatology.cdf coads-$i.asc; "
	  "  done "
	  "fi "
	  );
#else // HR
  system ("if ! test -f wind/hr-1-x.asc; then "
	  " mkdir wind; cd wind; "
	  " wget https://iridl.ldeo.columbia.edu/SOURCES/.HELLERMAN/.taux/data.nc -O data.nc; "
	  "  for i in `seq 1 1 12`; do "
	  "    gdal_translate -of AAIGrid -ot float32 -b $i -sds -q "
	  "    data.nc hr-$i-x.asc; "
	  "  done; "
	  " wget https://iridl.ldeo.columbia.edu/SOURCES/.HELLERMAN/.tauy/data.nc -O data.nc; "
	  "  for i in `seq 1 1 12`; do "
	  "    gdal_translate -of AAIGrid -ot float32 -b $i -sds -q "
	  "    data.nc hr-$i-y.asc; "
	  "  done; "
	  "fi "
	  );
#endif // HR
}

/**
We use two vector fields, one before and one after the current
time. */

vector wind1[], wind2[];

event acceleration (i++)
{
  int i = t/month;
  double deltaw = (t - i*month - month/2.)/month;
  while (i > 11) i -= 12;
  int i1 = deltaw > 0 ? i : i - 1;
  int i2 = deltaw > 0 ? i + 1: i;
  if (deltaw < 0.) deltaw += 1.;
  if (i1 < 0) i1 = 11;
  if (i2 > 11) i2 = 0;
  static int t1 = -1, t2 = -1;
  if (i1 != t1)
    load_wind (wind1, i1), t1 = i1;
  if (i2 != t2)
    load_wind (wind2, i2), t2 = i2;

  /**
  The wind stress is added directly as an acceleration, only in the
  topmost layer and only if the fluid layer thickness is larger than
  10 metres. We also interpolate linearly in time, between the times
  associated with `wind1` and `wind2`. */

  foreach_face() {
    point.l = nl - 1;
    if (hf.x[] > 10.) {
      double tauw = ((1. - deltaw)*(wind1.x[] + wind1.x[-1]) +
		     deltaw*(wind2.x[] + wind2.x[-1]))/2.;
#if COADS
      double n = Cw*rho_air*sqrt(sq(tauw.x) + sq(tauw.y));
#else // HR
      double n = 0.1 [0]; // conversion from dynes/cm^2 to kg/m/s^2
#endif
      ha.x[] += n*tauw/rho0;
    }
  }
}

/**
## Horizontal viscosity

We add a (small) Laplacian horizontal viscosity in each layer. It is
not clear whether this is really necessary i.e. how sensitive the
results are to this parameter. At low resolutions, horizontal
viscosity is most probably dominated by numerical diffusion due to
upwinding in the advection scheme. The simulation runs fine without
viscosity for a resolution of 1024 x 512 (we haven't tested other
resolutions) and the statistics (and dynamics) are undistinguishable
from those with viscosity. */

double nu_H = 10; // m^2/s

event viscous_term (i++)
{
  if (nu_H > 0.) {
    vector d2u[];
    foreach_layer() {
      double dry = 1.;
      foreach()
	foreach_dimension()
	d2u.x[] = 2.*(sq(fm.x[1])/(cm[1] + cm[])*u.x[1]*(h[1] > dry) +
		      sq(fm.x[])/(cm[-1] + cm[])*u.x[-1]*(h[-1] > dry) +
		      sq(fm.y[0,1])/(cm[0,1] + cm[])*u.x[0,1]*(h[0,1] > dry) +
		      sq(fm.y[0,-1])/(cm[0,-1] + cm[])*u.x[0,-1]*(h[0,-1] > dry))
	/(sq(Delta)*cm[]);
      foreach()
	foreach_dimension() {
	double n = 2.*(sq(fm.x[1])/(cm[1] + cm[])*(1. + (h[1] <= dry)) +
		       sq(fm.x[])/(cm[-1] + cm[])*(1. + (h[-1] <= dry)) +
		       sq(fm.y[0,1])/(cm[0,1] + cm[])*(1. + (h[0,1] <= dry)) +
		       sq(fm.y[0,-1])/(cm[0,-1] + cm[])*(1. + (h[0,-1] <= dry)))
	  /(sq(Delta)*cm[]);
	u.x[] = (u.x[] + dt*nu_H*d2u.x[])/(1. + dt*nu_H*n);
      }
    }
  }
}

/**
## Daily outputs 

We compute the kinetic energy in the top and bottom layer. */

event outputs (t += day)
{
  double ke = 0., keb = 0., vol = 0., volb = 0.;
  scalar etad[], m[], nu[];
  
  foreach(reduction(+:ke) reduction(+:vol) reduction(+:keb) reduction(+:volb)) {
    point.l = 0;
    keb += dv()*h[]*(sq(u.x[]) + sq(u.y[]));
    volb += dv()*h[];
    foreach_layer() {
      ke += dv()*h[]*(sq(u.x[]) + sq(u.y[]));
      vol += dv()*h[];
    }
    point.l = nl - 1;
    etad[] = h[] > dry ? eta[] : 0.;
    nu[] = h[] > dry ? norm(u) : 0.;
    m[] = etad[] - zbs[];
  }

  /**
  Various diagnostics. */
  
  if (i == 0) {
    fprintf (stderr, "t ke/vol keb/vol dt "
	     "mgH.i mgH.nrelax etad.stddev nu.stddev");
    for (int l = 0; l < nl; l++)
      fprintf (stderr, " d%s%d.sum/dt", h.name, l);
    fputc ('\n', stderr);
  }
  else
    fprintf (stderr, "%g %g %g %g %d %d %g %g", t/day, ke/vol/2., keb/volb/2., dt,
	     mgH.i, mgH.nrelax,
	     statsf (etad).stddev, statsf(nu).stddev);

  /**
  This computes the variation of the volume-averaged thickness of each
  layer. This is zero when the volume of each layer is conserved. */

  static double s0[NL] = {0}, t0 = 0.;
  foreach_layer() {
    double s = statsf(h).sum;
    if (t0 == 0.)
      fprintf (stderr, " 0");
    else
      fprintf (stderr, " %g", (s - s0[_layer])/(t - t0));
    s0[_layer] = s;
  }
  fputc ('\n', stderr);
  t0 = t;

  /**
  Animations of the free-surface height, norm of the velocity and
  surface vorticity. */
  
  #define BOX {{X0, max(Y0, - maxlat)}, {X0 + L0, min(Y0 + L0/dimensions().x, maxlat)}}
  output_ppm (etad, mask = m, file = "eta.mp4", n = clamp(N,1024,2048),
	      min = -0.8, max = 0.6, box = BOX, map = jet);
  output_ppm (nu, mask = m, file = "nu.mp4", n = clamp(N,1024,2048),
	      min = 0, max = 1.5, box = BOX, map = cool_warm);

  char name[80];
  sprintf (name, "u%d", nl - 1);
  vector utop = lookup_vector (name);
  if (utop.x.i >= 0) {
    vorticity (utop, nu);
    output_ppm (nu, mask = m, file = "omega.mp4", n = clamp(N,1024,2048),
		linear = false,
		//	      spread = 5,
		min = -0.5e-4, max = 0.5e-4, box = BOX, map = blue_white_red);
  }
}

/**
## Real-time display on GPUs */

#if _GPU && SHOW
event display (i++)
{
  scalar etad[], m[], nu[];
  foreach() {
    etad[] = h[] > dry ? eta[] : 0.;
    m[] = etad[] - zbs[];
    nu[] = h[] > dry ? norm(u) : 0.;
  }
  output_ppm (etad, mask = m, fp = NULL, fps = 30,
	      n = clamp(N,1024,2048),
	      min = -0.8, max = 0.6, box = BOX, map = jet);
  output_ppm (nu, mask = m, fp = NULL, fps = 30,
	      n = clamp(N,1024,2048),
	      min = 0, max = 1.5, box = BOX, map = cool_warm);
}
#endif // _GPU && SHOW

/**
## Fluxes through vertical cross-sections 

These diagnostics are based on Figure 2a and Table 6 of [Hurlburt & Hogan,
2000](#hurlburt2000). */

Flux fluxes[] = {
  { "florida", {{- 80.25, 27.}, {- 78.75, 27.}},
    "Florida Straits at 27N" },
  { "abaco",   {{- 77.2, 26.5}, {- 74.13, 26.5}},
    "East of Abaco Island at 26.5N" },
  { "hatteras",   {{- 76.15, 34.25}, {- 74.5, 34.25}},
    "Gulf Stream at Cape Hatteras" },
  { "38N",   {{- 74.2, 38}, {- 72.8, 38}},
    "Western boundary current at 38N" },
  { "NACwest",   {{- 52.5, 44}, {- 53.9, 43}},
    "N. Atlantic Current, west of Grand Banks" },
  { "NACeast",   {{- 49, 44}, {- 46, 44}},
    "N. Atlantic Current, east of Grand Banks at 44N" },
  { "south",   {{- 50, 42.8}, {- 50, 36}},
    "S. of Grand Banks to 36N" },
  { "68W",   {{- 68, 38.34}, {- 68, 33.48}},
    "Gulf Stream at 68W" },
  {NULL}
};

event fluxes1 (t += day)
  output_fluxes (fluxes, h, u);

/**
## Monthly snapshots 

We dump all fields (and the vorticity). This can be used for
post-processing and/or to restart the simulation. */

event snapshots (t += month)
{
  scalar omega[];
  char name[80];
  sprintf (name, "u%d", nl - 1);
  vector utop = lookup_vector (name);
  if (utop.x.i >= 0)
    vorticity (utop, omega);
  else
    reset ({omega}, 0);
  dump (zero = false);
}

/**
## Time averages 

We allocate new fields to store the time-averaged velocities,
geostrophic velocities, free-surface and their standard deviations. */

event init (i = 0)
{
  ua = new vector[nl];
  ud = new vector[nl];
  Ha = new scalar[nl];
}

event average (t = tspinup; t <= 60.*year; i++)
{
  double Dt = t - tspinup;
  foreach() {
    foreach_layer() {
      Ha[] = (Dt*Ha[] + dt*h[])/(Dt + dt);
      foreach_dimension() {
        ua.x[] = (Dt*ua.x[] + dt*u.x[])/(Dt + dt);
	ud.x[] = (Dt*ud.x[] + dt*sq(u.x[]))/(Dt + dt);
      }
    }

    coord ug = geostrophic_velocity (point);
    foreach_dimension() {
      uga.x[] = (Dt*uga.x[] + dt*ug.x)/(Dt + dt);
      ugd.x[] = (Dt*ugd.x[] + dt*sq(ug.x))/(Dt + dt);      
    }
      
    etam[] = (Dt*etam[] + dt*eta[])/(Dt + dt);
    eta2[] = (Dt*eta2[] + dt*sq(eta[]))/(Dt + dt);
  }
}

/**
We make movies of the averaged free-surface and standard deviation. */

event average_outputs (t = tspinup; t += 30*day)
{
  scalar etad[], m[];
  foreach() {
    point.l = nl - 1;
    etad[] = eta2[] - sq(etam[]) > 0. ? sqrt(eta2[] - sq(etam[])) : 0.;
    m[] = (h[] > dry ? eta[] : 0.) - zbs[];
  }
  output_ppm (etad, mask = m, file = "etad.mp4", n = clamp(N,1024,2048),
	      linear = false,
	      min = 0, max = 0.4, box = BOX, map = jet);
  output_ppm (etam, mask = m, file = "etam.mp4", n = clamp(N,1024,2048),
	      linear = false,
	      min = -0.6, max = 0.6, box = BOX, map = jet);
}
