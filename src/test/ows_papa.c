/**
# Ocean Weather Ship Papa test case from GOTM

See [the GOTM web site](https://gotm.net/cases/ows_papa) and
section 12.3.1 of [the GOTM manual](#gotm).

"This scenario is a classical scenario for the Northern Pacific, for
which long term observations of meteorological parameters and
temperature profiles are available. The station Papa at 145$^\circ$W,
50$^\circ$N has the advantage that it is situated in a region where
the horizontal advection of heat and salt is assumed to be
small. Various authors used these data for validating turbulence
closure schemes."

## Results

~~~gnuplot Evolution of the potential temperature (field data)
set term PNG size 800,400
set output 'Tobs.png'
set pm3d map
# jet colormap
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 \
 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392,   \
 0.625 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0,       \
 1 0.498 0 0 )
set cbrange [4:17]
set yrange [-250:0]
unset xlabel
unset key
set ylabel 'z (m)'
set xdata time
set timefmt "%Y/%m/%d"
set xrange ["1961/03/25" : "1962/03/25"]
set xtics rotate by 45
set xtics center offset 0,-0.5
set timefmt "%Y/%m/%d %H:%M:%S"
splot "< awk '{ if (NF == 4) { time = $1 \" \" $2; print \"\"; } \
                else if (NF == 2) print time, $0;}' tprof.dat" u 1:3:4
~~~

~~~gnuplot Evolution of the potential temperature (k-$\epsilon$ model)
set output 'T.png'
set timefmt "%s"
start = system("date -d'1961-03-25' +%s")
splot '../ows_papa/out' u (start + $1):2:5
~~~

Something strange seems to happen for KPP around 01/01/1962.

~~~gnuplot Evolution of the potential temperature (KPP model)
set output 'Tkpp.png'
splot '../ows_papa-kpp/out' u (start + $1):2:5
~~~

## References

~~~bib
@manual{gotm, 
  title = {{GOTM} Source Code and Test Case Documentation}, 
  author = {Lars Umlauf and Hans Burchard and Karsten Bolding}, 
  edition = {version 4},
  year = {2018}, 
  pdf = {https://gotm.net/manual/stable/pdf/a4.pdf}
}
~~~
*/

#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/gotm.h"

/**
We add temperature and salinity fields for each layer. */

event defaults (i = 0)
{
  T = new scalar[nl];
  S = new scalar[nl];
}

event cleanup (i = end) {
  delete ({T, S});
}

/**
This case depends on (0D) Coriolis terms. */

#define F0() (2.*(2.*pi/86164.)*sin(2.*pi*y/360.))

event coriolis (i++)
{
  foreach() {
    double cosomega = cos (F0()*dt);
    double sinomega = sin (F0()*dt);
    foreach_layer() {
      double ua = u.x[];
      u.x[] =   ua*cosomega + u.y[]*sinomega;
      u.y[] = - ua*sinomega + u.y[]*cosomega;
    }
  }
}

/**
The external GOTM global variables `a`, `g1` and `g2` control the
Jerlov-type extinction of the short-wave radiation flux. */

#include <gotm/observations/observations.h> // for a, g1, g2

/**
The starting date. */

int syear = 1961, smonth = 3, sday = 25;

int main()
{
  foreach_dimension()
    periodic (right);
  size (100e3);
  origin (- 145. - L0/2., 50. - L0/2.);
  G = 9.81;
  N = 1;
  nl = 250;
  CFL = CFL_H = HUGE;
  
  /**
  These are set to match those used in the GOTM reference case. */

  observations_a = 0.58, observations_g1 = 0.35, observations_g2 = 23.;
    
  //  eqstate_eq_state_mode = 1; // does not work with kpp

  /**
  ## The GOTM turbulence model setup

  Two versions can be used KPP or k-$\epsilon$. The number of
  parameters for k-$\epsilon$ seems way too large... */

  meanflow_maxitz0b = 1;	// default = 10

#if KPP // kpp
  turbulence_turb_method = 99;	// default = 0
  kpp_kpp_sbl = true;	        // default = 0
  kpp_kpp_bbl = true;	        // default = 0
  kpp_kpp_interior = true;      // default = 0
  kpp_ric = 0.3;      	        // default = 0
#else // k-epsilon
  turbulence_iw_model = 2;	// default = 0 // no effect whatsoever!!
  
  turbulence_turb_method = 3;	// default = 2
  turbulence_compute_c3 = 1;	// default = 0
  turbulence_length_lim = 1;	// default = 0
  turbulence_const_num = 0.0001;	// default = 0.0005
  turbulence_const_nuh = 0.0001;	// default = 0.0005
  turbulence_k_min = 1e-06;	// default = 1e-08
  turbulence_kb_min = 1e-10;	// default = 1e-08
  turbulence_epsb_min = 1e-14;	// default = 1e-12
  turbulence_gen_m = 1;	// default = 1.5
  turbulence_gen_n = -0.67;	// default = -1
  turbulence_cpsi1 = 1;	// default = 1.44
  turbulence_cpsi2 = 1.22;	// default = 1.92
  turbulence_cpsi3minus = 0.05;	// default = 0
  turbulence_sig_kpsi = 0.8;	// default = 1
  turbulence_sig_psi = 1.07;	// default = 1.3
  turbulence_ce3minus = -0.4;	// default = 0
  turbulence_my_length = 3;	// default = 1

  turbulence_scnd_method = 1;	// default = 0
  turbulence_kb_method = 1;	// default = 0
  turbulence_epsb_method = 1;	// default = 0
  turbulence_scnd_coeff = 7;	// default = 0
  turbulence_cc1 = 3.6;	// default = 0
  turbulence_cc2 = 0.8;	// default = 0
  turbulence_cc3 = 1.2;	// default = 0
  turbulence_cc4 = 1.2;	// default = 0
  turbulence_cc6 = 0.3;	// default = 0
  turbulence_ct1 = 3.28;	// default = 0
  turbulence_ct2 = 0.4;	// default = 0
  turbulence_ct3 = 0.4;	// default = 0
  turbulence_ct5 = 0.4;	// default = 0
  turbulence_ctt = 0.8;	// default = 0
  turbulence_alpha = 0.7;	// default = 0
  turbulence_nuhiw = 1e-05;	// default = 5e-05
#endif

  /**
  ## Forcing flux data

  These were derived from observations. */
  
  // This seems to be unreliable
  // #define GOTM_CASES "https://github.com/gotm-model/cases/raw/master/"
  // so we use a mirror here
  
  #define GOTM_CASES "http://basilisk.fr/src/test/"

  system ("test -f sst.dat || ( "
	  "wget -q " GOTM_CASES "ows_papa/heatflux.dat && "
	  "wget -q " GOTM_CASES "ows_papa/momentumflux.dat && "
	  "wget -q " GOTM_CASES "ows_papa/sprof.dat && "
	  "wget -q " GOTM_CASES "ows_papa/sst.dat && "
	  "wget -q " GOTM_CASES "ows_papa/swr.dat && "
	  "wget -q " GOTM_CASES "ows_papa/tprof.dat "
	  ")");
  
  run();
}

/**
## Initial conditions 

A rough function to read GOTM-formatted data files. */
 
void init_profile (Point point, const char * fname,
		   int sy, int sm, int sd, scalar s)
{
  FILE * fp = fopen (fname, "r");
  if (!fp) {
    perror (fname);
    exit (1);
  }
  int yy, m, d, n;
  do {
    if (fscanf (fp, "%d/%d/%d %*s %d %*d", &yy, &m, &d, &n) != 4) {
      fprintf (stderr, "%s: error: cannot read date\n", fname);
      exit (1);
    }
    if (yy < sy || (yy == sy && m < sm) || (yy == sy && m == sm && d < sd))
      for (int i = 0; i < n; i++)
	fscanf (fp, "%*f %*f");
    else
      break;
  } while (1);
  double zd[n], od[n];
  int i = 0;
  while (fscanf (fp, "%lf %lf", &zd[i], &od[i]) == 2)
    i++;
  if (i != n) {
    fprintf (stderr, "%s: error: not enough points: read %d out of %d\n",
	     fname, i, n);
    exit (1);
  }
  double zc = zb[];
  foreach_layer() {
    zc += h[]/2.;
    int i;
    for (i = 1; i < n && zd[i] >= zc; i++);
    if (i < n)
      s[] = (od[i-1]*(zd[i] - zc) + od[i]*(zc - zd[i-1]))/(zd[i] - zd[i-1]);
    else
      s[] = od[i-1];
    zc += h[]/2.;
  }
  fclose (fp);
}

/**
The initial temperature and salinity profiles are obtained from the
data. */
 
event init (i = 0)
{
  if (turbulence_turb_method != 99)
    turbulence_report_model();
  foreach() {
    zb[] = - 250.;
    foreach_layer()
      h[] = 250./nl;
    init_profile (point, "sprof.dat", syear, smonth, sday, S);
    init_profile (point, "tprof.dat", syear, smonth, sday, T);
  }
}

/**
## Timestep

Set to one hour, as in the GOTM test case. */
 
event timestep (t += 3600);

/**
## Surface fluxes

### Momentum */

vector tau[];

event input_momentum_flux (i += 3)
{
  static FILE * fp = fopen ("momentumflux.dat", "r");
  static coord hfp, hfn;    
  // the file is sampled at 3 hours intervals
  if (i == 0) {
    int y, m, d;
    while (fscanf (fp, "%d/%d/%d %*s %lf %lf", &y, &m, &d, &hfn.x, &hfn.y) == 5
	   && (y < syear || m < smonth || d < sday));
  }
  if (i % 3 == 0) {
    hfp = hfn;
    fscanf (fp, "%*s %*s %lf %lf", &hfn.x, &hfn.y);
  }
  coord hf;
  foreach_dimension()
    hf.x = hfp.x + (i % 3)/3.*(hfn.x - hfp.x);
  foreach()
    foreach_dimension()
      tau.x[] = hf.x;
  airsea_tau = tau;
}

/**
### Heat */

scalar heat_flux[];

event input_heat_flux (i += 1)
{
  static FILE * fp = fopen ("heatflux.dat", "r");
  static double hfp, hfn;    
  // the file is sampled at 3 hours intervals
  if (i == 0) {
    int y, m, d;
    while (fscanf (fp, "%d/%d/%d %*s %lf", &y, &m, &d, &hfn) == 4 &&
	   (y < syear || m < smonth || d < sday));
  }
  if (i % 3 == 0) {
    hfp = hfn;
    fscanf (fp, "%*s %*s %lf", &hfn);
  }
  double hf = hfp + (i % 3)/3.*(hfn - hfp);
  foreach()
    heat_flux[] = hf;
  airsea_heat_flux = heat_flux;
}

/**
### Short-wave radiation

We keep the option to compute these from solar zenith angle etc... but
it is not used by default. */

#include <gotm/airsea/solar_zenith_angle.h>
#include <gotm/airsea/albedo_water.h>
#include <gotm/util/time.h>

scalar swr_flux[];

event input_swr_flux (i += 1)
{
  static FILE * fp = fopen ("swr.dat", "r");
  static double hfp, hfn;    
  // the file is sampled at 3 hours intervals
  if (i == 0) {
    int y, m, d;
    while (fscanf (fp, "%d/%d/%d %*s %lf", &y, &m, &d, &hfn) == 4 &&
	   (y < syear || m < smonth || d < sday));
  }
  if (i % 3 == 0) {
    hfp = hfn;
    fscanf (fp, "%*s %*s %lf", &hfn);
  }
  double hf = hfp + (i % 3)/3.*(hfn - hfp);

#if 0  
  int julianday, year = 1959, month = 9, day = 14;
  time_julian_day (&year, &month, &day, &julianday);
  julianday += t/86400.;
  // fixme: assumes t start at midnight
  double secondsofday = fmod (t, 86400.);
  double hoursofday = secondsofday/3600.;
  
  int firstjan;
  month = 1, day = 1;
  time_julian_day (&year, &month, &day, &firstjan);
  int yearday = julianday - firstjan + 1;
#endif
  
  foreach() {
#if 0    
    realtype zenith_angle =
      airsea_solar_zenith_angle (&yearday, &hoursofday, &x, &y);
    realtype albedo = airsea_albedo_payne (&zenith_angle);
#endif
    swr_flux[] = hf; // *(1. - albedo);
  }
  airsea_swr_flux = swr_flux;
}

/**
## Outputs */

void profile (FILE * fp)
{
  foreach() {
    double z = zb[];
    foreach_layer()
      fprintf (fp, "%g %g %g %g %g %g %g %g\n", t, z + h[]/2.,
	       u.x[], u.y[], T[], S[],
	       meanflow_nn.a[point.l + 1], turbulence_nuh.a[point.l + 1]),
      z += h[];
  }
  fprintf (fp, "\n");  
}

event profiles (i += 12)
  profile (stdout);

event end (i = 365*24)
  profile (stderr);
