/**
# Equilibrium tide

This is directly adapted from
[equilibrium_tide.F](https://github.com/myroms/roms/blob/develop/ROMS/Utility/equilibrium_tide.F)
in the [ROMS source code](https://github.com/myroms).

This module computes the equilibrium tide defined as the shape the sea
surface (m) would assume if it were motionless and in equilibrium with
the tide generating forces on a fluid planet. It is used to compute
the tide generation force (TGF) terms for the pressure gradient.
                                                                   
## References

~~~bib
@article{arbic2004,
  title={The accuracy of surface elevations in forward global barotropic and baroclinic tide models},
  author={Arbic, Brian K and Garner, Stephen T and Hallberg, Robert W and Simmons, Harper L},
  journal={Deep Sea Research Part II: Topical Studies in Oceanography},
  volume={51},
  number={25-26},
  pages={3069--3101},
  year={2004},
  publisher={Elsevier}
}

@article{arbic2018,
  title={A primer on global internal tide and internal gravity wave continuum modeling in 
         {HYCOM} and {MITgcm}},
  author={Arbic, Brian K and Alford, Matthew H and Ansong, Joseph K and Buijsman, Maarten C and 
          Ciotti, Robert B and Farrar, J Thomas and Hallberg, Robert W and Henze, Christopher E and 
	  Hill, Christopher N and Luecke, Conrad A and others},
  journal={New frontiers in operational oceanography},
  DOI={10.17125/gov2018.ch13},
  year={2018}
}

@article{doodson1941,
  title={Admiralty manual of tides},
  author={Doodson, Arthur Thomas and Warburg, Harold Dreyer},
  journal={His Majesty's Stationery Office},
  address={London, UK},
  year={1941}
}

@article{egbert2017,
  title={Tidal prediction},
  author={Egbert, Gary D and Ray, Richard D},
  journal={Journal of Marine Research},
  volume={75},
  pages={189-237},
  year={2017}
}

@unpublished{roms2021,
  title={Wiki {ROMS}},
  author={Arango et al.},
  url={https://www.myroms.org/wiki/Tidal_Forcing#Astronomical_Tide_Generating_Forces},
  year={2021},
  month={September}
}
~~~

## Define equilibrium tide constituents structure */

typedef struct {
  double Afl;   // product: amp*f*love
  double amp;   // amplitude (m)
  double chi;   // phase at Greenwich meridian (degrees)
  double f;     // f nodal factor (nondimensional)
  double love;  // tidal Love number factor
  double nu;    // nu nodal factor (degrees)
  double omega; // frequency (1/s)
} Etide_t;

struct {
  Etide_t Q1;   // Q1 component, diurnal
  Etide_t O1;   // O1 component, diurnal
  Etide_t K1;   // K1 component, diurnal
  Etide_t N2;   // N2 component, semi-diurnal
  Etide_t M2;   // M2 component, semi-diurnal
  Etide_t S2;   // S2 component, semi-diurnal
  Etide_t K2;   // N2 component, semi-diurnal  
} Etide = {0};

const double deg2rad = pi/180.;

#include <time.h>

/**
## Computation of the tidal constituents

The tidal constituents are computed relative to the optional starting
date given as argument.

If *Lnodal* is true then lunar nodal correction is applied. */

void equilibrium_tide_constituents (const char * date = "2000-01-01 12:00:00",
				    bool Lnodal = false)
{
  /**
  Compute fundamental astronomical parameters. Adapted from Egbert and
  Ray, 2017, Table 1.

  Set Egbert and Ray time reference for their astronomical parameters
  (Table 1): days since 2000-01-01:12:00:00 */

  struct tm astro_tm;
  assert (strptime ("2000-01-01 12:00:00", "%Y-%m-%d %H:%M:%S", &astro_tm));
  time_t astro_t = mktime (&astro_tm);
  assert (astro_t != (time_t) -1);
  
  /**
  Terrestial time (in centuries) since tide reference date number.
  Recall that the length of a year is 365.2425 days for the now called
  Gregorian Calendar (corrected after 15 October 1582). 

  *We also assume that a day is 86400 seconds (need to check
  this). Note that the original code in ROMS seems to round T to the
  nearest day??, which we do not do here.* */
  
  struct tm date_tm;
  if (!strptime (date, "%Y-%m-%d %H:%M:%S", &date_tm)) {
    fprintf (stderr,
	     "%s:%d: error: date string '%s' is not formatted as '%%Y-%%m-%%d %%H:%%M:%%S'\n",
	     __FILE__, LINENO, date);
    exit (1);
  }
  time_t date_t = mktime (&date_tm);
  assert (date_t != (time_t) -1);
  double T  = difftime (date_t, astro_t)/(86400.*36524.25);
 
  /**
  Compute the harmonic constituents of the equilibrium tide at
  Greenwich. Adapted from Doodson and Warburg, 1941, Table 1.

  Nodal factors "f" and "nu" account for the slow modulation of the
  tidal constituents due (principally) to the 18.6-year lunar nodal
  cycle. */

  if (Lnodal) {

    /**
    Mean longitude of lunar node (Period = 18.6 years). */

    double N = -234.955 - 1934.1363*T; // degrees
    N *= deg2rad; // radians

    Etide.O1.f = 1.009 + 0.187*cos(N) - 0.015*cos(2.0*N);
    Etide.K1.f = 1.006 + 0.115*cos(N) - 0.009*cos(2.0*N);
    Etide.M2.f = 1.0 - 0.037*cos(N);
    Etide.S2.f = 1.0;
    Etide.K2.f = 1.024 + 0.286*cos(N) + 0.008*cos(2.0*N);
    
    Etide.O1.nu = 10.8*sin(N) - 1.3*sin(2.0*N);
    Etide.K1.nu = -8.9*sin(N) + 0.7*sin(2.0*N);
    Etide.M2.nu = -2.1*sin(N);
    Etide.S2.nu = 0.0;
    Etide.K2.nu = -17.7*sin(N) + 0.7*sin(2.0*N);    
  }
  else {
    Etide.O1.f = 1.0;
    Etide.K1.f = 1.0;
    Etide.M2.f = 1.0;
    Etide.S2.f = 1.0;
    Etide.K2.f = 1.0;
    
    Etide.O1.nu = 0.0;
    Etide.K1.nu = 0.0;
    Etide.M2.nu = 0.0;
    Etide.S2.nu = 0.0;
    Etide.K2.nu = 0.0;
  }
  Etide.Q1.f = Etide.O1.f;
  Etide.N2.f = Etide.M2.f;
  Etide.Q1.nu = Etide.O1.nu;
  Etide.N2.nu = Etide.M2.nu;

  /**
  Mean longitude of the moon (Period = tropical month). */

  double s = 218.316 + 481267.8812*T;

  /**
  Mean longitude of the sun (Period = tropical year). */

  double h = 280.466 + 36000.7698*T;

  /**
  Mean longitude of lunar perigee (Period = 8.85 years). */

  double p = 83.353 + 4069.0137*T;

  /**
  Compute tidal constituent phase "chi" (degrees) at Greenwich meridian. */

  Etide.Q1.chi = h - 3.0*s + p - 90.0;
  Etide.O1.chi = h - 2.0*s - 90.0;
  Etide.K1.chi = h + 90.0;
  Etide.N2.chi = 2.0*h - 3.0*s + p;
  Etide.M2.chi = 2.0*h - 2.0*s;
  Etide.S2.chi = 0.0;
  Etide.K2.chi = 2.0*h;

  /**
  Set tidal constituent frequency (1/s). */
  
  Etide.Q1.omega = 0.6495854E-4;
  Etide.O1.omega = 0.6759774E-4;
  Etide.K1.omega = 0.7292117E-4;
  Etide.N2.omega = 1.378797E-4;
  Etide.M2.omega = 1.405189E-4;
  Etide.S2.omega = 1.454441E-4;
  Etide.K2.omega = 1.458423E-4;

  /**
  Set tidal constituent amplitude (m). */

  Etide.Q1.amp =  1.9273E-2;
  Etide.O1.amp = 10.0661E-2;
  Etide.K1.amp = 14.1565E-2;
  Etide.N2.amp =  4.6397E-2;
  Etide.M2.amp = 24.2334E-2;
  Etide.S2.amp = 11.2743E-2;
  Etide.K2.amp =  3.0684E-2;

  /**
  Set tidal constituent Love number factors (1 + k2 - h2).  The Love
  number is defined as the ratio of the body tide to the height of
  the static equilibrium tide. */

  Etide.Q1.love = 0.695;
  Etide.O1.love = 0.695;
  Etide.K1.love = 0.736;
  Etide.N2.love = 0.693;
  Etide.M2.love = 0.693;
  Etide.S2.love = 0.693;
  Etide.K2.love = 0.693;

  /**
  Compute product of amp*f*love. */

  Etide.Q1.Afl = Etide.Q1.amp*Etide.Q1.f*Etide.Q1.love;
  Etide.O1.Afl = Etide.O1.amp*Etide.O1.f*Etide.O1.love;
  Etide.K1.Afl = Etide.K1.amp*Etide.K1.f*Etide.K1.love;
  Etide.N2.Afl = Etide.N2.amp*Etide.N2.f*Etide.N2.love;
  Etide.M2.Afl = Etide.M2.amp*Etide.M2.f*Etide.M2.love;
  Etide.S2.Afl = Etide.S2.amp*Etide.S2.f*Etide.S2.love;
  Etide.K2.Afl = Etide.K2.amp*Etide.K2.f*Etide.K2.love;
}

/**
## Computes surface elevation associated with the equilibrium tide 

Assumes that the coordinate system is longitude in degrees (x) and
latitude in degrees (y). 

Note that the tidal constituents must first be defined by calling
equilibrium_tide_constituents().

The time `t` given (in seconds) is relative to the date/time specified
when calling equilibrium_tide_constituents(). */

void equilibrium_tide (scalar tide, double t)
{
  if (!Etide.Q1.omega) {
    fprintf (stderr,
	     "%s:%d: error: the tidal constituents must first be defined by calling "
	     "equilibrium_tide_constituents()\n", __FILE__, LINENO);
    exit (1);
  }
  
  /**
  Compute astronomical equilibrium tide (m) with diurnal and semidiurnal
  constituents. The long-period constituents and high-order harmonics
  are neglected. */
  
  foreach()
    tide[] =
    sin(2.*deg2rad*y)*(Etide.Q1.Afl*cos(Etide.Q1.omega*t +
					deg2rad*(x + Etide.Q1.chi + Etide.Q1.nu)) + 
		       Etide.O1.Afl*cos(Etide.O1.omega*t +
					deg2rad*(x + Etide.O1.chi + Etide.O1.nu)) + 
		       Etide.K1.Afl*cos(Etide.K1.omega*t +
					deg2rad*(x + Etide.K1.chi + Etide.K1.nu))) + 
    sq(cos(deg2rad*y))*(Etide.N2.Afl*cos(Etide.N2.omega*t +
					 deg2rad*(2.0*x + Etide.N2.chi + Etide.N2.nu)) +  
			Etide.M2.Afl*cos(Etide.M2.omega*t +
					 deg2rad*(2.0*x + Etide.M2.chi + Etide.M2.nu)) +
			Etide.S2.Afl*cos(Etide.S2.omega*t +
					 deg2rad*(2.0*x + Etide.S2.chi + Etide.S2.nu)) +
			Etide.K2.Afl*cos(Etide.K2.omega*t +
					 deg2rad*(2.0*x + Etide.K2.chi + Etide.K2.nu)));
}
