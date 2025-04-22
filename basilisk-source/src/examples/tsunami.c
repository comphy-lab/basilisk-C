/**
# The 2004 Indian Ocean tsunami

The 2004 Indian Ocean tsunami was caused by a large-scale fault
rupture (> 1000 km) at the Indian–Australian and Eurasian–Andaman
plate boundaries. This example uses the fault model of [Grilli et
al, 2007](/src/references.bib#grilli2007) as initial conditions for a 
Saint-Venant solution of the subsequent tsunami. A similar setup is
discussed in [Popinet, 2011](/src/references.bib#popinet2011).

## Solver setup

The following headers specify that we use spherical coordinates and
the [Saint-Venant solver](/src/saint-venant.h) together with
[(dynamic) terrain reconstruction](/src/terrain.h) and the [Okada
fault model](/src/okada.h). */

#include "spherical.h"
#include "saint-venant.h"
#include "terrain.h"
#include "okada.h"

/**
We then define a few useful macros and constants. */

int maxlevel = 10;
#define MINLEVEL 5
double ETAE  = 1e-2; // error on free surface elevation (1 cm)
double HMAXE = 5e-2; // error on maximum free surface elevation (5 cm) 

/**
The maximum number of levels to use can be set as an argument to the
program. */

int main (int argc, char * argv[])
{
  if (argc > 1)
    maxlevel = atoi(argv[1]);

  /**
  Here we setup the domain geometry. We choose to use metre as length
  unit, so we set the radius of the Earth (required for the [spherical
  coordinates](/src/spherical.h)) in metres. The *x* and *y*
  coordinates are longitude and latitude in degrees, so we set the
  size of the box *L0* and the coordinates of the lower-left corner
  *(X0,Y0)* in degrees. */

  Radius = 6371220. [1];
  // the domain is 54 degrees squared
  size (54.);
  // centered on 94,8 longitude,latitude
  origin (94. - L0/2., 8. - L0/2.);

  /**
  *G* is the acceleration of gravity required by the Saint-Venant
  solver. This is the only dimensional parameter. We rescale it so that
  time is in minutes. */

  // acceleration of gravity in m/min^2
  G = 9.81*sq(60.);

  /**
  When using a tree (i.e. adaptive) discretisation, we want to start
  with the coarsest grid, otherwise we directly refine to the maximum
  level. Note that *1 << n* is C for $2^n$. */

#if TREE
  // 32^2 grid points to start with
  init_grid (1 << MINLEVEL);
#else // Cartesian
  // 1024^2 grid points
  init_grid (1 << maxlevel);
#endif

  /**
  We then call the *run()* method of the Saint-Venant solver to
  perform the integration. */

  run();
}

/**
We declare and allocate another scalar field which will be used to
store the maximum wave elevation reached over time. */

scalar hmax[];

/**
## Boundary conditions

We set the normal velocity component on the left, right and bottom
boundaries to a "radiation condition" with a reference sealevel of
zero. The top boundary is always "dry" in this example so can be left
alone. Note that the sign is important and needs to reflect the
orientation of the boundary. */

u.n[left]   = - radiation(0);
u.n[right]  = + radiation(0);
u.n[bottom] = - radiation(0);

/**
## Adaptation

Here we define an auxilliary function which we will use several times
in what follows. Again we have two *#if...#else* branches selecting
whether the simulation is being run on an (adaptive) tree or a
(static) Cartesian grid.

We want to adapt according to two criteria: an estimate of the error
on the free surface position -- to track the wave in time -- and an
estimate of the error on the maximum wave height *hmax* -- to make
sure that the final maximum wave height field is properly resolved.

We first define a temporary field (in the
[automatic variable](http://en.wikipedia.org/wiki/Automatic_variable)
*η*) which we set to $h+z_b$ but only for "wet" cells. If we used
$h+z_b$ everywhere (i.e. the default $\eta$ provided by the
Saint-Venant solver) we would also refine the dry topography, which is
not useful. */

int adapt() {
#if TREE
  scalar eta[];
  foreach()
    eta[] = h[] > dry ? h[] + zb[] : 0;

  /**
  We can now use wavelet adaptation on the list of scalars *{η,hmax}*
  with thresholds *{ETAE,HMAXE}*. The compiler is not clever enough yet
  and needs to be told explicitly that this is a list of *double*s,
  hence the *(double[])*
  [type casting](http://en.wikipedia.org/wiki/Type_conversion). 
  
  The function then returns the number of cells refined. */

  astats s = adapt_wavelet ({eta, hmax}, (double[]){ETAE,HMAXE},
			    maxlevel, MINLEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  return s.nf;
#else // Cartesian
  return 0;
#endif
}

/**
## Initial conditions

We first specify the terrain database to use to reconstruct the
topography $z_b$. This KDT database needs to be built beforehand. See the
[*xyz2kdt* manual](http://gerris.dalembert.upmc.fr/xyz2kdt.html)
for explanations on how to do this.

We then consider two cases, either we restart from an existing
snapshot or we start from scratch.

The next line tells the Saint-Venant solver to conserve water surface
elevation rather than volume when adapting the mesh. This is important
for tsunamis since most of the domain will be close to "lake-at-rest"
balance. */

event init (i = 0)
{
  terrain (zb, "~/terrain/etopo2", NULL);

#if MULTIGRID  
  if (restore (file = "restart"))
    conserve_elevation();
  else
#endif // MULTIGRID
  {
    conserve_elevation();
    
    /**
    The initial still water surface is at $z=0$ so that the water depth
    $h$ is... */
    
    foreach()
      h[] = max(0., - zb[]);
    
    /**
    The initial deformation is given by an Okada fault model with the
    following parameters. The *iterate = adapt* option will iterate this
    initialisation until our *adapt()* function above returns zero
    i.e. until the deformations are resolved properly. */
  
    fault (x = 94.57, y = 3.83,
	   depth = 11.4857e3, // fixme: this should be a length, but okada.h does not seem to be dimensionally consistent!!
	   strike = 323, dip = 12, rake = 90,
	   length = 220e3, width = 130e3,
	   U = 18,
	   iterate = adapt);
  }
}

/**
The 4 other fault segments are triggered at the appropriate times
(seconds converted to minutes). */

event fault2 (t = 272./60.) {
  fault (x = 93.90, y = 5.22,
	 depth = 11.4857e3,
	 strike = 348, dip = 12, rake = 90,
	 length = 150e3, width = 130e3,
	 U = 23,
	 iterate = adapt);
}

event fault3 (t = 588./60.)
{
  fault (x = 93.21, y = 7.41,
	 depth = 12.525e3,
	 strike = 338, dip = 12, rake = 90,
	 length = 390e3, width = 120e3,
	 U = 12,
	 iterate = adapt);
}

event fault4 (t = 913./60.)
{
  fault (x = 92.60, y = 9.70,
	 depth = 15.12419e3,
	 strike = 356, dip = 12, rake = 90,
	 length = 150e3, width = 95e3,
	 U = 12,
	 iterate = adapt);
}

event fault5 (t = 1273./60.)
{
  fault (x = 92.87, y = 11.70,
	 depth = 15.12419e3,
	 strike = 10, dip = 12, rake = 90,
	 length = 350e3, width = 95e3,
	 U = 12,
	 iterate = adapt);
}

/**
## Outputs

### At each timestep

We output simple summary statistics for *h* and *u.x* on standard
error. */

event logfile (i++) {
  stats s = statsf (h);
  norm n = normf (u.x);
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt speed tn\n");
  fprintf (stderr, "%g %d %g %g %g %g %g %g %g %ld\n",
	   t, i, s.min, s.max, s.sum, n.rms, n.max, dt, perf.speed, grid->tn);

  /**
  We also use a simple implicit scheme to implement quadratic bottom
  friction i.e.
  $$
  \frac{d\mathbf{u}}{dt} = - C_f|\mathbf{u}|\frac{\mathbf{u}}{h}
  $$
  with $C_f=10^{-4}$. */

  double C_f = 1e-4;
  foreach() {
    double a = h[] < dry ? HUGE : 1. + C_f*dt*norm(u)/h[];
    foreach_dimension()
      u.x[] /= a;

    /**
    That is also where we update *hmax*. */

    if (h[] > dry && h[] + zb[] > hmax[])
      hmax[] = h[] + zb[];
  }
}

/**
### Snapshots

Every 60 minutes, the $h$, $z_b$ and *hmax* fields are interpolated
bilinearly onto a *n x n* regular grid and written on standard
output. */

event snapshots (t += 60; t <= 600) {
#if !BENCHMARK
#if !_MPI
  printf ("file: t-%g\n", t);
  output_field ({h, zb, hmax}, stdout, n = 1 << maxlevel, linear = true);
#endif
  
  /**
  We also save a snapshot file we can restart from. */

#if TREE
  dump (file = "dump");
#endif
#endif // !BENCHMARK
}

/**
After completion of the simulation, doing

~~~bash
make tsunami/plots
~~~

will generate the [inline plots](/src/test/README#inline-plots) such
as this one:

~~~gnuplot Maximum wave elevation (metres) reached over 10 hours.
# We first split the large standard output file into its subfiles (as
# defined by the "file:" keyword, see tsunami.c

! awk '{ if ($1 == "file:") file = $2; else print $0 > file; }' < out

set term pngcairo enhanced size 1024,1024 font ",8"
set output 'maximum.png'

# this sets the color palette to "jet"
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1, \
    	              0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, \
		      0.625 1 0.9333 0, 0.75 1 0.4392 0, \
		      0.875 0.9333 0 0, 1 0.498 0 0 )

# here we plot the value of hmax ($5) but only for wet cells 
# (i.e. h = $3 > 1e-3), dry cells take a 'no data' value i.e. 1e1000
# We use only the last output file (i.e. 't-600' minutes = 10 hours)

unset key
set size ratio -1
set xrange [67:105]
set yrange [-10:25]
set pm3d map
set logscale cb
set cbrange [0.1:10]
splot 't-600' u 1:2:($3 > 1e-3 ? $5 : 1e1000)

# we trim the image
! mogrify -trim maximum.png

# we remove the subfiles of tsunami/out
! rm -f t-*
~~~

### Movies

This is done every minute (*t++*).

We use the *mask* option of *output_ppm()* to mask out the dry
topography. Any part of the image for which *m[]* is negative
(i.e. for which *etam[] < zb[]*) will be masked out. */

#if !BENCHMARK
event movies (t++) {
  scalar m[], etam[];
  foreach() {
    etam[] = eta[]*(h[] > dry ? 1. : 0.);
    m[] = etam[] - zb[];
  }
  output_ppm (etam, mask = m, min = -2, max = 2, n = 512, linear = true,
	      fps = 30, file = "eta.mp4");

  /**
  After completion this will give the following animation
  
  ![Animation of the wave elevation. Dark blue is -2 metres and
  less. Dark red is +2 metres and more.](tsunami/eta.mp4)
  
  We also use the *box* option to only output a subset of the domain
  (defined by the lower-left, upper-right coordinates). */
  
  output_ppm (etam, mask = m, min = -2, max = 2, n = 512, linear = true,
	      fps = 30, box = {{91,5},{100,14}}, file = "eta-zoom.mp4");

  /**
  ![Animation of the wave elevation. Dark blue is -2 metres and
  less. Dark red is +2 metres and more.](tsunami/eta-zoom.mp4)
  
  And repeat the operation for the level of refinement...*/

#if TREE  
  scalar l = etam;
  foreach()
    l[] = level;
  output_ppm (l, min = MINLEVEL, max = maxlevel, n = 512, file = "level.mp4");
#endif
  
  /**
  ![Animation of the level of refinement. Dark blue is 5 and dark red
  is 10.](tsunami/level.mp4)
  
  ...and for the process id for parallel runs. */
  
#if (_OPENMP || _MPI) && !_GPU
  foreach()
    etam[] = tid();
  double tmax = npe() - 1;
  output_ppm (etam, max = tmax, n = 512, file = "pid.mp4");
#endif // _OPENMP || _MPI
}

/**
![Animation of the OpenMP process id.](tsunami/pid.mp4)

### Tide gauges

We define a list of file names, locations and descriptions and use the
*output_gauges()* function to output timeseries (for each timestep) of
$\eta$ for each location. */

Gauge gauges[] = {
  // file   lon      lat         description
  {"coco", 96.88,  -12.13, "Cocos Islands, Australia"},
  {"colo", 79.83,    6.93, "Colombo, Sri Lanka"},
  {"hani", 73.18,    6.77, "Hanimaadhoo, Maldives"},
  {"male", 73.52,    4.18, "Male, Maldives"},
  {"gana", 73.17,   -0.68, "Gan, Maldives"},
  {"dieg", 72.38,    -7.3, "Diego Garcia, UK"},
  {"rodr", 63.42,  -19.67, "Rodriguez I., Mauritius"},
  {"loui", 57.5,   -20.15, "Port Louis, Mauritius"},
  {"lare", 55.3,   -20.92, "La Reunion, France"},
  {"hill", 115.73, -31.82, "Hillarys, Australia"},
  {"sala", 54,         17, "Salalah, Oman"},
  {"laru", 55.53,   -4.68, "Pointe La Rue, Seychelles"},
  {"lamu", 40.9,    -2.27, "Lamu, Kenya"},
  {"zanz", 39.18,   -6.15, "Zanzibar, Tanzania"},
  {"chen", 80.3,     13.1, "Chennai, India"},
  {"para", 86.7,    20.26, "Paradip, India"},
  {"visa", 83.28,   17.68, "Visakhapatnam, India"},
  {"koch", 76.26,    9.96, "Kochi, India"},
  {"morm", 73.8,    15.42, "Mormugao, India"},
  {"okha", 69.08,   22.47, "Okha, India"},
  {"tuti", 78.15,     8.8, "Tuticorin, India"},
  {"taru", 99.65,   6.702, "Tarutao, Thailand"},
  {"tapa", 98.425,  7.765, "Tapaonoi, Thailand"},
  {NULL}
};

event gauges1 (i++) output_gauges (gauges, {eta});

/**
As before gnuplot processes these files to produce this image:

~~~gnuplot Comparison between observed and simulated timeseries (hours) of wave elevations (metres) for a selection of tide gauges.
reset
set term svg enhanced size 625,800 font ",10"
set multiplot layout 5,1 scale 1,1.1
set xrange [3:8]
set key bottom right
set title 'Hanimaadhoo, Maldives'
plot 'hani' u ($1/60.):2 w l t 'modelled', \
     '../hanires.txt' u 1:($2/100.) w lp t 'observed'
unset key
set title 'Male, Maldives'
plot 'male' u ($1/60.):2 w l t 'modelled', \
     '../maleres.txt' u 1:($2/100.) w lp t 'observed'
set title 'Gan, Maldives'
plot 'gana' u ($1/60.):2 w l t 'modelled', \
     '../ganares.txt' u 1:($2/100.) w lp t 'observed'
set title 'Diego Garcia'
plot 'dieg' u ($1/60.):2 w l t 'modelled', \
     '../diegres.txt' u 1:($2/100.) w lp t 'observed'
set title 'Columbo, Sri Lanka'
set xrange [2.5:8]
plot 'colo' u ($1/60.):2 w l t 'modelled', \
     '../colores.txt' u 1:($2/100.) w lp t 'observed'
unset multiplot
~~~

### Google Earth KML file

We also generate images and a [Keyhole Markup Language]() file which
can be imported into [Google Earth](https://www.google.com/earth/) to
superpose the evolving wave height field on top of Google Earth
data. */

event kml (t += 15)
{
  static FILE * fp = fopen ("eta.kml", "w");
  if (t == 0)
    fprintf (fp,
	     "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
	     "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n"
	     "  <Folder>\n");
  fprintf (fp,
	   "    <GroundOverlay>\n"
	   "      <TimeSpan>\n"
	   "        <begin>2004-12-26T%02d:%02d:00</begin>\n"
	   "        <end>2004-12-26T%02d:%02d:00</end>\n"
	   "      </TimeSpan>\n"
	   "      <Icon>\n"
	   "	    <href>eta-%g.png</href>\n"
	   "      </Icon>\n"
	   "      <LatLonBox>\n"
	   "	    <north>35</north>\n"
	   "	    <south>-19</south>\n"
	   "	    <east>121</east>\n"
	   "	    <west>67</west>\n"
	   "      </LatLonBox>\n"
	   "    </GroundOverlay>\n",
	   8 + ((int)t)/60, ((int)t)%60,
	   8 + ((int)t + 15)/60, ((int)t + 15)%60,
	   t);
  fflush (fp);
  scalar m[], etam[];
  foreach() {
    etam[] = eta[]*(h[] > dry ? 1. : 0.);
    m[] = etam[] - zb[];
  }
  char name[80];
  sprintf (name, "eta-%g.png", t);
  output_ppm (etam, file = name, mask = m,
	      min = -2, max = 2, n = 1 << maxlevel, linear = true);
  if (t == 600)
    fprintf (fp, "</Folder></kml>\n");
}
#endif // !BENCHMARK

/**
## Adaptivity

And finally we apply our *adapt()* function at every timestep. */

event do_adapt (i++) adapt();
