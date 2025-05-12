/**
# Estuary of the river Guadalquivir

The river Guadalquivir is an old river with a lot of history. Despite
its low depth the river has been always navigable. In fact, it was
used by the vikings to atack Sevilla in the IX century.  Today, the
river is still navigable up to Sevilla where a port is in use. Ships
with large drafts must ride the tide to sail up to the port of
Sevilla. 

![The Estuary of the Guadalquivir](estuary.jpg)

The river discharges in the Atlantic sea. The Doñana
National Park is located on the right riverbank. The village of
Sanlucar de Barrameda is on the left bank and is famous for its
wines.

# Simulation 

The simulation of the estuary is performed using an explicit scheme for
the Saint-Venant shallow water equations with adaptation.  That means to
use the most general model [hydro.c](/src/layered/hydro.h) assuming a
single layer and that the vertical momentum equation reduces to the
hydrostatic relationship. These assumptions are the default one. So
nothing has to be set in addition. Details of the layered model can be
found in [Popinet, 2020](/Bibliography#popinet2020). */

#include "grid/multigrid.h"
#if !SV
#include "layered/hydro.h"
#else
#include "saint-venant.h"
#endif
#include "terrain.h"

/**
*riemann.h* contains kurganov() which is required in *discharge.h*. */

#if !SV
#include "riemann.h"
#endif
#include "discharge.h"

#define MINLEVEL 4
int maxlevel = 8;
scalar river[];

/**
The sea goes up and down sinusoidally because of a semidiurnal moon
tide M2 of amplitude equal to 0.7 m. The period of
this component is of 12.42 h (or 745.2 min). We assume that the origin
of times is in phase with the tide. We have decided to measure the
lengths in meters and the times in minutes. There is no problem with
it meanwhile we were consistent with the system of units elected. */

#define OM  (2*pi/745.2) // frequency (The M2 period is 12.42h)
#define PHASE 0. // Time delay (in case the time origin is shifted)
#define AMP 0.7 // Amplitude

u.n[bottom] = - radiation(AMP*sin(OM*t + PHASE));
u.n[top]    = + radiation(AMP*sin(OM*t + PHASE));
u.n[left]   = - radiation(AMP*sin(OM*t + PHASE));

/**
In the main() we define the left bottom corner of the domain (with
origin()) and the size of the computational domain of 21000 x 21000
$m^2$. Note we are using a cartesian UTM proyection (ED50 H 30N).*/

int main()
{
  L0 = 21000 [0];
  DT = HUGE [0];
  origin (189000, 4071000);
  G = 9.81*sq(60); //The time is measured in minutes 
  init_grid (1 << maxlevel);
  run();
}

/** 
In the init event we initialize the bathymetry *zb* using the database
*guadalquivir*. The bathymetry must go beyond the simulation domain.
Un particular must be availabe even for the coarsest ghost cells. 
If not a fake ground is set in the boundary (run the present test with 
minlevel equal to 2 and see). The folder where the database locates 
must be well specified. Also, conserve_elevation() serves to keep 
the depth of the fluid instead of the fluid volume when 
coarsening/refining the grid. If not, gravity waves are induced in every change of grid. */

event init (i = 0)
{
  terrain (zb, "/home/basilisk/terrain/guadalquivir", NULL);

  if (restore (file = "dump"))
    conserve_elevation();
  else {
    conserve_elevation();
    
    foreach() {
      h[] = max(0., - zb[]);
      river[] = (y > 4087935  && y < 4088357);
    }
    boundary ({h});
  }
  output_ppm (zb, file = "zb.png");
}

/** 
The bathymetry is initially described with a grid size equal to
$\Delta = L0/2^N \simeq 82 m$ where $N=8$ is the level of refinement
used.

![Bathymetry](guadalquivir/zb.png)

The course of the river is very well described because of the high
number of points in the database at that location. Not the same
happens with some dry parts. The bathymetry also shows that the low
curse of the Guadalquivir goes through a marshland.

An inflow enters through the right. The inflow depends on time. Over a
constant flow rate $Qo$ = 150 $m^3/s$ a flood having the form of a
gaussian bell of height Qp = 2000 $m^3/s$ is set. The peak of flood
enters the domain with a delay of 500 minutes after the first peak of
the tide. If the riverbed is not well distinguished from the
riverbank, the inflow could spread along the right boundary as it trys
to compute a water elevation $\eta$ compatible with the inflow. To
avoid this unwanted mechanism we force the inflow to enter through the
region for which *river* is equal to 1. That is, through the narrow
gap going from $y_{low}$ = 4087935 up to $y_{up}$ = 4088357. */

#define gaussian(t,to,b) (exp(-sq((t-to)/(2.*b))))
#define Qo (150*60)
#define Qp (2000.*60)

double etar;
event inflow (i++) {
  double qlaw = Qo + Qp*gaussian (t, 500, 100);
  etar = eta_b (qlaw,right, river, 1);

  h[right] = (river[] ? max (etar - zb[], 0.) : 0.);
  eta[right] = (river[] ? max (etar - zb[], 0.) : 0.) + zb[];
}

/**
Friction slow down the fluid according to,

$$ \frac{d \mathbf{u}}{d t} = - \frac{\lambda}{4 r_h} |\mathbf{u}| \mathbf{u} 
= - \frac{c_f}{r_h} |\mathbf{u}| \mathbf{u}
$$

where $\lambda$ is the Darcy friction factor. $c_f = \lambda/4$ is the
Fanning friction factor.  The hydraulic radius is in this case $r_h =
h$. We set $c_f = 10^{-3}$. The above equation is linearized and
solved implicitly,

$$ 
\frac{\mathbf{u}^{n+1} - \mathbf{u}^n}{\Delta t} 
= - \frac{c_f}{h} |\mathbf{u}|^n \mathbf{u}^{n+1}
$$
 */

#define fa 1e-3 // The Fanning friction factor

event friction (i++) {
  foreach() {
    double a = h[] < dry ? HUGE : 1. + fa*dt*norm(u)/(h[]);
    foreach_dimension()
      u.x[] /= a;
  }
  boundary ((scalar *){u});
}

/**
## Adaption

We adapt the grid after the water elevation. */

#if TREE
event adapt (i++) {
  scalar etap[];
  foreach() 
    etap[] = h[] > dry ? eta[] : 0.;
  boundary ({etap});
  adapt_wavelet ({etap}, (double[]){3e-4}, maxlevel, MINLEVEL);
}
#endif

/**
## Outputs

We save each 100 time step a dump file in case we need to restart the
simulation. */

event save_dump( i += 100) {
  dump("dumpfile");
}

/**
Also, we save in ASCII ESRI format the inundation area at instant t = 600.*/

event inundation(t = 600) {
  scalar hp[];
  foreach()
    hp[] = h[] > dry ? h[] : nodata;
  
  FILE * fp = fopen ("inundation.asc", "w");
  output_grd(hp, fp);
  fclose(fp);
}

/**
Also, we keep a time record (with a time step of 1 min) of the
amount of water in the domain as well as we monitorize the relevant
variables near the "Salina of Monte Algaida". */

event logfile(t++) {
  double vol = 0;
  foreach(reduction(+:vol))
    vol += h[]*sq(Delta);
  foreach_point (203000, 4089170) //Near "salina del monte algaida"
    fprintf(stderr, "%g %g %g %g %g\n",
	  t, vol, h[], eta[], sign(u.x[])*sqrt(sq(u.x[])+ sq(u.y[])));
  fflush(stderr);
}

/**
~~~gnuplot  top: velocity at Salina del monte Algaida, middle: water elevation low: Time dependece of the tides and the flood
om=2.*pi/745.2
gaussian(x)=exp(-((x-500)/(2.*100))**2)
set multiplot layout 3, 1
set arrow from 0.,0. to 1200., 0. nohead front lc rgb "black" lw 2  dashtype "-"
plot 'log' u 1:($5/60) w l lw 2 t 'velocity (m/s)', '../guadalquivir-sv/log' u 1:($5/60) w l lw 2 t 'velocity (SV)'
plot 'log' u 1:4 w l lw 2 t 'eta', '../guadalquivir-sv/log' u 1:4 w l lw 2 t 'eta (SV)'
set xrange [0 : 1200]
plot sin(om*x) lw 2 t 'tide', gaussian(x) lw 2 t 'flood'
set xlabel "t (min)"
unset multiplot
~~~

Finally,  we record a couple of movies illustrating the time evolution 
of the elevation $\eta$ and the norm of the velocity $|u|$. 
*/

#include "view.h"
event snapshot(t++; t < 1200) {
  char str[40];
  sprintf (str, "t = %g ", t);
  
  scalar etap[], vmax[];
  foreach() {
    etap[] = h[] > dry ? eta[] : nodata;
    vmax[] = h[] > dry ? norm(u) : nodata;
  }
  boundary({etap, vmax});

  /** 
  To center the graph in the canvas we set the translation tx = -X0/L0
  - 0.5 and ty = -Y0/L0 - 0.5 where X0 and Y0 son the coordinates of
  the left bottom corner. */

  view (fov = 22, tx = -X0/L0 - 0.5, ty = -Y0/L0 - 0.5);
  draw_string (str, 2, size = 25, lw = 1); 
  squares("etap", linear = false);
  save ("eta.mp4");
  clear();
  squares("vmax", linear = false);
  vectors("u", scale = 2);
  draw_string (str, 2, size = 25, lw = 1);
  save ("vel.mp4");
}

/**
![Velocity at the Salina](guadalquivir/vel.mp4)

![Elevation at the Salina](guadalquivir/eta.mp4)
*/
