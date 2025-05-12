/**
# The Toce fluvial case


## Compilation
To compile, type in the terminal: 
qcc ToceFluvial.c -o ToceFluvial.x -lm ~/basilisk/src/kdt/kdt.o

## Introduction
The entire river is 50 meters long and 11 meters large. The DTM is at a resolution of 5 cm, note that this includes the reproduction of houses. The imposed inlet condition is an hydrograph. The hydrograph is composed of a brutal rising stage from 0 to 210 $l.s^{-1}$ then by a slower and continuous descent phase which reaches up to $60\ l.s^{-1}$ at the end of the experiment.


## Include and def */

#include "b-flood/saint-venant-topo.h"
#include "b-flood/manning.h"
#include "b-flood/inflow.h"
#include "b-flood/hydro.h"
#include "terrain.h"
#include "input.h"

#define MAXLEVEL 10 // 1024 Volume
#define MINLEVEL 6
#define EMAXE 5e-4

int nf, nc, ncell;
timer t0;

/**
## Parameters */

int main()
{
  n = 0.0162;
  L0 = 43.4;
  X0 = 0.1;
  Y0 = 0;
  G = 9.81;
  dry = 1e-6;
  N = 1 << MAXLEVEL; // Which is equivalent to N = 2^(maxlevel)
  run();
}

// The case lasts 180 seconds
event stop (t = 180)
  fprintf (stderr, "# The end\n");
  
/**
## Initial condition */

event init (i = 0)
{
  // We record the topography in a special way (k-tree) so that we can
  // use adaptive refinement on it.  
  terrain (zb, "./Topo/topo-house", NULL);
  
  // We have to set a time step because there is no water at the beginning.
  DT = 1;
  
  t0 = timer_start();
  double dx = L0/N;

  // We set the interval through which the water flow will come in.
  foreach()
    if (y <= 8.4 && y >= 5 && x <= 8*dx)
      river[] = 1;
  boundary ({river});
}

/**
## Inflow boundary condition

On the left, we impose the flow rate using the "discharge" function. 

~~~gnuplot Hydrograph
reset
set xlabel 'T(s)' 
set ylabel 'Inflow (cumsec)' 
set xtics; set ytics
set key bot
plot './hydro.asc' w l t'hydrograph'
~~~
*/

u.t[left] = dirichlet(0);
u.n[left] = neumann(0);

double l1, eta1, vola = 0;
event bound_left( i++ )
{
  l1 =  hydrograph ("hydro.asc", 1);
  if (l1 < 0)
    l1 = 0;
  eta1 = eta_b(l1,left,1);
  
  h[left] = dirichlet( max( river[] == 1 ? eta1 - zb[] : 0 , 0) );
  eta[left] = dirichlet( max( river[] == 1 ? eta1 : zb[]  , zb[]) );

  double vol = 0;
  foreach(reduction(+:vol)){
    if( h[] > dry )
      vol += h[] * sq(Delta);
  }
  double Q = 0;
  if( dt > 0 )
    Q =  (vol- vola)/dt;
  vola = vol;
 
  static FILE * fpp = fopen ("eta.dat", "w");
   if (i == 0)
     fprintf (fpp, "#t Q Qread EtaImp\n");
  fprintf (fpp, "%lf %lf %lf %lf \n", t, Q, l1, eta1);
  fflush (fpp);  
}

/**
## Outflow boundary conditions

We set a free exit condition on the right edge of the domain. */

u.t[right] = neumann(0);
u.n[right] = neumann(0);
h[right] = neumann(0);
eta[right] = neumann(0);

/**
We take care to reduce the maximum time step just before the
water inlet starts. */

event timestepreduce (t = 17) {
  DT = 0.01;
}

/**
## Adaptive refinement 

The adaptive refinement starts just before the water inlet starts.
*/

event adapt (t = 17; i++)
{
  scalar nh[];
  foreach()
    nh[] = h[] > dry ? h[] : 0;
  boundary ({nh});
  astats s = adapt_wavelet ({nh}, (double[]){EMAXE}, MAXLEVEL, MINLEVEL);
  nc += s.nc;
  nf += s.nf;
}

/**
## Logfile */

event logfile (i += 10)
{
  timing tr = timer_timing (t0, i, 0, NULL);
  stats s = statsf (h);
  scalar v[];
  foreach()
    v[] = norm(u);
  
  stats no = statsf (v);
  if( i == 0 )
    fprintf (stderr, "#1t 2treal 3i 4h.min 5h.max 6h.sum 7u.min 8u.max 9dt\n");  
  fprintf (stderr, "%g %g %d %g %g %g %g %g %g\n",
	   t,tr.real, i, s.min, s.max, s.sum,  no.min, no.max, dt);
}


/**
## Gauges 

We record the water level at the same location than the experiment*/
Gauge gauges[] = {  
  // gauge   x         y       zb
  {"S2",0.280560,7.585584, "#S2"},
  {"S4",4.255610,7.745400, "#S4"},
  {"S6S",6.753537,9.831888, "#S6S"},
  {"S6D",7.552486,7.696863, "#S6D"},
  {"S8D",10.743772,7.788068, "#S8D"},
  {"P1",2.113886,7.866009, "#P1"},
  {"P2",3.931183,8.893101, "#P2"},
  {"P3",4.527122,6.587255, "#P3"},
  {"P4",7.086251,8.768313, "#P4"},
  {"P5",10.163515,10.218714, "#P5"},
  {"P8",14.950365,11.804039, "#P8"},
  {"P9",15.684258,13.802726, "#P9"},
  {"P10",18.190768,10.429219, "#P10"},
  {"P13",19.814693,11.984270, "#P13"},
  {"P18",23.848138,14.744612, "#P18"},
  {"P19",25.469928,16.186758, "#P19"},
  {"P21",30.442363,18.357176, "#P21"},
  {"P23",36.466741,18.750427, "#P23"},
  {"P24",39.024308,22.827062, "#P24"},
  {"P25",40.443587,21.238121, "#P25"},
  {"P26",40.965033,26.182320, "#P26"},
  {NULL}
};
event gauges1 (t += 1) output_gauges (gauges, {h});
/**

## Gauges results

~~~gnuplot Gauges P1, P2 and P3
set term svg enhanced size 800,400 font ",10"
set ylabel 'Level (m)' 
set xtics; set ytics
set multiplot
set size .95,.35
set key bot
set origin 0.02, 0.6
plot 'P1' w lp t'P1 - numerical', 'P1.ref' w lp t'P1 - measured'
set origin 0.02, 0.3
plot 'P2' w lp t'P2 - numerical', 'P2.ref' w lp t'P2 - measured'
set origin 0.02, .01
set xlabel 'T(s)' 
plot 'P3' w lp t'P3 - numerical', 'P3.ref' w lp t'P3 - measured'
unset multiplot
~~~ 
 
~~~gnuplot Gauges P4, P5 and P8
set term svg enhanced size 800,400 font ",10"
set ylabel 'Level (m)' 
set xtics; set ytics
set multiplot
set size .95,.35
set key bot
set origin 0.02, 0.6
plot 'P4' w lp t'P4 - numerical', 'P4.ref' w lp t'P4 - measured'
set origin 0.02, 0.3
plot 'P5' w lp t'P5 - numerical', 'P5.ref' w lp t'P5 - measured'
set origin 0.02, .01
set xlabel 'T(s)' 
plot 'P8' w lp t'P8 - numerical', 'P8.ref' w lp t'P8 - measured'
unset multiplot
~~~ 

~~~gnuplot Gauges P9, P10 and P13
set term svg enhanced size 800,400 font ",10"
set ylabel 'Level (m)' 
set xtics; set ytics
set multiplot
set size .95,.35
set key bot
set origin 0.02, 0.6
plot 'P9' w lp t'P9 - numerical', 'P9.ref' w lp t'P9 - measured'
set origin 0.02, 0.3
plot 'P10' w lp t'P10 - numerical', 'P10.ref' w lp t'P10 - measured'
set origin 0.02, .01
set xlabel 'T(s)' 
plot 'P13' w lp t'P13 - numerical', 'P13.ref' w lp t'P13 - measured'
unset multiplot
~~~ 
 
~~~gnuplot Gauges P18, P19 and P21
set term svg enhanced size 800,400 font ",10"
set ylabel 'Level (m)' 
set xtics; set ytics
set multiplot
set size .95,.35
set key bot
set origin 0.02, 0.6
plot 'P18' w lp t'P18 - numerical', 'P18.ref' w lp t'P18 - measured'
set origin 0.02, 0.3
plot 'P19' w lp t'P19 - numerical', 'P19.ref' w lp t'P19 - measured'
set origin 0.02, .01
set xlabel 'T(s)' 
plot 'P21' w lp t'P21 - numerical', 'P21.ref' w lp t'P21 - measured'
unset multiplot
~~~ 

~~~gnuplot Gauges P23, P24 and P25
set term svg enhanced size 800,400 font ",10"
set ylabel 'Level (m)' 
set xtics; set ytics
set multiplot
set size .95,.35
set key bot
set origin 0.02, 0.6
plot 'P23' w lp t'P23 - numerical', 'P23.ref' w lp t'P23 - measured'
set origin 0.02, 0.3
plot 'P24' w lp t'P24 - numerical', 'P24.ref' w lp t'P24 - measured'
set origin 0.02, .01
set xlabel 'T(s)' 
plot 'P25' w lp t'P25 - numerical', 'P25.ref' w lp t'P25 - measured'
unset multiplot
~~~ 


~~~gnuplot Gauges P26, S2 and S4
reset
set term svg enhanced size 800,400 font ",10"
set ylabel 'Level (m)' 
set xtics; set ytics
set multiplot
set size .95,.35
set key bot
set origin 0.02, 0.6
plot 'P26' w lp t'P26 - numerical', 'P26.ref' w lp t'P26 - measured'
set origin 0.02, 0.3
plot 'S2' w lp t'S2 - numerical', 'S2.ref' w lp t'S2 - measured'
set origin 0.02, .01
set xlabel 'T(s)' 
plot 'S4' w lp t'S4 - numerical', 'S4.ref' w lp t'S4 - measured'
unset multiplot
~~~ 

~~~gnuplot Gauges S6D, S6S and S8D
reset
set term svg enhanced size 800,400 font ",10"
set ylabel 'Level (m)' 
set xtics; set ytics
set multiplot
set size .95,.35
set key bot
set origin 0.02, 0.6
plot 'S6D' w lp t'S6D - numerical', 'S6D.ref' w lp t'S6D - measured'
set origin 0.02, 0.3
plot 'S6S' w lp t'S6S - numerical', 'S6S.ref' w lp t'S6S - measured'
set origin 0.02, .01
set xlabel 'T(s)' 
plot 'S8D' w lp t'S8D - numerical', 'S8D.ref' w lp t'S8D - measured'
unset multiplot
~~~ 

*/

/**
## Movies

We record movies of the water depth, the velocity and the level of
refinement. 


![Evolution of the water level](ToceFluvial/height.mp4)( width="600px" )

![Evolution of the Froude number. Dark blue is zero, dark red is two.](ToceFluvial/froude.mp4)( width="600px" )

![Evolution of the refinement level](ToceFluvial/level.mp4)( width="600px" )

![Evolution of the water level close to town #1](ToceFluvial/height-town1.mp4)( width="600px" )

![Evolution of the refinement level close to town #1](ToceFluvial/level-town1.mp4)( width="600px" )

![Evolution of the water level close to town #2](ToceFluvial/height-town2.mp4)( width="600px" )

![Evolution of the refinement level close to town #2](ToceFluvial/level-town2.mp4)( width="600px" )
*/

event movies ( t= 17; t += 0.1 )
{
  scalar m[], v[], fr[];
  
  foreach(){
    v[] = norm(u);
    m[] = (h[] > dry ) - 0.5;
    fr[] = h[] > dry ? v[]/sqrt(G*h[]) : 0;
  }    
  boundary ({m,v,fr});
  
  output_ppm (h, mask = m, min = 0, max = 0.4, n = N, linear = true, file = "height.mp4");
  output_ppm (h, mask = m, min = 0, max = 0.4, n = N, linear = true,
	      box = {{0,5},{12.5,12.5}}, file = "height-town1.mp4"); 
  output_ppm (h, mask = m, min = 0, max = 0.4, n = N, linear = true,
	      box = {{15,9},{25,19}}, file = "height-town2.mp4");
  
  

  output_ppm (v,  mask = m, min = 0, max = 3, n = N, linear = true, file = "vel.mp4");

  output_ppm (fr, mask = m, min = 0, max = 2, n = N, linear = true, file = "froude.mp4");
  
  scalar l = m;
  foreach()
    l[] = level;
  
  output_ppm (l, min = MINLEVEL, max = MAXLEVEL, n = N, file = "level.mp4");
  output_ppm (l, min = MINLEVEL, max = MAXLEVEL, n = N,
	      box = {{0,5},{12.5,12.5}}, file = "level-town1.mp4");
  output_ppm (l, min = MINLEVEL, max = MAXLEVEL, n = N,
	      box = {{15,9},{25,19}}, file = "level-town2.mp4");
}


event snap (t = 65)
{
  static FILE * rph = fopen ("raster_h-t65.asc", "w");
  output_grd (h, rph);

  scalar m[], v[];
    
  foreach(){
    v[] = norm(u);
    m[] = (h[] > dry ) - 0.5;
  }    
  boundary ({m,v});
  
  static FILE * rpu = fopen ("raster_v-t65.asc", "w");
  output_grd (v, rpu);
}

event ending (t = end)
{
  static FILE * rph = fopen ("raster_h.asc", "w");
  output_grd (h, rph);

  scalar m[], v[];
    
  foreach(){
    v[] = norm(u);
    m[] = (h[] > dry ) - 0.5;
  }    
  boundary ({m,v});
  
  static FILE * rpu = fopen ("raster_v.asc", "w");
  output_grd (v, rpu);
}

/**
Figures 

Here we record the photos representing the dynamics of adaptive refinement and
water level used in the scientific article.
*/
event figures (t = 17; t <= 67; t += 2)
{
  scalar m[], etam[];
  foreach() {
    etam[] = eta[]*(h[] > dry);
    m[] = etam[] - zb[];
  }
  boundary ({m, etam});
  char name[80];
  sprintf (name, "h-%g.png", t);
  output_ppm (h, min = 0, max = 0.4, mask = m, file = name, n = 1024,
	      linear = true,
	      opt = "-fill white -opaque black");
    
  sprintf (name, "level-%g.png", t);
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, min = MINLEVEL, max = MAXLEVEL, file = name, n = 1024,
	      linear = false,);
}
/**
## References

~~~bib
Todo
~~~

## Link to the homepage

* [Homepage](Readme)
*/