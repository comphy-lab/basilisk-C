/**
#The Toce urban case
 
## Compilation
To compile, type in the terminal: 
qcc ToceUrban.c -o ToceUrban.x -lm ~/basilisk/src/kdt/kdt.o

## Introduction 
The second case of validation is reproducing "The Model city flooding experiment benchmark". This model is build on the first $15 \ m$ of the fluvial Toce model. The authors added $20$ buildings distributed in 4 aligned rows, 9 gauge stations are distributed around the buildings and at the entrance of the flow. The imposed entry condition is again a hydrograph. The flow rate goes from 0 to a maximum of $130\  l.s^{-1}$ in 4 seconds and then progressively decreases to $30\ l.s^{-1}$ in 50 seconds. This condition reproduces the typical water flow of a flash flood. The experiment lasts 60 seconds.
*/



/**
## Include and def
*/
#include "b-flood/saint-venant-topo.h"
#include "b-flood/manning.h"
#include "b-flood/inflow.h"
#include "b-flood/hydro.h"
#include "terrain.h"
#include "input.h"

#define MAXLEVEL 10 // 1024 Volume
#define MINLEVEL 6
#define EMAXE 1e-4

int nf,nc,ncell;
timer t0;


/**
## Parameters  
*/
int main() {
  n = 0.0162;
  L0 = 14.9;
  X0 = 0.1;
  Y0 = 0;
  G = 9.81;
  dry = 1e-6;
  N = 1 << MAXLEVEL; // Which is equivalent to N = 2^(maxlevel)
  run();
}

event stop(t = 60)
  fprintf(stderr,"# The end\n");


/**
## Initial condition
   
*/
event init (i = 0){
  terrain (zb, "./Topo/topo-a", NULL);
  DT = 0.05;
  t0 = timer_start();
  double dx = L0/N;
 foreach(){
    if( y <= 8.4 && y >= 5 && x <= 32*dx)
      river[] = 1;
  }
  boundary({zb,river});

  static FILE * rpz = fopen ("raster_zb-stag.asc", "w");
  output_grd (zb, rpz);
}

/**
## Inflow boundary condition

On the left, we impose the following flow rate using the "discharge" function, 

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
event bound_left( i++ ){
  l1 =  hydrograph ( "hydro.asc", 1);
  if( l1 < 0 )
    l1 = 0;
  l1 /= 1000.; // convertir les litres en metres cube
  eta1 = eta_b (l1,left,1);
  
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
 
  static FILE * fpp = fopen("eta.dat","w");
  if( i == 0)
     fprintf(fpp,"#t Q Qread EtaImp\n");
  fprintf(fpp,"%lf %lf %lf %lf \n",t,Q,l1,eta1);
  fflush(fpp);
  
}

/**
## Outflow boundary condition
We set a free exit condition on the right edge of the domain.
*/

u.t[right] = neumann(0);
u.n[right] = neumann(0);
h[right] = neumann(0);
eta[right] = neumann(0);


/**
## Adaptative refinement
   
*/
int adapt_H() {
  scalar nh[];
  foreach()
    nh[] = h[] > dry ? h[] : 0;
  boundary ({nh});
  astats s = adapt_wavelet ({nh}, (double[]){EMAXE} ,MAXLEVEL,MINLEVEL);
  nc += s.nc;
  nf += s.nf;
  return s.nf;
}
event do_adapt_H (i++){
  if( t > 1 )
    adapt_H();
}
//////// OUTPUT
/**
## Logfile
*/

event logfile (i += 10) {
  timing tr = timer_timing (t0, i, 0, NULL);
  stats s = statsf (h);
  scalar v[];
  foreach()
    v[] = norm(u);
  stats no = statsf (v);
  if( i == 0 )
    fprintf (stderr, "#1t 2treal 3i 4h.min 5h.max 6h.sum 7u.min 8u.max 9dt\n");  
  fprintf (stderr, "%g %g %d %g %g %g %g %g %g\n", t,tr.real, i, s.min, s.max, s.sum, 
	   no.min, no.max, dt);
}
/** Gauges 

*/
Gauge gauges[] = {  
  // gauge   x         y       zb
  //{"./Gauges/P1",-0.572923,7.618307, "P1"}, 
  {"P2",0.280560,7.585584, "P2"}, 
  {"P3",3.279013,7.324738, "P3"}, 
  {"P4",3.153249,7.834554, "P4"}, 
  {"P5",3.534738,7.748334, "P5"}, 
  {"P6",3.916227,7.662114, "P6"}, 
  {"P7",4.128212,7.533788, "P7"}, 
  {"P8",4.044554,7.874098, "P8"}, 
  {"P9",4.256538,7.745772, "P9"}, 
  {"P10",4.638028,7.659552, "P10"}, 
  {NULL}
};
event gauges1 (t += 0.2) output_gauges (gauges, {h});
/**
Gauges Results 

~~~gnuplot P1
reset
set xlabel 'T(s)' 
set ylabel 'Level (m)' 
set xtics; set ytics
set key bot
plot 'P1' w lp t'P1- numerical', 'P1.ref' w lp t'P1 - measured'
~~~

~~~gnuplot P2
reset
set xlabel 'T(s)' 
set ylabel 'Level (m)' 
set xtics; set ytics
set key bot
plot 'P2' w lp t'P2- numerical', 'P2.ref' w lp t'P2 - measured'
~~~

~~~gnuplot P3
reset
set xlabel 'T(s)' 
set ylabel 'Level (m)' 
set xtics; set ytics
set key bot
plot 'P1' w lp t'P3- numerical', 'P3.ref' w lp t'P3 - measured'
~~~

~~~gnuplot P4
reset
set xlabel 'T(s)' 
set ylabel 'Level (m)' 
set xtics; set ytics
set key bot
plot 'P4' w lp t'P4- numerical', 'P4.ref' w lp t'P4 - measured'
~~~

~~~gnuplot P5
reset
set xlabel 'T(s)' 
set ylabel 'Level (m)' 
set xtics; set ytics
set key bot
plot 'P5' w lp t'P5- numerical', 'P5.ref' w lp t'P5 - measured', 'P5Byun.ref' w lp t'P5 - Kim'
~~~

~~~gnuplot P6
reset
set xlabel 'T(s)' 
set ylabel 'Level (m)' 
set xtics; set ytics
set key bot
plot 'P6' w lp t'P6- numerical', 'P6.ref' w lp t'P6 - measured'
~~~

~~~gnuplot P7
reset
set xlabel 'T(s)' 
set ylabel 'Level (m)' 
set xtics; set ytics
set key bot
plot 'P7' w lp t'P7- numerical', 'P7.ref' w lp t'P7 - measured'
~~~

~~~gnuplot P8
reset
set xlabel 'T(s)' 
set ylabel 'Level (m)' 
set xtics; set ytics
set key bot
plot 'P8' w lp t'P8- numerical', 'P8.ref' w lp t'P8 - measured'
~~~

~~~gnuplot P9
reset
set xlabel 'T(s)' 
set ylabel 'Level (m)' 
set xtics; set ytics
set key bot
plot 'P9' w lp t'P9- numerical', 'P9.ref' w lp t'P9 - measured'
~~~

~~~gnuplot P10
reset
set xlabel 'T(s)' 
set ylabel 'Level (m)' 
set xtics; set ytics
set key bot
plot 'P10' w lp t'P10- numerical', 'P10.ref' w lp t'P10 - measured'
~~~
*/

/**
## Movies

We record movies of the water depth, the velocity and the level of
refinement.

![Evolution of the water level](ToceUrban/height.mp4)( width="600px" )*/
event movies ( t += 0.1 ) {
  scalar m[], v[], fr[];
  
  foreach(){
    v[] = norm(u);
    m[] = (h[] > dry ) - 0.5;
    fr[] = h[] > dry ? v[]/sqrt(G*h[]) : 0;
  }    
  boundary ({m,v,fr});
  output_ppm (h, mask = m, min = 0, max = 0.2, n = N, linear = true, file = "height.mp4");
  output_ppm (v, mask = m, min = 0, max = 4, n = N, linear = true, file = "vel.mp4");
  output_ppm (fr, mask = m, min = 0, max = 2, n = N, linear = true, file = "froude.mp4");

  scalar l = m;
  foreach()
    l[] = level;
  output_ppm (l, min = MINLEVEL, max = MAXLEVEL, n = N, file = "level.mp4");
}



event ending(t = end){
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

 
event figures (t += 1)
{
  scalar m[], etam[];
  foreach() {
    etam[] = eta[]*(h[] > dry);
    m[] = etam[] - zb[];
  }
  boundary ({m, etam});
  char name[80];
  sprintf (name, "h-%g-1.png", t);
  output_ppm (h, mask = m, min = 0, max = 0.2, file = name, n = 1024,
	      linear = true,
	      opt = "-fill white -opaque black");
  sprintf (name, "h-%g-2.png", t);
  output_ppm (h, mask = m, min = 0, max = 0.2, file = name, n = 1024,
	      linear = true);
  
  
  sprintf (name, "level-%g.png", t);
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, min = MINLEVEL, max = MAXLEVEL, file = name, n = 1024,
	      linear = false,);
}
/**
## Link to the homepage
* [Homepage](Readme)	    		  
   		  
*/