/**
# Karamea flood tutorial on basilisk

This page is a step-by-step description of how to setup, run and
visualise a flood simulation with basilisk. This is the basilisk
equivalent of [this tutorial on
Gerris](http://gerris.dalembert.upmc.fr/karamea_flood_tutorial.html).
We will simulate the lower reaches of the Karamea river, including
part of the sea on the western boundary.

## Prerequisites

For this tutorial, we need three input files :

* A Digital Terrain Model of the area:
 [topo.asc](http://gerris.dalembert.upmc.fr/karamea/topo.asc)

* A flow-rate file
 defining the evolution in time of the river flow at the eastern
 boundary: [flow.asc](http://gerris.dalembert.upmc.fr/karamea/flow.asc)

* A tide file defining the sea level at the western boundary:
 [tide.asc](http://gerris.dalembert.upmc.fr/karamea/tide.asc)

Download and put them in the same repository than your C file. You
need to convert the topo.asc file in a topography file that can be
read by basilisk. [The whole process is explained
here.](http://gerris.dalembert.upmc.fr/karamea_flood_tutorial.html#Creating_the_topography)
   
## Solver Setup

The following headers specify that we use the Saint-Venant
solver together with [(dynamic) terrain
reconstruction](/src/terrain.h). The "inflow" library is a collection of functions for imposing
boundary conditions and the "hydro" library specify how to compute flow. 
We will use [inputs](/src/input.h) only when
restarting a simulation from a snapshot.
*/

#include "b-flood/saint-venant-topo.h"
#include "b-flood/inflow.h"
#include "b-flood/hydro.h"
#include "terrain.h"
#include "output.h"

/**
We then define a few useful constants.
*/
#define MAXLEVEL 8
#define MINLEVEL 6
#define HMAXE 1e-1

int nf,nc;
double Z0,Cf,tmax,etaleft,dx;


/**
## Parameters
    
The size of the area is exactly 5600 m. We take care to fix our
coordinate system as the one of the topography file. We fix the unity
of the gravity constant to $m.s^{-2}$, so unity of variables will be
meters for length and seconds for time. We fix the "dry" value at 1
mm. We also fix the total time of simulation at 24 hours.
*/

int main(){ 
  Z0 = 0.001;
  L0 = 5600.;
  X0 = 1523880;
  Y0 = 5430180;
  G = 9.81;
  Cf = 0.01;
  N = 1 << MAXLEVEL; // Which is equivalent to N = 2^(maxlevel)
  dry = 1e-3;
  tmax = 24*3600;  
  run();
}


/**
This function allows us to bring the river to the east only into its bed.
*/

int source(double x,double y) {
  if ( x >= X0 + L0 - 4*dx){
    if ( y >= 5430347 &&  y <= 5431922 )
      return 1;
  }
  return 0;
}



/**
## Initial condition

We initialize the topography by using the "terrain" function. 
We also include a RESTART option. This
option will be used if a crash occurs.
*/




event init (i = 0){
  terrain (zb, "./topo", NULL);

  dx = L0/N;
  DT = 100;
  foreach(){
    /**
     We first fix the water elevation to the altitude 0. 
     */
    h[] = max(0-zb[],0);
    
    /**
    The "river" scalar is used to fix the conditions at the edges. 
    It is equal to 1 where the flow can enter and zero where it is prohibited.
    We use the source function to define it properly.
    */
    if ( source(x,y) ){
    	river[] = 1;
  	}
  }
  boundary ({h,river});
  
  
  /**
  We save the RASTER files of the topography and the scalar "river". 
  */
  static FILE * grdz = fopen("./raster_river.asc","w");
  output_grd (river,grdz);
  
  grdz = fopen("./raster_zb.asc","w");
  output_grd (zb,grdz);
}


/**
## Bondaries conditions

On the right boundary, we impose an inflow equal to the hydrograph
data "flow.asc". In this file, the inflow is gaven in $m^3.s^{-1}$
(cumsec), so we have to fix a boundary condition as a total inflow
condition, which is not trivial. Thankfully, the function "discharge" do
this for us. This function takes as arguments the total inflow wanted,
the considered boundary and a scalar fixing the area where the inflow
occurs (the River[] scalar). "discharge" returns the water elevation
needed to impose the inflow.

By using the "discharge" function, you must fix a zero Neumann condition
on the velocity  :
*/

h[bottom] = 0;
u.n[right] = neumann(0);
u.t[right] = neumann(0);

/**
Now we can fix the water elevation on the right border with a
Dirichlet condition. The water elevation will change in respect to the
time, so we must do this in an event. Note that the "hydrograph" function
reads ASCI datas and returns a linear interpolation in respect to the
time. We fix the second argument of the function "hydrograph" to 3600
in order to convert hours in seconds.
*/

double eta1;
event discharge ( i++ ; t <= tmax ){

    double qimp = hydrograph ( "./flow.asc", 3600);
    eta1 = eta_b ( qimp , right, 1);

    
    h[right] = max ((river[] == 1 ? eta1 : zb[]) - zb[], 0);
    eta[right] = max ((river[] == 1 ? eta1 : zb[]), zb[]);


    // For the tides
    etaleft = hydrograph ( "./tide.asc" , 3600 );

   
    u.n[left] = - radiation( etaleft );
    
    /**
    We record the altitude of the ocean on the west boundary.
    */
    double etatide = 0;
    if(t >= 0 ){
      // We fix all scalars in this loop at the position X0 and YO+L0/2
      Point point = locate ( X0 , Y0 + L0/2. ); 
      etatide = h[] + zb[]; // We register the altitude of the water.
    }
    
    /**
    We record the altitude of water at the east boundary (karamea river).
    */
    double hr;
    if( t > 0 ){
    	Point point = locate ( X0 + L0 - dx ,5431450 );
    	hr = h[] + zb[];
    }
    /**
    We compute the real inflow entering the simulation by the east.
    */
    double qmes = compdischarge( xi = X0 + L0 - dx, yi = Y0, yf = Y0 + L0/2.);
    
    static FILE * fpeta = fopen("./check.dat","w");
    if(i == 0)
      fprintf(fpeta,"#1t 2etatide 3etaflow 4etariver 5qimp 6qmes 7hright \n");
    fprintf(fpeta,"%lf %lf %lf %lf %lf %lf %lf \n",t,etatide, etaleft,eta1,qimp,qmes,hr);
    

}

/**
## Friction term

Here we add the Graam friction term. We fix the roughness of the soil
equal to 1 mm. This term is treated implicitly.
*/

event friction (i++){
  foreach(){
    if (h[] > dry && Z0 > 1e-6){    
      double a = h[]/Z0;
      double r = a > 2.718 ? 2.5 * (log (a) - 1. + 1.359/a) : 0.46*a;
      double DU = norm(u)/sq(h[]*r);
  	     
      foreach_dimension()
         	u.x[] /= 1 + dt*DU;
	    
    }
  }
}

/**
## Adaptative refinement

Finding a good refinement criterium is none trivial. Here, we refine
according to the water elevation with a treshold of 5 cm. We have two
#if…#else branches selecting whether the simulation is being run on an
(adaptive) quadtree or a (static) Cartesian grid.
*/

int adapt_H() {
#if QUADTREE
  scalar h1[];
  foreach(){
    h1[] = h[];
  }
  boundary ({h1});

  astats s = adapt_wavelet ({h1}, (double[]){HMAXE} ,MAXLEVEL,MINLEVEL);
  if (s.nf != 0 || s.nc != 0){
    nc += s.nc;
    nf += s.nf;
  }   
  return s.nf;
#else // Cartesian
  return 0;
#endif
}
/**
   We refine quadtrees at each time step.
*/
event do_adapt_H (i++) adapt_H();


/**
## Outputs 
### Log output

Here we print in the terminal some information on the computation
process : timestep, maximum velocity, number of cells refined and
coarsened, etc.

*/
event logfile (t += 60) {

  stats s = statsf (h);
  scalar v[];
  foreach(){
    v[] = norm(u);
  }
  boundary({v});
  stats n =  statsf (v);
  if (i == 0)
    fprintf (stderr, "#t  i h.min h.max h.sum u.max dt\n");
  fprintf (stderr, "%.2lf %d %g %g %g %g %g\n", t/3600., i, s.min, s.max, s.sum, 
	   n.max, dt);
}

/**
### Movies

This is done every 60 seconds (t+=60). The static variable fp is NULL
when the simulation starts and is kept between calls (that is what
static means). The first time the event is called we set fp to a
ppm2mpeg pipe. This will convert the stream of PPM images into an mpeg
video.

We use the mask option of output_ppm() to mask out the dry
topography. Any part of the image for which m[] is negative (i.e. for
which etam[] < zb[]) will be masked out.
*/
event movies ( t += 60 ; t <= tmax) {
  if ( t > 0 ){
    scalar m[];
    foreach() {
      m[] = (h[] > dry) - 0.5;
    }
    boundary ({m});

    output_ppm (h, mask = m, min = 0, max = 5, n = 2^N, linear = true, file = "height.mp4");
    
    scalar v[];
    foreach()
      v[] = norm(u);
    boundary({v});
    output_ppm (v, mask = m, min = 0, max = 10, n = 2^N, linear = true, file = "vel.mp4");
    
    scalar fr = v;
    foreach(){
      if( h[] > dry )
      	  fr[] = norm(u)/sqrt(G*h[]);
      else 
      	  fr[] = 0;
    }

    output_ppm (fr, mask = m, min = 0, max = 1.1, n = 2^N, linear = true, file = "froude.mp4");
 
    scalar l = fr;
    foreach()
      l[] = level;
    output_ppm (l, min = MINLEVEL, max = MAXLEVEL, n = 2^N, file = "level.mp4");
 }
} 

/**
### Rasters

We take snapshots of simulation every hour. 
*/
event raster ( t += 3600  ) { 
  scalar m[];
  foreach() {
    m[] = (h[] > dry) - 0.5;
  }
  boundary ({m});
  
  double th = (t+1)/3600;
  char name[50];
  sprintf(name,"./raster-h-exp-t-%.0lf.asc",th);
  
  FILE * grdh = fopen(name,"w");
  output_grd (h,grdh, mask = m);
 
}



/**
 
## Compiling

We tell to the compiler to use the object kdt.o. This object is needed
to use the function "terrain". We also add the math header and
an optimisation option :

~~~bash 
./

~~~

## Running

Now we can run the simulation. 
   
~~~bash
./karamea.x

~~~
## Results


Check if the tidal condition was correctly treated :

~~~gnuplot
plot './check-exp.dat' u ($1/3600):2 w l
replot './check-exp.dat' u ($1/3600):3 w l
replot './tide.asc'
~~~

Check if the right inflow condition was correctly treated :

~~~gnuplot
plot './check-exp.dat' u ($1/3600):6 w l
replot './check-exp.dat' u ($1/3600):5 w l
replot './flow.asc'

~~~


## Link to the homepage
* [Homepage](Readme)


*/
