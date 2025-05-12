/**
## IMPORTANT NOTE : 
We do not have the right to share the topography or plant cover that has been 
graciously loaned to us by the IGN. 
We do not have the right to share the rain RADAR rasters that has been 
graciously loaned to us by METEO-FRANCE. 
We share this script anyway for information purposes.
You will not be able to reproduce our results without buying these files.
*/

#include "b-flood/saint-venant-topo.h"
#include "b-flood/manning.h"
#include "b-flood/rain.h"
#include "b-flood/infiltration.h"
#include "b-flood/vthreshold.h"
#include "b-flood/inflow.h"
#include "terrain.h"
#include "input.h"

#define SPELEVEL 10
#define MAXLEVEL 9 // dx = 13.6 m
#define MINLEVEL 6 // dx = 235 m
#define EMAXE 0.05

double maxrain = 80 /3.6*1e-6;
int nf,nc,ncell;
timer t0;

int main() {
  threshold = 10;	
  L0 = 7000.00;
  X0 = 1022400;
  Y0 = 6278404;
  G = 9.81;
  dry = 1e-4;
  N = 1 << MAXLEVEL; // Which is equivalent to N = 2^(maxlevel)
  run();
}
/**
## Stop 
*/
event stopprog( t = 5*3600 ) // 18h00 -> 23h00
  fprintf(stderr,"stop\n rainmax = %lf",maxrain*1e3);


int spec(double x, double y){
  double R = 30;
  double x0 = 1024228, y0 = 6280726; // Cannes Town Hall
  if( sq(x - x0) + sq(y - y0 ) <= sq(R) )
    return 1;
  
  x0 = 1024023, y0 = 6280923; // Fire station
  if( sq(x - x0) + sq(y - y0 ) <= sq(R) )
    return 1;

  x0 = 1024252, y0 = 6280890; // Other Town Hall
  if( sq(x - x0) + sq(y - y0 ) <= sq(R) )
    return 1;

  x0 = 1024530, y0 = 6281037; // Central police station
  if( sq(x - x0) + sq(y - y0 ) <= sq(R) )
    return 1;

  x0 = 1024947, y0 = 6281920; // police station
  if( sq(x - x0) + sq(y - y0 ) <= sq(R) )
    return 1;

  x0 = 1023557, y0 = 6282754; // Other Town Hall
  if( sq(x - x0) + sq(y - y0 ) <= sq(R) )
    return 1;

  x0 = 1024628, y0 = 6283539; // Cannet Town Hall
  if( sq(x - x0) + sq(y - y0 ) <= sq(R) )
    return 1;

  x0 = 1024399, y0 = 6283146; // Other Cannet Town Hall
  if( sq(x - x0) + sq(y - y0 ) <= sq(R) )
    return 1;

  return 0;
}

/**
## Initial condition
*/
scalar veget[];
event init (i = 0){
	
  terrain (zb, "./TopoW/topo", NULL);
 
  
  DT = 60;
  foreach(){
    //We fill the mediteranean sea
    eta[] = 0;
    h[] = max(eta[] - zb[], 0);
  }
  
  // the vegetation
  input_grd(veget, file = "./Topo/veget-good.asc");   
  foreach(){
    // high vegetation zone
    if( veget[] > 0.1 ){
      nmanning[] = 0.1;
      // Loamy sand
      Ch[] = 0.061;
      K[] = 3e-5/3.6; //cm/h in m/s
      Sat[] = 0.05;
    }
      
    // roads and city
     else if( zb[] < 50 && zb[] >= 0 ){
      nmanning[] = 0.03; 
      
      //Sandy Clay
      Ch[] = 0.21;
      K[] = 0.15e-5/3.6;
      Sat[] = 0.05; //5% free
    }
    
    // in the sea
    else if( zb[] < 0 ){
      nmanning[] = 0.03; 
      
      //fully saturated 
      Ch[] = 0;
      K[] = 0;
      Sat[] = 0; //0% free
    }
    else{
      nmanning[] = 0.06; // other
      //Sandy Clay
      Ch[] = 0.21;
      K[] = 0.15e-5/3.6;
      Sat[] = 0.05; //5% free
    }
  }

  refine( level < SPELEVEL && spec(x,y));


  FILE * fpp = fopen("raster_manning.asc","w");
  output_grd(nmanning, fpp);

  fpp = fopen("raster_zb.asc","w");
  output_grd(zb,fpp);
      
  fpp = fopen("raster_K.asc","w");
  output_grd(K, fpp);
    
  t0 = timer_start();
}

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
event do_adapt_H (i++) adapt_H();

/**
## Boundary conditions
*/

u.n[top] = neumann(0);

u.n[left] = zb[] < 0 ? - radiation(0) : neumann(0);
u.n[right]  = + radiation(0);
u.n[bottom] = - radiation(0);

/**
## Rain

*/
scalar train[];
double maxrain;
event comprain( i++ ){
  foreach(){

    train[] += rain[] * dt;
    if( train[] > maxrain )
      maxrain = train[]; 
  }
}
  
event rrain ( t+=60){
  read_lameeau(name = "./ASC/LAME_EAU.COLL.20151003", start = 18, end = 23, step = 60, linear = true, timelin = true);
}


//////// OUTPUT
/**
## Logfile
*/
event logfile (t += 72) {
  timing tr = timer_timing (t0, i, 0, NULL);
  stats s = statsf (h);
  scalar v[];
  foreach()
    v[] = norm(u);
  stats no = statsf (v);
  if( i == 0 )
    fprintf (stderr, "#1t 2treal 3i 4h.min 5h.max 6h.sum 7u.min 8u.max 9dt\n");  
  fprintf (stderr, "%g %g %d %g %g %g %g %g %g\n", t/3600.,tr.real, i, s.min, s.max, s.sum, 
	   no.min, no.max, dt);
}

/**
## Flood print
We store the maximum level of water in the scalar hmax[].
 */

scalar hmax[];
scalar alea[];
scalar u2[], u2max[];
event comphmax(i += 10) {
  foreach() {
    if (h[] > hmax[]) 
      hmax[] =  h[];
  }
  
  foreach(){
    double vel = norm(u);
    if( (h[] >= 1 && vel >= 0.2) || vel >= 1 )
      alea[] = 4;
    else if((h[] >= 1 || (h[] >= 0.5 && vel >= 0.5)) && (alea[] < 4))
      alea[] = 3;
    else if(( h[] >= 0.3 || vel >= 0.2) && alea[] < 3)
      alea[] = 2;
    else if ( alea[] < 2 )
      alea[] = 1;
    
    u2[] = 1000*h[]*sq(vel);
    if( u2[]  > u2max[] )
      u2max[] = u2[];
  }
}

/**
## Rasters
*/

event raster( t += 1800  ) { 
  
  double seuil = 0.1; 
  scalar m[];
  foreach() {
    m[] = (h[] > dry + seuil) - 0.5;
  }
  boundary ({m});
  
  double th = (t+1)/3600;
  char name[50];
  sprintf(name,"Raster-t-h-%.1lf.asc",th);


  FILE * grdh = fopen(name,"w");
  output_grd (h,grdh, mask = m);
}
/**
## Movies

We record movies of the water depth, the velocity and the level of
refinement.*/

event movies ( t += 30 ) {
  
    scalar m[], v[];
    
    foreach(){
      v[] = norm(u);
      m[] = (h[] > dry ) - 0.5;
      if( zb[] < 0 ) m[] = -0.5;
    }    
    boundary ({m,v});


    output_ppm (h, mask = m, min = 0, max = 2, n = N, linear = true, file = "height.mp4");

    output_ppm (Vol, mask = m, min = 0, max = 2, n = N, linear = true, file = "VolInf.mp4");

    output_ppm (v, mask = m, min = 0, max = 8, n = N, linear = true, file = "vel.mp4");

    output_ppm (u2, mask = m, n = N, linear = true, file = "u2.mp4");
    

    scalar l = m;
    foreach()
      l[] = level;
    output_ppm (l, min = MINLEVEL, max = SPELEVEL, n = N, file = "level.mp4");

    output_ppm (rain, min = 0, max = 120e-3/3600., n = N, linear = true, file = "rain.mp4");
    
    output_ppm (train, min = 0, max = 200e-3, n = N, linear = true, file = "totalrain.mp4");
    
  
}

/**
## Fractal dim 

*/
int ifractal = 0;
event fracdim( t += 60){
  int ncell = 0;
  foreach(reduction(+:ncell))
    ncell++;
  static FILE * ffrac = fopen("fdim.dat","w");
  if( ifractal == 0){
    fprintf(ffrac,"# t \t i \t Ncell \t Nrefined \t Ncoarsen");
    ifractal++;
  }
  fprintf(ffrac,"%.1lf \t %i \t %i \t %i \t %i\n",t/3600.,i,ncell,nf,nc);
  nf = nc = 0;
  }

event ending(t = end){
  static FILE * rptrain = fopen ("raster_train.asc", "w");
  output_grd (train, rptrain);

  static FILE * rpvolinf = fopen ("raster_volinf.asc", "w");
  output_grd (Vol, rpvolinf);

  static FILE * rphmax = fopen ("raster_hmax.asc", "w");
  output_grd (hmax, rphmax);

  static FILE * rpalea = fopen ("raster_alea.asc", "w");
  output_grd (alea, rpalea);

  static FILE * rpu2max = fopen ("raster_u2max.asc", "w");
  output_grd (u2max, rpu2max);

}
/**
## Link to the homepage
* [Homepage](Readme)
*/