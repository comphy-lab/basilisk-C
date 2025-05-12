
/**
# **Wave**
## Set up for the code 
* Download [bash.sh](http://basilisk.fr/sandbox/aaykin/wave_code/bash.sh) and wave.c
* Run:. */

/**
## Code

### Variables and Libraries used */

/**
# **Wave** 
## Set up for the code 
* Download [bash.sh](http://basilisk.fr/sandbox/aaykin/wave_code/bash.sh) and wave.c
* Run:. */

/**
## Code

### Libraries used */

#define zoom 1
#define P 0
#define multiple_run 1

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/check_eta.h"
#if P
#include "layered/perfs.h"
#endif




/**
### Initialisation of variables, main part of the code and exit condition
*/
double amp=0.01,k=1.0,H=1.0,om;
int t_end=10,counter=0,POINT,LEVEL=4,il=1;




event init (i = 0) {
 
  G=9.81;
  om = sqrt(G*k*tanh(k*H));
  foreach(serial)
  {
    foreach_layer()
      {
	h[] = (H - ((amp*k)/om)*sin(k*(x-pi/2))*sinh(k*H))/nl;
	u.x[] = -(amp*k)*sin(k*(x-pi/2))*cosh(k*H);
      }
  }
}

int main() 
{
  G=9.81;
  om = sqrt(G*k*tanh(k*H));
  periodic(right);


#if zoom
  L0 = 2*pi/k;
#else
  L0=(2*pi/k)*10;
#endif
  for (nl = 1.0;nl <= 1; nl += 1)
    {
      printf("%i",nl);
      printf("\n");
#if multiple_run
      for (POINT = 2;POINT <= pow(2,LEVEL); POINT *= 2)
	{
	  init_grid(POINT); 
	  printf("%i",POINT);
	  run();
	}
#else
      POINT = pow(2,LEVEL);
      init_grid(POINT); 
      FILE * fp = fopen ("error.txt", "a");
      fprintf(fp,"\n\n");
      printf("%i",POINT);
      run();  
#endif
    }
}

/**
### Code's events and functions
*/

/**
* function z0
 */

double z0 (int q)
{
  return (H - ((amp*k)/om)*sin(-t_end*om+k*q)*sinh(k*H));
}

/**
* output's events
 */


event error (t = t_end)
{	 
  static FILE * fp = fopen ("error.txt", "a");
  double err_avg=0.0,xmax=0.0;
  
  foreach(serial)
  { 
    foreach_layer()
      {
	err_avg += fabs(eta[]-z0(x));
	fprintf (fp,"%i %g %i %g %g\n",nl,t,POINT,x, fabs(eta[]-z0(x)));
      }
  }
  if(il!=nl)
    {
      fprintf(fp,"\n\n");
      il+=1;
      printf("%i",il);
    }

  fprintf(fp,"\n\n");
}
event end (t = t_end){}


event plot (t<=t_end;t+=1.0) {
  static FILE * fp = popen ("gnuplot -persist 2> /dev/null", "w");
  fprintf (fp,"set term gif animate size 1920,1080;set output 'surface/surface_t_%.2lf.png';set size ratio 0.5\n",t);
  fprintf (fp,"\nset grid\n");
  fprintf (fp,"set title 'Wave at t= %.2lf '\n"
	   "G = 9.81 \n"
	   "k = 1.0 \n"
	   "H = 1.0 \n"
	   "om = sqrt(G*k*tanh(k*H)) \n"	 
	   "amp = 0.01 \n"	  
	   "set xrange[0:L0] \n"
	   "set samples 1000 \n"
	   "set yrange[0.5:1.5] \n"
	   "z(x) = H - ((amp*k)/om)*sin(-t*om+k*(x-pi/2))*sinh(k*H); t= %.2lf ;\n"
	   //"z(x) = H - cos(x); t= %.2lf ;\n"
	   "p[%g:%g]  '-' u 1:2 t'Numerical result' w lp ,"
	   "z(x) t 'Analytical result'\n",
	   t,t,X0,X0+L0);
  foreach(serial)
  {
    fprintf (fp,"%g %g %g\n", x,  eta[], t);
    fprintf (fp,"e\n\n");
    fflush (fp); 
  }
}

/**
* Animation of the surface analytically and numerically according to time
 */

event plot (t<t_end;t+=0.016) {
  static FILE * fp = popen ("gnuplot -persist 2> /dev/null", "w");
#ifdef gnuX
#else
  if(t==0) fprintf (fp,"set term gif animate size 1920,1080;set output 'wave.gif';set size ratio 0.5\n");
#endif
  fprintf (fp,"\nset grid\n");
  fprintf (fp,"set title 'Wave at t= %.2lf '\n"
	   "G = 9.81 \n"
	   "k = 1.0 \n"
	   "H = 1.0 \n"
	   "om = sqrt(G*k*tanh(k*H)) \n"	 
	   "amp = 0.01 \n"	  
	   "set xrange[0:L0] \n"
	   "set yrange[0.5:1.5] \n"
	   "set samples 1000 \n"
	   "z(x) = H - ((amp*k)/om)*sin(-t*om+k*(x-pi/2))*sinh(k*H); t= %.2lf ;\n"
	   //"z(x) = H - cos(x); t= %.2lf ;\n"
	   "p[%g:%g]  '-' u 1:2 t'Numerical result' w lp ,"
	   "z(x) t 'Analytical result'\n",
	   t,t,X0,X0+L0);
  foreach(serial)
    fprintf (fp,"%g %g %g\n", x,  eta[], t);
  fprintf (fp,"e\n\n");
}

/**
   Plotting of the wave a fixed time

   ~~~gnuplot Wave at fixed time
   #!/bin/bash
   set terminal png size 1920,1080
   set output 'surface.png'
   G = 10
   k = 1.0
   H = 1.0
   om = sqrt(G*k*tanh(k*H))
   t = '10.0'
   n = t*100
   amp = 0.001
   z(x) = H - ((amp*k)/om)*sin(-t*om+k*x)*sinh(k*H)
   set xlabel 'x'
   set ylabel 'Z0'
   set grid
   set zrange[-2:2]
   set samples 1000
   plot for[in=n:n] for [H in t] "surface.txt" i in u 1:2 w lp t 'Surface time ='.H, z(x) w l t 'Analytical solution'
   ~~~
*/

/**
   Plotting of the error according to the level

   ~~~gnuplot Wave at fixed time
   #!/bin/bash
   set terminal png size 1920,1080
   set output 'error.png'
   set xlabel 'Point'
   set ylabel 'Error'
   set xlogscale
   set grid
   plot "error.txt" u 1:2 w lp t 'Error'
   ~~~
*/