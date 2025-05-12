/** 
# Free time test case */

/**
## Include and parameters
This test uses the multilayer solver and the surface tension additive. */
#include "grid/multigrid1D.h"
#include "../hydro.h"
#include "layered/remap.h"
#include "../nh_PhiS.h"
#include "../tension.h"
#include "final_shape.h"

/**
## Geometry */
#define L 60.       // size of the box
#define T_END 1000.
#define H_ERR 1e-10
#define LEVEL 7      //horizontal resolution
#define layers 1     //vertical resolution

#define alpha 0.01   //alpha parameter
#define Rc 50.       //Curvature ray
#define theta 60.    //contact angle
#define We 1.        //Weber number (1/sigma)
#define Re 100.      //Reynolds number (1/nu)
double sigma_inf;
 
/**
## Main function */
int main()
{
  L0 = L;
  origin (0.);
  N = 1 << LEVEL;
  nl = layers;
  gradient = NULL ;
  G = 0;
  nu = 1/Re;
  sigma_inf = 1/We;
  
  run();
}

/**
## Initialisation */
FILE * fp = NULL;
scalar Hi[];
event init (i = 0)
{
/** We set the boundary conditions. Symetric conditions are imposed on the left : a dirichlet zero-condition on the normal velocity and a neumann zero-derivative for the height. On the right, we impose a neumann condition on velocity, a pressure at a given curvature and a contact angle.
*/
  eta[left] = neumann(0.);
  u.n[left] = dirichlet(0.);
  
  u.n[right] = neumann(0.);
  phi[right] = dirichlet(- sigma_inf/Rc);
  eta[right] = contact_angle ((theta*pi/180.), L0/N);
  int ne = 0;
  
/** Initial conditions are fixed from the file final_shape.h, corresponding to the equilibrum.
Surface tension is obtained as a function of the thickness.
*/
  foreach(){
    double H = shape[ne++][1];
    Hi[] = H;
    foreach_layer()
      h[] = H/nl;
    eta[] = H;
    sigma[] = sigma_inf * (1 + alpha/(2*eta[]));
  }
  boundary({eta, sigma}); 
  fp = fopen ("change", "w");
}

/** At each time step, surface tension is updated  as a function of the thickness. */
event pressure (i++)
{
  foreach() 
    sigma[] = sigma_inf * (1 + alpha/(2*eta[]));
  boundary({sigma});
}

/**
## Movie and outputs */

event logfile (i++; t <= T_END)
{
  /**
  At every timestep, we check whether the volume fraction field has
  converged. */
  double dH = change (eta, Hi)*N/L0;
  if (i > 1 && dH < H_ERR)
    return 1; /* stop */

  /**
  And we output the evolution of the maximum velocity. */
  scalar un[];
  foreach() {
    un[] = 0.;
    foreach_layer()
      un[] = max(un[],norm(u));
  }  
  fprintf (fp, "%g %g %g\n", t, normf(un).max, dH);
}

/**
At the end, we save the final shape of the interface :*/ 
event final_shape (t = end) {
  FILE * fps = fopen ("final_shape", "w");
  FILE * fpss = fopen ("final_shape2", "w");
  fprintf(fps, "static double shape[][2] = {\n");
  foreach() {
    fprintf(fps, "{ %g, %g },\n", x, eta[]);
    fprintf(fpss, "%g %g %g\n", x, eta[], phi[1]);
  }
  fprintf(fps, "};");
}

/**# Movie
To see the movement*/
event movie (i+=10)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fprintf (fp, "set term pngcairo font ',9' size 800,250;"
	     "set size ratio 0.3\n");  
  fprintf (fp,
	   "set output 'plot%04d.png'\n"
	   "set title 't = %.2f'\n"
	   "p [%g:%g][-8:8]'-' u 1:2:3 w filledcu lc 3 t '',"
     "'./final_shape' u 1:2 w l",
	   i/10, t, X0, X0 + L0);
  fprintf (fp, "\n");
  foreach_leaf() {
    double H = 0.;
    foreach_layer()
      H += h[];
    fprintf (fp, "%g %g %g", x, H, -H);
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  fflush (fp);
}

event moviemaker (t = end)
{
  system ("rm -f movie.mp4 && "
	  "ffmpeg -r 25 -f image2 -i plot%04d.png "
	  "-vcodec libx264 -vf format=yuv420p -movflags +faststart "
	  "movie.mp4 2> /dev/null && "
	  "rm -f plot*.png");
}

/**
# Results
~~~gnuplot Final shape
set terminal svg enhanced size 640,640 font ",8"
set xlabel "z"
set ylabel "H"
set size ratio 0.3
plot 'final_shape2' u 1:2 w p ls 1
~~~

~~~gnuplot Change
set xlabel "t"
set ylabel "dh"
set size ratio 1
set logscale y
plot [0:500]'change' u 1:3 w p ls 1
~~~
*/
