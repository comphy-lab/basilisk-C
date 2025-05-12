/** 
# Free time test case */

/**
## Include and parameters
This test uses the multilayer solver and the surface tension additive. */
#include "grid/multigrid1D.h"
#include "../hydro.h"
#include "layered/remap.h"
#include "../nh_PhiS.h"
#include "../solute.h"
#include "../surface.h"
#include "../tension.h"
#include "final_shape.h"

/**
## Geometry */
#define L 150.           // size of the box
#define T_END 6500
#define LEVEL 7          //horizontal resolution
#define layers 50        //vertical resolution


#define Y (2.e-3)        //alpha parameter
#define Rc (1/(2.*Y))    //Curvature ray (imposed by Y and hf)
#define theta 72.        //contact angle
#define Oh 0.1           //Ohnesorge number (1/nu)
#define d_sig (0.2)      //adimensionned Delta sigma
#define K (Y/d_sig)      //adsorption rate (adimensionned delta*lambda)
double sigma_i;          //Surface tension for an infinite height

#define activate 1       //Influence of convection on c

/**
## State laws
Adsorption at the surface is a function of the height and the molar fraction. 
Adsorbed quantity is computed and the solute is immediately diffused in the height.
*/
void adsorption (scalar c, scalar M, scalar h)
{
  foreach() {
    double total = (activate == 1 ? M[] : 0.);
    foreach_layer()
      total += (activate == 1 ? c[] : 1.) * h[];
    M[] = K * (total/eta[]) * (1 + K/eta[]);
    total -= M[];
    foreach_layer()
      c[] = total/eta[];
  }
  boundary({c,M});
}

/**
Surface tension evolution is linearized.*/

void tension (scalar sigma, scalar M)
{
  foreach()
    sigma[] = sigma_i + d_sig * (M[]/K - 1);
  boundary({sigma});
}

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
  nu = Oh;
  sigma_i = 1.;
  
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
  u.n[left] = neumann(0.);
  c[left] = neumann(0.);
  
  u.n[right] = neumann(0.);
//  phi[right] = dirichlet(- sigma_i/Rc);
  eta[right] = 2*eta[] - eta[-1] + sq(L0/N)/Rc;
//  eta[right] = contact_angle ((theta*pi/180.), L0/N);
  Neum_b = true;
  int ne = 0;

/** Initial conditions are fixed from the file final_shape.h, corresponding to the equilibrum.
Surface tension is obtained as a function of the thickness.
*/
  foreach(){
    double H = shape[(ne)][1];
//    double H = ne<128 ? shape[0][1] : shape[ne-128][1];
    ne++;    
    Hi[] = H;
    foreach_layer() {
      h[] = H/nl;
      c[] = 1.;
    }
    eta[] = H;
  }
  boundary({eta,c ,M}); 
  adsorption(c, M, h);
  tension(sigma, M);
  fp = fopen ("change", "w");
}

/** At each time step, the solute is adsorbed and imeditely diffused before the event remap (so as to have the new surface concentration for viscosity_event).
Surface tension is updated  as a function of the thickness. */
event remap (i++) {
  adsorption(c, M, h);
  tension(sigma, M);
}

/**
## Movie and outputs */

event logfile (i++; t <= T_END)
{
  /**
  At every timestep, we check the height convergence.*/
  double dH = change (eta, Hi)*N/L0;
  /**
  And we output the evolution of the maximum velocity. */
  scalar un[];
  foreach() {
    un[] = 0.;
    foreach_layer()
      un[] = max(un[],norm(u));
  }  
  fprintf (fp, "%g %g %g %g\n", t, normf(un).max, dH, statsf(eta).min);
}

/**
At the end, we save the final shape of the interface :*/ 
event final_shape (t = end) {
  FILE * fps = fopen ("final_shape", "w");
  FILE * fpss = fopen ("final_shape2", "w");
  fprintf(fps, "static double shape[][2] = {\n");
  foreach() {
    double U = 0;
    foreach_layer()
      U += u.x[];
    fprintf(fps, "{ %g, %g },\n", x, eta[]);
    fprintf(fpss, "%g %g %g %g %g %g %g %g\n", x, eta[], phi[1], sigma[], M[]/(K*(1 + K/eta[])), U/nl, u.x[0,0,0], u.x[0,0,nl-1]);
  }
  fprintf(fps, "};");
}

/**# Movie
To see the movement*/
event movie (i+=5)
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
	   i/5, t, X0, X0 + L0);
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

~~~gnuplot Height
set xlabel "t"
set ylabel "Hmin"
set xrange [6000:6500] 
plot 'change' u 1:4 w p ls 1
~~~
*/