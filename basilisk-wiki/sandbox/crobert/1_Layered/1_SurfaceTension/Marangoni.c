/**  
# Marangoni flow vs gravity

We test here the ability of the multilayer solver to reproduce correctly the Marangoni flow given a surface tension gradient. 
This example is inspired of Guyon, Hulin, Petit (4.8.3).

A tank of finite dimensions containing a liquid film of mean height $h_0$ is considered. 
A gradient of surface tension $\sigma = 1 + b*x$, $b>0$ is imposed on the surface of the liquid (phisically, such a gradient could be obtained for a thin film with a temperature gradient).

Marangoni flows tend to accumulate liquid on the left side of the tank whereas gravity tends to flatten the surface, which trigger a reversed flow at the bottom of the layer. The competition between those two forces can be analitically studied and, at the equilibrium, the interface is tilted.
For unidirectional flow and small deformations, we get :

$$
\nu \frac{\partial v_x}{\partial y^2} = g \frac{\mathrm{d}h}{\mathrm{d}x}
$$

After integration with equilibrium condition other a vertical slice and the non-slip condition, we get :

$$
v_x = \frac{g}{\nu} \frac{\mathrm{d}h}{\mathrm{d}x} \left( \frac{y^2}{2} - \frac{h_0 y}{3}\right)
$$

and

$$
\frac{\mathrm{d}h}{\mathrm{d}x} = \frac{3}{2} \frac{b}{\rho g h_0}
$$

![Height](Marangoni/movie.mp4)

*/


/**
## Include and parameters
This test uses the multilayer solver and the surface tension additive. */
#include "grid/multigrid1D.h"
#include "../hydro.h"
#include "layered/remap.h"
#include "../nh_PhiS.h"
#include "../tension.h"

#define _g 100.
#define Re 1.

#define sigma0 10.
double h_0 = 0.4; //Averaged height
#define b 0.1  //Surface tension gradient

#define L 8
#define T_END 5.
#define DELTA_T (T_END/200.)
#define H_ERR 1e-10

/**
## Main function
The boundary conditions are set to dirichlet conditions. The test is done with viscosity and gravity.
However, Laplace pressure is killed to zero in order to avoid boundary effects (by defaults the contact angle would be 90). */
int test_number = 0;
int main()
{
  L0 = L;
  origin(-L0/2);
  G = _g;
  nu = 1/Re;
  N = 64;
  nl = 16;
  Laplace = 0;
  
  
  run();
  test_number++;

  for (h_0 = 0.2; h_0 <= 0.8; h_0 += 0.1){
    run();
    test_number++;
  }  
  h_0 = 0.4;
  
  test_number = 100;
  for (nl = 1; nl <= 8; nl *= 2){
    run();
  }
  nl = 16;
}

/**
## Initialisation 
The initial shape of the interface is a slightly tilted, corresponding to the theoretical equilibrium shape. 
Surface tension is initialised with the surface tension gradient b. */
scalar Hi[];
double maxerror;
event init (i = 0)
{
  u.n[right] = dirichlet(0);
  u.n[left] = dirichlet(0);
  foreach() {
    double H = h_0 + 1.5*b/(G*h_0)*x;
    Hi[] = H;
    foreach_layer() {
      h[] = H/nl;
      sigma[] = sigma0 + b*x;
    }
  }
  boundary({h, sigma, Hi});
  maxerror = 0;
}

/**
## Outputs 
*/

event logfile (t += DELTA_T; t = 0; t <= T_END)
{
  /**
  At every timestep, we check whether the height field has
  converged. */
  double dH = change (eta, Hi)*N/L0;
  if (i > 1 && dH < H_ERR)
    return 1; /* stop */
  maxerror = max(maxerror, dH);
  
  char name[80];
  sprintf (name, "height-%g", h_0);
  static FILE * fph = fopen (name, "w");
  if (t > DELTA_T)
    fprintf(fph, "%g %g %g\n", t, statsf(eta).sum, dH);
}

event error (t = end) {
  /**
  We recompute the reference solution. */
  scalar Href[];
  foreach()
    Href[] =  h_0 + 1.5*b/(G*h_0)*x;
  boundary({Href});
  
  /**
  We compute the shape error *ec* and output these in the *log* file. */
  scalar ef[];
  foreach() {
    ef[] = (eta[] - Href[]);
  }
  norm ne = normf (ef);
  fprintf (stderr, "%g %g %g\n", ne.avg, ne.rms, ne.max);

  /**
  We save the velocity profile at the center of the domain for the different nl. */
  if (test_number > 99 || test_number == 0) {
    char name[80];
    sprintf (name, "center-%d", nl);
    static FILE * fpc = NULL;
    fpc = fopen (name, "w");
    int J = 0;
    foreach() {
      double z = zb[];
      if (J == N/2)
        foreach_layer() {
          z += h[]/2;
  	      fprintf (fpc, "%g %g %g\n", z, u.x[], 1.5*b/nu*(sq(z)/(2.*h_0) - z/3.));
          z += h[]/2;
        }
      J += 1;
    }
    fflush(fpc);
    fclose(fpc);
  }
  else {}
    /**
    At the end, we save the final shape of the interface :*/ 
    char nameBis[80];
    sprintf (nameBis, "final_shape-%g", h_0);
    FILE * fpf = NULL;
    fpf = fopen (nameBis, "w");
    foreach() {
      fprintf(fpf, "%g %g %g\n", x, eta[], Href[]);
    }
    fflush(fpf);
    fclose(fpf);
    
  static FILE * fpe = fopen ("error", "w");
  fprintf(fpe, "%d %d %g %g %g %g %g %g\n", N, nl, h_0, b, G, nu, maxerror, ne.max);

}

/**
# Results
~~~gnuplot Final shapes using multilayers

darkgray="#666666"
purple="#9E5DB6"
blue="#5082DC"
turquoise="#008C7D"
forest="#149632"
orange="#FF780F"
raspberry="#FA0F50"
set style line 1 pt 7 ps 0.7

set terminal @SVG enhanced size 640,640 font ",8"
set output '_final_shape.svg'
set key 
set border
set tics
set xlabel "x"
set ylabel "y"
set size ratio 0.2

plot \
  'final_shape-0.4' u 1:2 w l lc rgb raspberry t 'multilayer', \
  'final_shape-0.4' u 1:3 w l lc rgb forest t 'model'
~~~

~~~gnuplot Error on the final shapes with the height functions
set output '_delta_final_shape-0.svg'
set xlabel "x"
set ylabel "y"
set size ratio 1

plot 'final_shape-0.2' u 1:($2 - $3) w l t 'h_0 = 0.2' , \
  'final_shape-0.4' u 1:($2 - $3) w l t 'h_0 = 0.4' , \
  'final_shape-0.6' u 1:($2 - $3) w l t 'h_0 = 0.6' , \
  'final_shape-0.8' u 1:($2 - $3) w l t 'h_0 = 0.8'
~~~

~~~gnuplot multi changes
set output '_multi_changes.svg'
set xlabel "t"
set ylabel "dH"
plot [0:1] 'height-0.2' u 1:3 w p t 'h_0 = 0.2', \
  'height-0.4' u 1:3 w p t 'h_0 = 0.4', \
  'height-0.6' u 1:3 w p t 'h_0 = 0.6', \
  'height-0.8' u 1:3 w p t 'h_0 = 0.8'
~~~

~~~gnuplot Numerical and analytical velocity profiles at the center depending of the number of layers.
set output '_velocity.svg'
set xlabel "u"
set ylabel "z"
set key left top
plot 'center-1' u 2:1 w l t 'nl = 1', \
  'center-2' u 2:1 w l t 'nl = 2', \
  'center-4' u 2:1 w l t 'nl = 4', \
  'center-8' u 2:1 w l t 'nl = 8', \
  'center-16' u 2:1 w l t 'nl = 16', \
  'center-16' u 3:1 w l lc rgb forest t 'model'
~~~
When the number of layers increase, the velocity approaches very well the analytical velocity.
*/