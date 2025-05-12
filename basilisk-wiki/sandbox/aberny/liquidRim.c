/**
# Circular Hole Expansion

This is a 2D axisymmetric simulation of a liquid sheet with a hole in the 
middle. The idea is to evaluate the retracting speed of the liquid sheet. The 
simulation is axisymmetric, so we will not take into account the 3D effects.

We expect that the sheet will retract at a constant velocity, the Taylor-
Culick speed $u_c$. This velocity was first describe by Culick and Taylor. 
$u_c = \sqrt{2\sigma/\rho H}$.

We will setup a dimensionless simulation. All the parameters are equal to 1, 
except the viscosity which is equal to the Ohnesorge number 
$Oh=\mu/\sqrt{\rho\sigma H}$

*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

#define LEVEL 11
#define L0 40
#define R1 0.5
#define R0 5.

double tEnd = 10;

#define rhoInt 1.
#define rhoRatio 100.
#define muRatio 50.


/**
We allow the fluids to exit to domain from the different position, except the bottom part (where there is the axis of symmetry)*/

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

u.n[top] = neumann(0.);
p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);

int main(int argc, char *argv[]) {

  /**
  The domain will be $[-10, 10]\times[0:20]$ and will be resolved with
  $1024\times 1024$ grid points at the beginning.*/

  size (L0);
  origin (-L0/2., 0.);
  init_grid(1 << 9);

  /**
  By default the Ohnesorge number is equal to 0.1, but we can change that
  with an input argument.*/

  double Oh = 0.1;
  if (argc >=2) {
    Oh = atof(argv[1]); 
    tEnd = atof(argv[2]);
  }
  /**
  We define the parameters of the fluids. The internal fluid (the water) is 
  the fluid 1. The external fluid (the air) is the fluid 2.*/

  double muInt = Oh;

  rho1 = rhoInt, mu1 = muInt;
  rho2 = rhoInt/rhoRatio, mu2 = muInt/muRatio;
  f.sigma = 1.;
    
  run();
}

/**
We define the geometry with boolean operations. We define a
rectangular area by intersecting 3 lines. One is on the left (Ll), one on the 
right (Lr) and one on the bottom (Lb). We also define a circle.  Then, our 
geometry will be the union of the rectangle with the circle. Of course, we 
have to pay attention to the sign of the area, to be sure that the fluid 1 is 
where we want. */

double geometry (double x, double y) {
  double Lr = x - R1+L0/2.;
  double Lb = y - (R0+R1);
  double Circle = sq(x+L0/2.)+sq(y-(R0+R1))-sq(R1);
  double rectangleHalf = max(Lr, -Lb);
  return max(-Circle, -rectangleHalf);
}

event init (t = 0) {

  double iteration = 0;
  /**
  We initialise the geometry of the interface.*/  
  do {
  fraction(f, geometry(x,y));
  iteration++;
  } while (adapt_wavelet ({f}, (double[]){1e-3}, LEVEL).nf !=0 && 
    iteration<=10);
  static FILE * fp = fopen ("InitialState.ppm", "w");
  static FILE * fp2 = fopen("initial.gfs", "w");
  output_ppm(f, fp, min=0, max=1);
  output_gfs(fp2);
}

/**
We use an adaptive mesh with a maximum level of refinement of 10. We adapt 
the mesh with respect to the interface the velocity.*/

event adapt (i++) {
  double uemax = 0.1;
  double intemax = uemax/10;
  adapt_wavelet ({f,u}, (double[]){intemax,uemax,uemax,uemax}, LEVEL, 5);
}

/**
We output the tracer to have an idea of the evolution of the instability. */

event interface (i+= 20) {
  static FILE * fp = popen ("ppm2mpeg > liquidRim.mpeg", "w");
  static FILE * fp2 = popen ("ppm2mpeg > grid.mpeg", "w");
  output_ppm(f, fp, 256, min=0, max=1, box = {{-L0/2.,0},{-L0/2.+3,L0}});
  scalar l[];
  foreach()
    l[] = level;
  output_ppm(l,fp2, 256, min=5, max=LEVEL, box = {{-L0/2.,0},{-L0/2.+3,L0}});
}

double yPrev = -1, tPrev = -1, uYMax;

event extractPosition (i++) {

  /**
  We define a new vector field, h. */
  
  vector h[];

  /**
  We reconstruct the height function field and take the corresponding 
  minimum along the y axis. */
  
  heights (f, h);
  double yMin = +HUGE;;
  foreach()
    if (h.y[] != nodata) {
      double yi = y + height(h.y[])*Delta;
      if (yi < yMin) {
        yMin = yi;
        uYMax = u.y[];
      }
    }
  double xMin = +HUGE;;
  double xMax = -HUGE;;
  foreach()
    if (h.x[] != nodata) {
      double xi = x +height(h.x[])*Delta;
      if (xi < xMin)
        xMin = xi;
      
      if (xi > xMax)
        xMax = xi;
    }

  /**
  We also output the velocity of the end of the bulge.*/

  double veloTip = tPrev >= 0 ? (yMin - yPrev)/(t - tPrev) : 0.;
  double deltaT = t-tPrev;
  fprintf(stderr, "%g %g %g %g %g %g %g %d %d %d\n", t, deltaT, yMin, veloTip, 
    uYMax, xMin, xMax, mgp.i, mgpf.i, mgu.i);
  fflush (stderr);
  tPrev = t, yPrev = yMin;
}

/**
We output the interface of the fluid to track the evolution of the liquid rim. 
*/

event plotInterface (t += tEnd/10; t<= tEnd) {

  char name[80];
  sprintf (name, "interface-%f.txt", t);

  char gfs[80];
  sprintf (gfs, "output-%f.gfs", t);

  FILE* fp = fopen (name,"w");
  FILE* fp2= fopen (gfs, "w");

  output_facets (f, fp);
  output_gfs(fp2);
}

/**
We output, in the standard output file, the step with the corresponding
time. */

event logfile(i++) {
  printf ("i = %d t = %g\n", i,t);
  fflush(stdout);
}

// event gfs ( i++ ) {
//   char name[80];
//   if ((i>= 270 && i<= 280)|| (i>= 305 && i<= 315) || (i>= 975 && i<= 985) || (i>= 1065 && i<= 1080) || (i>= 2405 && i<=2420)) {
    
//     vector h[];
//     heights (f, h);

//     scalar hauteurY[];
//     scalar hauteurX[];
//     foreach() {
//       if (h.x[] != nodata){
//         hauteurX[] = height(h.x[])*Delta;
//       }
//       if (h.y[] != nodata) {
//         hauteurY[] = height(h.y[])*Delta;
//       }
//     }
//     sprintf (name, "output-%05ld.gfs",i);
//     FILE* fp = fopen (name, "w");
//     output_gfs (fp);
//   }
// }

event end (t = tEnd) {
  output_gfs (file = "final.gfs");
  static FILE * fp = fopen("Final.ppm", "w");
  output_ppm(f, fp, 512, min=0, max=1);
  dump (file = "dump");
}

/**
## Results

We expect the retraction velocity to be constant after a transient regime.

~~~gnuplot Evolution of the tip of the bulge
set xlabel 'time'
set ylabel 'y position'
f(x) = a*x + b
fit [4:] f(x) 'log' u 1:3 via a,b
plot [0:] 'log' u 1:3 every 10 t "Basilisk data",\
  f(x) t sprintf("%.2f t + %.2f", a, b)
~~~

We can also observe the evolution of the size of the liquid bulge.

~~~gnuplot Evolution of the size of the liquid bulge
set xlabel 'time'
set ylabel 'bulge size'
plot 'log' u 1:($7-$6) title "Bulge size"
~~~

If we supperpose the interface at the different time, we have:

~~~gnuplot Evolution of the interface
set xlabel 'x'
set ylabel 'y^*'
f(x) = a*x + b
fit [4:] f(x) 'log' u 1:3 via a,b
set size ratio -1
plot 'interface-0.000000.txt' w l title "t = 0", \
     'interface-1.000000.txt' u 1:($2 - a*1) w l t "t = 1", \
     'interface-2.000000.txt' u 1:($2 - a*2) w l t "t = 2", \
     'interface-3.000000.txt' u 1:($2 - a*3) w l t "t = 3", \
     'interface-4.000000.txt' u 1:($2 - a*4) w l t "t = 4", \
     'interface-5.000000.txt' u 1:($2 - a*5) w l t "t = 5", \
     'interface-6.000000.txt' u 1:($2 - a*6) w l t "t = 6", \
     'interface-7.000000.txt' u 1:($2 - a*7) w l t "t = 7", \
     'interface-8.000000.txt' u 1:($2 - a*8) w l t "t = 8", \
     'interface-9.000000.txt' u 1:($2 - a*9) w l t "t = 9", \
     'interface-10.000000.txt' u 1:($2 - a*10) w l t "t = 10"
~~~

Finally, we can also observe the tip velocity compare to the culick velocity, 
by using the viscous time as a time scale. This scale is: $\tau_{vis}=mu 
H/2\sigma$

~~~gnuplot Evolution of the tip velocity
set xlabel 't^*'
set ylabel 'u^*'
set yrange [0:1]
plot 'log' u ($1/(0.1/2)):(($4)/sqrt(2)) title "retraction velocity",\
  'log' u($1/(0.1/2)):(($5)/sqrt(2)) title "velocity in interface cell"
~~~
*/