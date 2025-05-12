/**

## Numerical Setup 

The Navier-Stokes flow solver is used in three dimensions. Further, a
switch is placed for 3D visualization (see below). 

![The ring breaks up.](crosswind/l2pi210.mp4)

~~~gnuplot Traject in the x-y plane
set xlabel 'x/R_o'
set ylabel 'y/R_o'
plot 'datapi210' u 2:3:1 palette  
~~~
*/
#include "grid/octree.h" 
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h" //Trace performance
		      
#define BVIEW 1 //0 is off.
#if BVIEW
#include "view.h"
#endif

double Uo = 1, Ro = 1;                // Normalized values
double Pi1 = 6, Re = 5000, Pi2 = 10;  // Dimensionless numbers  
double to, nu, Uc;                    // Dependend variables

#define RADIUS (sqrt(sq(y) + sq(z)))

u.n[left] = dirichlet ((t <= to)*(RADIUS <= Ro)*Uo); // Injection

u.n[right]  = neumann (0);   // free Outflow
p[right]    = dirichlet (0);

int maxlevel = 8;
int main() {
  L0 = 15*Ro;
  X0 = Y0 = Z0 = -L0/2;
  nu = Uo*Ro/Re;
  to = Pi1*Ro/Uo;
  const face vector muc[] = {nu, nu, nu};
  mu = muc;
  Uc = Uo/Pi2;
  periodic (top);
  run();
}


event init (t = 0) {
  refine (x < X0 + Ro    && RADIUS < 2*Ro   && level < maxlevel - 1);
  refine (x < X0 + Ro/2. && RADIUS < 1.1*Ro && level < maxlevel);
  foreach()
    u.y[] = Uc;
  boundary({u.y});
}


event movie (t += 0.1) {
  scalar omg[], lev[];
  vorticity (u, omg);
  boundary ({omg});
  foreach()
    lev[] = level;
  output_ppm (omg, file = "o.mp4", n = 300, linear = true, min = -1, max = 1);
  output_ppm (lev, file = "l.mp4", n = 300, min = 1, max = maxlevel);
}


#include "lambda2.h"
double val = -1;

#if BVIEW
event bviewer (t += 0.1) {
  scalar l2[];
  lambda2 (u, l2);
  boundary ({l2});
  view (theta = 0.4, phi = 0.3);
  isosurface ("l2", val);
  translate (z = Z0)
    cells();
  char video[99];
  sprintf (video, "l2pi2%g.mp4", Pi2);
  save (video);
}
#endif


event track_and_trace (t += 0.5) {
  scalar l2[];
  lambda2 (u, l2);
  double xp = 0, yp = 0, zp = 0, tot = 0;
  foreach (reduction (+:xp) reduction (+:yp) reduction (+:zp) reduction (+:tot)) {
    if (l2[] < val) {
      tot +=   l2[]*dv();
      xp  += x*l2[]*dv();
      yp  += y*l2[]*dv();
      zp  += z*l2[]*dv();
    }
  }
  if (tot != 0) {
    char file[99];
    sprintf (file, "datapi2%g", Pi2);
    static FILE * fp = fopen (file, "w");
    fprintf (fp, "%g\t%g\t%g\t%g\n", t, xp/tot, yp/tot, zp/tot);
  }
}


event adapt (i++) {
  double uc = Uo*0.04;
  adapt_wavelet ((scalar*){u}, (double[]){uc, uc, uc}, maxlevel);
}


event stop (t = 6*to);
