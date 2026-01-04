/**
#  I want to check if the slip implementation works with non flat surfaces. 

I copy the poiseuille45.c test case and adapt it as a coutte flow with navier slip boundary condition as in sandbox/twitkamp/MovingEmbed/embed_slip.c
*/

#include "embed_navier.h"
#include "navier-stokes/centered.h"
#include "view.h"

#define WIDTH 0.5
#define EPS 1e-14

u.y[embed] = dirichlet(((y < x - 1e-14 + 0.75 && y > x - 1e-14 + 0.25 ) ||  (y < x - 1e-14 - 0.25 && y > x - 1e-14 - 0.75 ) ) ? sqrt(2)/2 : 0);
u.x[embed] = dirichlet(((y < x - 1e-14 + 0.75 && y > x - 1e-14 + 0.25 ) ||  (y < x - 1e-14 - 0.25 && y > x - 1e-14 - 0.75 ) )  ? sqrt(2)/2 : 0);

double slip;
int main()
{
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (1. [0]);
  DT = HUGE [0];
  
  origin (-0.5, -0.5);
  periodic (right);
  periodic (top);
  N = 64;
  stokes = true;
  TOLERANCE = 1e-7;
	
  for (double s = 0.2; s <= .2; s += 0.2){
    slip = s;
    run();
  }
}

scalar un[];
scalar lambdax[], lambday[];
event init (t = 0) {
  
  mu = fm;

  /**
  The channel geometry is defined using Constructive Solid Geometry. */  

  solid (cs, fs, difference (union (y - x - EPS, x - y - 0.5 + EPS),
			     y - x - 0.5 + EPS));

  /**
  We use "third-order" [face flux interpolation](/src/embed.h). Not sure why needed yet.*/

  for (scalar s in {u})
    s.third = true;

  /**
  We initialize the reference velocity. */
  

  foreach(){
    lambdax[] = (y > x - 0.25 ? (y < x + 0.25 ? slip : 0.) : 0.);
    lambday[] = (y > x - 0.25? (y < x + 0.25 ? slip : 0.) : 0.);
  }

  u.x.lambda = lambdax;
  u.y.lambda = lambday;

  foreach()
    un[] = u.y[];

}

/**
We check for a stationary solution. */

event logfile (t += 0.1; i <= 1000) {
  double du = change (u.y, un);
  if (i > 0 && du < 1e-6)
    return 1; /* stop */
}

/**
We compute the tangential velocity and display the solution using bview. */
scalar utan[], unorm[];
scalar e[], theoretical[];
event profile (t = end) {
  printf ("\n");
	foreach(){
    utan[] = (u.x[] + u.y[])/sqrt(2.);
    unorm[] = (u.x[] - u.y[])/sqrt(2.);
  }
  foreach(){
		fprintf(stderr, "pipj, %d %d\n", point.i, point.j);
		if (point.i + point.j == N)
fprintf (stdout, "%g %g %g %g %g\n", x, y, unorm[], utan[], p[]);
	}
  
  //foreach() {
  //  if (cs[] > 0. && y > x){
  //    theoretical[] = 2*(y - x) ;
  //    e[] = utan[] - theoretical[];
  //  }
  //  else e[] = nodata;
  // }
  // norm n = normf (e);
  // fprintf (stderr, "%d %.3g %.3g %.3g %d %d %d %d %d\n",
  //   N, n.avg, n.rms, n.max, i, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
   
    draw_vof ("cs", "fs");
    squares ("utan", linear = true, spread = -1);
    save ("utan.png");
    dump();
}

/**
![Velocity field](couette45_slip/utan.png)

Slip velocity seems very good. 

~~~gnuplot tangential velocity as a function of channel normal position (x = y).
set xlabel 'x'
set ylabel 'Tangential velocity'
set xrange [-0.5:0.5] 
set yrange [0.3:1.1]
plot 'out' u 1:4 title 'Simulation data', \
     0.2/(0.2 + sqrt(2)/2 * 0.5) title 'Analytical slip velocity'

~~~
*/