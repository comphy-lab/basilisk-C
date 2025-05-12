/**
# A more real case for imposing flow rates

In this example, we use the discharge() function in order to impose
different flow rates on different rivers situated on the same boundary of
a Saint-Venant simulation. Wa also use adaptative quadtrees.*/

#include "saint-venant.h"
#include "discharge.h"

/**
The domain is 400 metres squared, centered on the origin. Time is in
seconds. The minimum level of refinement is $2^4$ and the maximum one
is $2^8$. We fix the maximum error on the water elevation to a low
value ($1e^{-3}$). Therefore, the grid is refined to its maximum value
where it has presence of fluid. */

#define MINLEVEL 4
#define MAXLEVEL 8
#define HMAXE 1e-3

int main(){
  L0 = 20.;
  X0 = - L0/2.;
  Y0 = - L0/2.;
  G = 9.81;
  N = 1 << MAXLEVEL; 
  run();
}

/**
## Initial conditions

We chose a complicated river bed with height "valleys" flowing in four
lakes (see figure below).  We set the field *river[]* with the values
1 and 2 in order to distinguish the rivers on the same boundary. */

event init (i = 0){
  /**
  We start with a dry riverbed, so that the problem does not have a
  natural timescale the Saint-Venant solver can use. We set a maximum
  timestep to set this timescale. */
  DT = 1e-2;
  foreach() {
    
    /** We define a topography with two valleys at each side of the
	domain and one lake at each quarter of the map. See the map
	below.*/
    
    zb[] = 5e-4*(x*y*x*y - x*x - y*y + 2000*(cos(x/10.*2*3.14)+cos(y/10.*2*3.14)));

    /**
    ~~~gnuplot Topography 
    reset
    set isosample 20,20
    set pm3d depthorder
    f(x,y) = 5e-4*(x*y*x*y - x*x - y*y + 2000*(cos(x/10.*2*3.14)+cos(y/10.*2*3.14))) 
    splot [-10:10][-10:10] f(x,y) with pm3d t'topography' 
    ~~~ 
    */

    /**The scalar river[] is defined by default in the
    [discharge.h](http://basilisk.fr/src/discharge.h) library. We use
    it to store the location of each rivers. Its value is 1 for the
    river at north-east and south-west of the map and 2 for the rivers
    at north-west and south-east.*/
    
    river[] = x*y > 0 ? 1 : 2;
    
    /**
    ~~~gnuplot river[] 
    reset
    set isosample 20,20
    set pm3d map
    f(x,y) = (x*y < 0)+1
    splot [-10:10][-10:10] f(x,y) t'river[]' 
    ~~~
    */
    
  }
  boundary ({zb, river});
}

/**
## Adaptation

We use the adap_wavelet function to adapt the grid only on the wet
cells. You can find a full explaination of this process in the
[tsunami
example](http://basilisk.fr/src/examples/tsunami.c#adaptation)*/

void adapt_H() {
#if QUADTREE
  scalar h1[];
  foreach(){
    h1[] = h[];
  }
  boundary ({h1});
  adapt_wavelet ({h1}, (double[]){HMAXE} ,MAXLEVEL,MINLEVEL);
#endif
}

/**
We refine at each time step.*/

event do_adapt_H (i++) adapt_H();

/**
## Boundary conditions

We impose inflow on all boundaries. In addition, on each boundary, the
tangential velocity is set to zero and a neumann condition is fixed on
the normal velocity. */

u.n[top] = neumann(0);
u.t[top] = dirichlet(0);

u.n[bottom] = neumann(0);
u.t[bottom] = dirichlet(0);

u.n[right] = neumann(0);
u.t[right] = dirichlet(0);

u.n[left] = neumann(0);
u.t[left] = dirichlet(0);


/**
We impose different flow rates : $1 m^3.s^{-1}$ for the rivers in the
north-east quarter, $1.5 m^3.s^{-1}$ for the ones at north-west, $2
m^3.s^{-1}$ for the ones at south-west, and 0 at south-east. The
process is more detailed in the [multi-river test
case](http://basilisk.fr/src/test/multiriverinflow.c).*/

double eta1, eta2, eta3, eta4, eta5, eta6;

event inflow (i++) {
  // Top boundary
  eta1 = discharge (0.5, top, 1); // N-E
  eta2 = discharge (1, top, 2); // N-W  
  h[top] = max ( (river[] == 1 ? eta1 : eta2) - zb[] , 0) ;
  eta[top] = max ((river[] == 1 ? eta1 : eta2), zb[]);

  // Left boundary
  eta4 = discharge (1.5, left, 1); // S-W
  eta3 = discharge (1, left, 2); // N-W
  h[left] = max ((river[] == 1 ? eta4 : eta3) - zb[] , 0);
  eta[left] = max ((river[] == 1 ? eta4 : eta3), zb[]);

  // Bottom boundary
  eta5 = discharge (1.5, bottom, 1); //S-W
  h[bottom] = max ((river[] == 1 ? eta5 : zb[]) - zb[] , 0);
  eta[bottom] = max ((river[] == 1 ? eta5 : zb[]), zb[]);

  // Right boundary
  eta6 = discharge (0.5, right, 1); //N-E
  h[right] = max ((river[] == 1 ? eta6 : zb[]) - zb[] , 0);
  eta[right] = max ((river[] == 1 ? eta6 : zb[]), zb[]);
}

/**
## Outputs

We compute the evolution of the water volumes in each quarter of the map. */

event volume (i += 10) {
  double volume1 = 0, volume2 = 0, volume3 = 0, volume4 = 0;
  foreach() {
    double dv = h[]*sq(Delta); 
    if (y > 0 && x > 0)
	volume1 += dv; // N-E
    else if(y > 0 && x < 0 )
      volume2 += dv; // N-W
    else if( y < 0 && x < 0 )
      volume3 += dv; // S-W
    else  volume4 += dv; // S-E-
  }
  fprintf (stderr, "%g %g %g %g %g\n",
	   t, volume1, volume2, volume3, volume4);
}
  
/**
We use gnuplot to produce an animation of the water surface. */

event animation (t <= 1; t += 0.01) {
  double dx = 2.*L0/(1 << MAXLEVEL), dy = dx;
  printf ("set title 't = %.3f'\n"
	  "sp [%g:%g][%g:%g][-5:5] '-'"
	  " u 1:2:($3+$4-.05) t 'free surface' w l lt 3,"
	  " '' u 1:2:4 t 'topography' w l lt 2\n",
	  t, X0, -X0, Y0, -Y0);
  for (double x = X0;  x <= X0 + L0; x += dx) {
    for (double y = Y0; y <= Y0 + L0; y += dy)
      printf ("%g %g %g %g\n",
	      x, y, interpolate (h, x, y),  interpolate (zb, x, y));
    putchar ('\n');
  }
  printf ("e\n"
	  "pause %.5lf \n\n", 0.);
}

/**
We record a movie for the field *level* which is storing the level of
refinement.*/

event movies (t+=0.01) {
  static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, fp1, min = MINLEVEL, max = MAXLEVEL, n = 512);
}

/**
## Results

~~~gnuplot Evolution of water volumes in both rivers
reset
set key top left
set xlabel 'Time'
set ylabel 'Volume'
f(x)=a*x + b
fit [0.1:0.8] f(x) './log' u 1:2 via a,b
g(x)=c*x + d
fit [0.1:0.8] g(x) './log' u 1:3 via c,d
h(x)=k*x + l
fit [0.1:0.8] h(x) './log' u 1:4 via k,l
i(x)=m*x
fit [0.1:0.8] i(x) './log' u 1:5 via m
title1=sprintf("f(x) = %1.3f*t + %1.3f", a,b)
title2=sprintf("g(x) = %1.3f*t + %1.3f", c,d)
title3=sprintf("h(x) = %1.3f*t + %1.3f", k,l)
title4=sprintf("i(x) = %1.3f*t", m)
plot [0:0.8] './log' u 1:2 t 'North-east (q=1)', \
     f(x) t title1, \
     './log' u 1:3 t 'North-west (q=2)', \
     g(x) t title2, \
      './log' u 1:4 t 'South-west (q=3)', \
     h(x) t title3, \
     './log' u 1:5 t 'South-east (q=0)', \
     i(x) t title4
~~~

~~~gnuplot Animation of the free surface.
reset
set view 40,15
set xlabel 'x'
set ylabel 'y'
set hidden3d
set term gif animate
set output 'movie.gif'
load './out'
~~~

![[Animation](morerivers/level.mpg) of the level of refinement. Dark blue is 4
and dark red is 8.](morerivers/level.png)
*/
