/**
# Une vague se propage en eau peu profonde de profondeur constante et arrive sur une plage en pente droite:
la plage n'est pas parallele à la vague initiale.
Résolution de Saint Venant simple en 2D. 

avec au temps t=0, $h=1+aexp(-x^2)$ et un fond $z=0$ pour $x<x_0$ qui remonte avec un angle $\alpha$ ensuite.


On utilise le solveur de Saint-Venant et on utilise l'adaptativite: 
*/
#include "saint-venant.h"

double tmax,Cf;
#define LEVEL 9

/** CL
*/
u.n[left] = neumann(0);
u.t[left] = neumann(0);
h[left] = h[];
u.n[top] = neumann(0);
u.t[top] = neumann(0);
h[top] = h[];

/** fond plat puis en pente avec un angle de côté, une vague arrive 
*/
event init (i = 0)
{ double xi;
  foreach(){
     xi=x+y/4.+1;   // fond incline
     zb[] =   (xi>-1)*(xi+1)/7.;    
     xi=x+7;
     h[]=fmax(1+.2*exp(-2*xi*xi)-zb[],0);
     u.x[]=.2*exp(-2*xi*xi);     
         }
     boundary (all); 
}
/** en fait 0
*/
event friction (i++) {
  // Poiseuille friction (implicit scheme)
  foreach()
  {double ff = h[] < dry ? HUGE : (1. + Cf*dt/h[]/h[]);
   foreach_dimension()
   u.x[] /= ff;
    }
}
  
/**
adapatif
*/     
event adapt (i++) {
  astats s = adapt_wavelet ({h}, (double[]){1e-3}, LEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
/** film et plot
*/
event movies (t <= tmax;t+=.1) {
  static FILE * fp = popen ("ppm2mpeg > eta.mpg", "w");
  scalar m[], etam[];
  
  foreach() {
    etam[] = eta[]*(h[] > dry);
    m[] = etam[] - zb[];
  }
  boundary ({m, etam});
  output_ppm (etam, fp, mask = m, min = 0.9, max = 1.2, n = 512, linear = true);}


 
event cut(t <= tmax;t+=.5) {
	double dx=8*L0/pow(2.,LEVEL),dy=dx;
	fprintf (stderr,"--- %lf \n",t); 
	  printf("set hidden3d; sp[][][0:1.25]  '-' u 1:2:($3+$4) not w l linec 3 ,''u 1:2:4 not w l linec 2\n");
    for(double x=-L0/2.;  x<L0/2.;x+=dx)
     {for(double y=-L0/2.;y<L0/2.;y+=dy)
      { 
     printf (" %g %g %g  %g \n", x, y, interpolate (h, x, y),  interpolate (zb, x, y) );}
     printf ("\n");
     } 
     printf ("e\n");
}
/**
paramètres
*/
 
int main() {
  L0 = 20.;
  X0 = -L0/2.;
  Y0 = -L0/2.;
  G = 1;
  N = 1 << LEVEL;  // ie 2^LEVEL
  tmax = 40;
  // friction nulle
  Cf = 0.000000;
  
  run();
}



/**
Pour exécuter, le plus simple:

~~~bash  
qcc -g -O2 -DTRASH=1 -Wall  -o slope_inc slope_inc.c -lm
./slope_inc | gnuplot
~~~

ce qui fait apparaître dans la fenêtre la simulation au cours du temps

si on veut reproduire:

![plage](/sandbox/M1EMN/Exemples/Img/slope_inc.gif)

 alors il faut faire 

~~~bash  
qcc -g -O2 -DTRASH=1 -Wall  -o slope_inc slope_inc.c -lm
./slope_inc > dmp
echo "set term gif animate; set output 'slope_inc.gif'" > slope_inc.out
cat dmp >> slope_inc.out
rm dmp
cat slope_inc.out | gnuplot
~~~

*/
