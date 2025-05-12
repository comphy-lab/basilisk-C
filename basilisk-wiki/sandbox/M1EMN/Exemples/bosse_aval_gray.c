/**
#Effondrement d'un tas de sable dans une gouttière

Résolution des équations de Saint Venant 2D avec un frottement constant de Coulomb, un tas initial disposé sur une pente incurvée en forme de gouttière est relâché au temps t=0, il s'effondre et finit par s'arrêter sur le sol plat.


Grâce à:
*/
#include "saint-venant.h"
/** 
on résout le problème de fluide parfait (sans la friction µ qui sera traîtée plus loin)$$
 \left\{\begin{array}{l}
         \frac{\partial }{\partial t}  h \; +\; \frac{\partial }{\partial x} Q_x
\;+\; \frac{\partial }{\partial y} Q_y=0\\
	 \frac{\partial }{\partial t} Q_x+  
 \frac{\partial }{\partial x}  \dfrac{Q_x^2}{h}+
 \frac{\partial }{\partial y}  \dfrac{Q_xQ_y}{h}+
 \frac{\partial }{\partial x}g\dfrac{h^2}{2}
= - gh \frac{\partial }{\partial x} Z-\mu gh\frac{Q_x}{|Q|}\\
	\frac{\partial }{\partial t} Q_y+  
 \frac{\partial }{\partial x}  \dfrac{Q_xQ_y}{h}+
 \frac{\partial }{\partial y}  \dfrac{Q_y^2}{h}+
 \frac{\partial }{\partial y}g\dfrac{h^2}{2}
= - gh \frac{\partial }{\partial y} Z-\mu g h \frac{Q_y}{|Q|}\\
        \end{array}\right., 
$$
*/
double tmax;
#define MAXLEVEL 7
#define MINLEVEL 4
#define LEVEL 6
/**
Conditions aux limites dirichlet en vitesse
*/
u.n[left] =0;
u.t[left] =0;
u.n[top] =0;
u.t[top] =0;
u.n[right] =0;
u.t[right] =0;
/** 
Une gouttière pour la topo zb et un tas de grains h relâchés
*/
event init (i = 0)
{ double wg2=1,Ltas=0.5,pmax=.25,xtas=1; 
   foreach(){
     zb[] = fmin(2,fmax(0,(2+(y*y/2/1.2)*(y*y<wg2)+(y*y>wg2)*wg2/2/1.2 -x*tan(40*pi/180))));
     h[]=(sqrt((x-xtas)*(x-xtas)+y*y)<Ltas)? pmax*(1 - sqrt((x-xtas)*(x-xtas)+y*y)/Ltas): 0 ; 
         }
}
/** 
à chaque pas de temps, on résout le problème de Riemann, puis on met le frottement de  Coulomb
$$
\frac{\partial u}{\partial t} = - \mu g \frac{u}{|u|} 
$$
*/
event friction (i++) {
	scalar q[];
    foreach() q[]=h[]*sqrt(u.x[]*u.x[]+u.y[]*u.y[]);
  
  foreach()
  { double mu;//,I;
// on peut tester des dépendances de mu en I  
//  int nbrg=32;
//    	I=(1./nbrg)*5./2.*q/pow(h+.0000000001,2.5) 
//     	mu = (0.4 + 0.26/(0.4/I +1)) 
//      mu = (0.4 + 0.26*exp(-0.136/I))  
        mu = .45;
   double ff = h[] < dry ? HUGE : (1. + dt*mu*G/(sqrt(u.x[]*u.x[]+u.y[]*u.y[])+.0000000001));
   foreach_dimension()
     u.x[] /= ff;
    }
}
/** 
à chaque pas de temps, on adapte
*/  
event adapt (i++) {
  astats s = adapt_wavelet ({h,zb}, (double[]){1e-3}, MAXLEVEL, MINLEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
/** 
Film
*/
event movies (t <= tmax;t+=.01) {
  static FILE * fp = popen ("ppm2mpeg > eta.mpg", "w");
  scalar m[], etam[];
  
  foreach() {
    etam[] = eta[]*(h[] > dry);
    m[] = etam[] - zb[];
  }
  boundary ({m, etam});
  output_ppm (h, fp, mask = m, min = 0, max = 0.1, n = 1024, linear = true);}
/** 
et sortie directe gnuplot
*/ 
 event splot(t <= tmax;t+=.05) {
   double dx=L0/pow(2.,LEVEL),dy=dx;
     fprintf (stderr,"--- %lf \n",t); 
     printf("set hidden3d;  unset border ; unset ytics ;unset xtics ; unset ztics;\n");
   //  printf("set view 30,%lf \n",t/tmax*360);
     printf("sp[:][:][0:3]  '-' u 1:2:($3+$4-.01) not w l linec -1,''u 1:2:4 not w l linec 2\n");
    for(double x=X0;  x<X0+L0;x+=dx)
     {for(double y=Y0;y<Y0+L0;y+=dy)
      { 
     printf (" %g %g %g  %g \n",x,y,interpolate (h, x, y),interpolate (zb, x, y) );}
     printf ("\n");
     } 
     printf ("e\n");
}
/** 
principal
taille du domaine et paramètres
*/ 
int main() {
  L0 = 5.;
  X0 = 0;
  Y0 = -2.5;
  G = 1;
  N = 1 << MAXLEVEL;
  tmax = 8;
  DT = HUGE;
  run();
}

/** 
Ensuite on compile et on lance:

~~~bash
qcc -fopenmp  -g -O2 -DTRASH=1 -Wall  bosse_aval_gray.c -o bosse_aval_gray -lm
./bosse_aval_gray | gnuplot
~~~

si on veut garder un film

~~~bash
./bosse_aval_gray > v.out
echo "set term gif animate; set output 'bosse.gif'" > dmp
cat v.out >> dmp
mv dmp v.out 
cat v.out | gnuplot
gifsicle --colors 256 --optimize --delay 1   --loopcount=0 bosse.gif > bosse_aval_gray.gif  
rm bosse.gif
rm v.out
~~~


ce qui produit:

![gouttière et granulaire](/sandbox/M1EMN/Exemples/Img/bosse_aval_gray.gif)

# Bibliographie
 * M. WIELAND, J. M. N. T. GRAY AND K. HUTTER
Channelized free-surface flow of cohesionless granular avalanches in a chute
with shallow lateral curvature
J. Fluid Mech. (1999), vol. 392, pp. 73–100.  
*/
