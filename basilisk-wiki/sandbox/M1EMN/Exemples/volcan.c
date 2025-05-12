/**
# écoulement de lave visqueuse le long d'un volcan

## problème

On suppose de manière simplifiée et abusive que la lave est un fluide visqueux, 
on va donc mettre un frottement de type Poiseuille en terme de frottement. 
Le volcan est un cône.

On utilise le solveur de Saint-Venant en 2D pour 
résoudre la partie non visqueuse.
 
## Code
*/
#include "saint-venant.h"
/**  quelques déclarations
*/
double tmax,Cf;
#define MAXLEVEL 8
#define MINLEVEL 6
#define LEVEL 6
/**
Conditions aux limites Neumann
*/
u.n[left]  = neumann(0);
u.n[right] = neumann(0);
u.t[left]  = neumann(0);
u.t[right] = neumann(0);
h[left]    = neumann(0);
u.n[top]   = neumann(0);
u.t[top]   = neumann(0);
h[top]     = neumann(0);
h[right]   = neumann(0);
/** 
Initialement un cône pour la topo zb et une hauteur de lave nulle
*/
event init (i = 0)
{ 
   foreach(){
     zb[] = fmax((1.-sqrt(x*x+y*y))+0*(x-2)/10,0);
     h[]= 0;
         }
}
/** 
à chaque pas de temps, on fait sortir de la lave
*/
event source (i++) {
  // source near the summit
  // source pendant un temps fini
 foreach(){
     h[]= ((x-.3)*(x-.3)+y*y < (0.15)*(0.15) ? h[]+0.3*dt*(t<5): h[]);}
} 
/** 
à chaque pas de temps, on met le frottement de Poiseuille
$$
\frac{\partial u}{\partial t} = - 3 \nu  \frac{u}{h^2} 
$$
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
à chaque pas de temps, on adapte le maillage
*/  
event adapt (i++) {
  astats s = adapt_wavelet ({h}, (double[]){1e-3}, MAXLEVEL, MINLEVEL);
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
  output_ppm (etam, fp, mask = m, min = 0, max = 1.005, n = 512, linear = true);}
/** 
et sortie directe gnuplot
*/ 
 event splot(t <= tmax;t+=.125) {
   double dx=L0/pow(2.,LEVEL),dy=dx;
     fprintf (stderr,"--- %lf \n",t); 
     printf("set hidden3d;  unset border ; unset ytics ;unset xtics ; unset ztics;\n");
     printf("set view 60,%lf \n",t/tmax*360);
     printf("sp[:][:][0:2]  '-' u 1:2:($3+$4-.01) not w l linec 1,''u 1:2:4 not w l linec 2\n");
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
taille du domaine et valeur de la friction
*/
int main() {
  L0 = 4.;
  X0 = -1.25;
  Y0 = -2;
  G = 1;
  N = 1 << MAXLEVEL;
  tmax = 20;
  Cf = .1;
  DT = HUGE;
  run();
}

/** 
## Run

Ensuite on compile et on lance:

~~~bash
qcc -fopenmp  -g -O2 -DTRASH=1 -Wall  volcan.c -o volcan
./volcan | gnuplot
~~~

## Plot

Si on veut générer le gif animé

~~~bash
./volcan > v.out
echo "set term gif animate; set output 'v.gif'" > dmp
cat v.out >> dmp
mv dmp v.out 
cat v.out | gnuplot
gifsicle --colors 256 --optimize --delay 1   --loopcount=0 v.gif > volcan.gif 
rm v.out v.gif
~~~

ce qui produit:


![lave visqueuse sur un cône](/sandbox/M1EMN/Exemples/Img/volcan.gif)

 
*/
