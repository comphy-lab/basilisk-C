/**
# Breaking of a waVe on a slope

Une vague se propage en eau peu profonde de profondeur constante et arrive sur une plage en pente droite:
"how many seas must a white dove sail before she sleeps in the sand?"

Pour répondre à la question:
Résolution de Saint Venant simple en 1D avec topographie 
$$
 \left\{\begin{array}{l}
         \partial_t h+\partial_x Q=0\\
	\partial_t Q+ \partial_x \left[\dfrac{Q^2}{h}+g\dfrac{h^2}{2}\right]
= - gh \partial_x Z_b
        \end{array}\right. 
$$
avec au temps t=0, $\eta=0$  et une vitesse nulle, le fond initialement en -1 remonte avec un angle $\alpha$ ensuite. Des vagues rentrent de la gauche.

On utilise le solveur de Saint-Venant et on utilise une grille simple cartésienne: 

## Code */
#include "grid/cartesian1D.h"
#include "saint-venant.h"
/**   declarations ... 
*/
double a,tmax;
/** condition de sortie à gauche, mais surtout condition de caractéristique à gauche
*/
u.n[left] = neumann(0);
u.n[left] = - radiation(a*sin(0.25*t));
/*This can be used to implement open boundary conditions at low Froude numbers.
 The idea is to set the velocity normal to the boundary so that the water level relaxes towards its desired value (ref).
 `radiation(ref) (sqrt (G*max(h[],0.)) - sqrt(G*max((ref) - zb[], 0.)))`
forme du fond, 
*/ 
event init (i = 0)
{
  foreach(){
    zb[] =  -1+  (x>10)*(x-10)/(25);
    u.x[]=0 ;
    h[]=u.x[]+max(-zb[],0);
    }
}
/**
sortie pour gnuplot
*/
event field (t<tmax;t+=.5) {
    printf("p[-15:][-1.1:.25]'-' u 1:3 t'u'w l,''u 1:4 t'z' w l,''u 1:($2+$4) t'eta' w l\n");
    foreach()
    printf (" %g %g %g %g %g \n", x, h[], u.x[], zb[], t);
    printf ("e\n\n");
}
 
event output (t += 5) {
     static FILE * fp1 = fopen ("vals.txt", "w");
  foreach()
    fprintf (fp1, " %g %g %g %g %g \n", x, h[], u.x[], zb[], t);
     fprintf (fp1, "\n\n"); 
    
}
/**
fin des procédures
 
Position du domaine `X0`, taille du domaine `L0`, on est sans dimension `G=1`, amplitude de la vague: `a`. 
*/
int main()
   {  
    X0 = -15.;
    L0 = 60.;
    G = 1;
    N = 1024; 
    a=0.125;
    tmax=300;	
    run();
    }
/**
## run

Pour exécuter, le plus simple:

~~~bash  
qcc -g -O2  -Wall  -o slope_deferle1D slope_deferle1D.c -lm
./slope_deferle1D   | gnuplot  
~~~

ce qui fait apparaître dans la fenêtre la simulation au cours du temps

un plot à différents temps, le temps monte de bas en haut.

~~~gnuplot
set key left
c=.01
p[][-.1:1.5]'vals.txt' u 1:($2+$4+$5*c) i 0 t't=0' w l linec 3,\
         '' u 1:($2+$4+$5*c) i 1 t't=5'w l linec 3, \
         '' u 1:($2+$4+$5*c) i 4 t't=20'w l linec 3,\
         '' u 1:($2+$4+$5*c) i 7 t't=35'w l linec 3,\
         '' u 1:($2+$4+$5*c) i 11 t't=55'w l linec 3,\
         '' u 1:($2+$4+$5*c) i 15 t't=75'w l linec 3,\
         '' u 1:($2+$4+$5*c) i 19 t't=95'w l linec 3
~~~

# Links

* [http://basilisk.fr/sandbox/M1EMN/Exemples/slope.c]()
* [http://basilisk.fr/sandbox/M1EMN/Exemples/slope_inc.c]()

plus fort: 

* [http://basilisk.fr/src/examples/shoal.c]()

*/
