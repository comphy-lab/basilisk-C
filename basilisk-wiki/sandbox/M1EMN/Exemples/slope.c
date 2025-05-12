/**
# Breaking of a wave on a slope

Une vague se propage en eau peu profonde de profondeur constante et arrive sur une plage en pente douce à droite:
"how many seas must a white dove sail before she sleeps in the sand?"

Pour répondre à la question:
Résolution de Saint Venant simple en 1D avec topographie !
$$
 \left\{\begin{array}{l}
         \partial_t h+\partial_x Q=0\\
	\partial_t Q+ \partial_x \left[\dfrac{Q^2}{h}+g\dfrac{h^2}{2}\right]
= - gh \partial_x Z_b
        \end{array}\right. 
$$
avec au temps t=0, $h=1+aexp(-x^2)$  et une vitesse $aexp(-x^2)$ (pour obtenir une petite onde (solution linéarisée) qui va vers le rivage)et un fond $z=0$ pour $x<x_0$ qui remonte avec un angle $\alpha$ ensuite.

On utilise le solveur de Saint-Venant et on utilise une grille simple cartésienne: 

## Code */
#include "grid/cartesian1D.h"
#include "saint-venant.h"
/**   declarations ... 
*/
double a,alpha,tmax;
/** condition de sortie à gauche
*/
u.n[left] = neumann(0);
h[left] = neumann(0);
/**
forme du fond, vitesse = hauteur, donc la vague se propage à droite, plus $a$ est petit, moins elle se casse et moins une vague sort à gauche.
*/ 
event init (i = 0)
{
  foreach(){
    zb[] =   (x>10)*(x-10)*alpha;
    u.x[]=a*exp(-(x+12)*(x+12)) ;
    h[]=fmax(1+u.x[]-zb[],0);
    }
}
/**
sortie pour gnuplot

la première pour un plot à la volée
*/
event field (t<tmax;t+=1) {
    printf("p[-15:][-.1:1.25]'-' u 1:3 t'u'w l,''u 1:4 t'z' w l,''u 1:($2+$4) t'eta' w l\n");
    foreach()
    printf (" %g %g %g %g %g \n", x, h[], u.x[], zb[], t);
    printf ("e\n\n");
}

/**
 
la seconde  pour un plot  plus tard (les deux pourraienet être astucieusement combinés)
*/

event printdatas (t +=10){
    char s[80];
    sprintf (s, "shapes.txt");
    static FILE * fp = fopen (s, "w");
    foreach()
      fprintf (fp, "%g %g %g %g %g\n", x, h[],u.x[], zb[], t);
    fprintf (fp, "\n");
    fflush(fp);
}
event end (t = tmax) {
    printf("#fin");    
}
/**
fin des procédures
*/

/** position du domaine `X0`, taille du domaine `L0`, on est sans dimension `G=1`, amplitude de la vague: `a`. 
*/
int main()
   {  
    X0 = -15.;
    L0 = 60.;
    G = 1;
    N = 1024; 
    a=0.01;
    alpha= 1./(25);
    tmax=150;	
    run();
    }
/**
## run

Pour exécuter, le plus simple:

~~~bash  
qcc -g -O2 -DTRASH=1 -Wall  -o slope slope.c -lm
./slope | gnuplot
~~~

ce qui fait apparaître dans la fenêtre la simulation au cours du temps

si on veut reproduire:

![how many seas must a white dove sail before she sleeps in the sand?](/sandbox/M1EMN/Exemples/Img/slope.gif)

alors il faut faire 

~~~bash  
qcc -g -O2 -DTRASH=1 -Wall  -o slope slope.c -lm
./slope > dmp
echo "set term gif animate; set output 'slope.gif'" > slope.out
cat dmp >> slope.out
rm dmp
cat slope.out | gnuplot
~~~

Un tracé en pose longue:

~~~gnuplot
set xlabel "x"
p[-10:35]'shapes.txt' u 1:($2+$4) t'h' w l
~~~

*/
