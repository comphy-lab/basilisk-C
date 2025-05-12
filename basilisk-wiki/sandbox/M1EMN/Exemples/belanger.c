/** 
# Saint Venant Shallow Water Bélanger relation

## Scope
 
We test the Bélanger relation for velocity and wtaer height before and after a steady standing jump. The file is in french.

Nous vérifions les relation de Bélanger de part et d'autre d'un ressaut. Les relations sont établies pour un ressaut stationnaire de vitesse $W=0$.
 
 ~~~gnuplot hydrolic jump configuration
 set arrow from 4,0 to 4,1 front
 set arrow from 4,1 to 4,0 front
 set arrow from 8,0 to 8,3. front
 set arrow from 0.1,.1 to 9.9,0.1 front
 set arrow from 9.1,.1 to 0.1,0.1 front
 set arrow from 1,0.75 to 3,0.75 front
 set arrow from 1,0.5 to 3,0.5 front
 set arrow from 1,0.25 to 3,0.25 front
 set arrow from 6,0.5 to 7,0.5 front
 set arrow from 6,1.50 to 7,1.5 front
 set arrow from 6,0.25 to 7,0.25 front
 set label  "L0" at 6,.15 front
 set label  "depth h1" at 3.2,1.1 front
 set label  "depth h2" at 8.2,3.2 front
 set label  "water" at 1.,.9 front
 set label  "air" at 1.4,2.05 front
 set xlabel "x"
 set ylabel "h"
 p [0:10][0:6]0 not,(x<5? 1 :3) w filledcurves x1 linec 3 t'free surface'
 reset
 ~~~
 
 
  
## Démonstration 
 
 Dans le repère qui se déplace avec le ressaut à vitesse constante, le ressaut est fixe. On écrit les équations de conservation de part et d'autre du ressaut.
Conservation de la quantité de mouvement
$$U_1^2h_1+g\frac{h_1^2}{2}= U_2^2h_2+g\frac{h_2^2}{2}.
$$ 
 Conservation de la masse:
$$ U_1h_1=U_2h_2, $$
 
à partir de ces relations de conservation nous allons écrire l'expression de $U_2/U_1$ et $h_2/h_1$ en fonction du nombre de Froude $F_1^2=\frac{U_1^2}{gh_1}$.
 Partant de la première multipliée par $h_2$ dans lequel $U_2h_2$ est remplacé par
 $U_1 h_1$ compte tenu de la deuxième
 $$
 U_1^2h_1h_2+g\frac{h_2h_1^2}{2}= U_1^2h_1^2+g\frac{h_2^2}{2}h_2
 $$
 d'où en mettant les vitesses à gauche et la gravité à droite et  comme on reconnaît une identité remarquable avec $(h_2^2-h_1^2)$:
 $$
 U_1^2h_1(h_2- h_1) = g\frac{(h_2-h_1)(h_2+h_1)}{2}
 \;\;\text{ ou }\;\;
 U_1^2  = g\frac{ (h_2+h_1)h_2}{2h_1}.
 $$
 Le nombre de Froude est donc
 $$
 F_1^2=\frac{U_1^2}{gh_1}  = \frac{ (h_2+h_1)h_2}{2h_1^2}
 $$
 si $h_2>h_1$ le nombre $F_1$ est supérieur à 1.
 On pourrait faire de même pour l'indice 1, on trouve alors
 $$
 F_2^2=\frac{U_2^2}{gh_2}  = \frac{ (h_2+h_1)h_1}{2h_2^2},
 $$
 cette fois,
 si $h_2>h_1$ le nombre $F_2$ est inférieur à 1. On constate que
 $$
 \frac{F_2}{F_1}=\big(\frac{h_1}{h_2}\big)^{3/2}.
 $$
 Reprenons l'équation de $F_1^2$ et exprimons la comme un polynôme en $h_2/h_1$ soit:
 $$\big(\frac{h_2}{h_1}\big)^{2} + \big(\frac{h_2}{h_1}\big) - 2 F_1^2=0$$
 cette équation du second degré a une racine positive qui est:
 $$\big(\frac{h_2}{h_1}\big) = \frac{-1+\sqrt{1+8 F_1^2}}{2}.$$
 On a donc maintenant toutes les quantités avant et après le ressaut fixe:
 $$\big(\frac{h_2}{h_1}\big)= \big(\frac{U_1}{U_2}\big)=
 \big(\frac{F_1}{F_2}\big)^{2/3}
 = \frac{-1+\sqrt{1+8 F_1^2}}{2}.$$
 
 
Remarque: ces équations classiques sont établies pour un ressaut fixe.
On peut faire bouger le ressaut en ajoutant une vitesse d'ensemble $W$ qui sera la vitesse d'avancée du ressaut. 
En effet les équations de Saint Venant sans frottement sur fond plat sont invarianets par changement de référentiel Galiléen.


# Code 
 */
#include "grid/cartesian1D.h"
#include "saint-venant.h"

double tmax,Fr,deltah,x_rs,W;

int main(){
  X0 = 0;
  L0 = 10;              
  N = 512;
  G = 1;
  tmax = 5;
  Fr = 2.5;      // Nombre de Froude initial à gauche avec indice 1
  W = 0;         //vitesse du ressaut, ici 0, on peut faire varier W
  x_rs= L0/2;    // position nitiale du saut dans [X0;X0+L0]
  DT = HUGE;
  FILE * fp = fopen ("WT.txt", "w");
  fclose(fp);    // fichier pour la position du ressaut 
  run();
}
/** 
sortie et entrée libres
*/
h[left] = h[];
h[right] = h[];
u.n[left] = neumann(0);
u.n[right] = neumann(0) ;

/**
 initialisation: Bélanger avec $h_1=1$ et  $h_2$ est calculé 
 (attention compiler avec `-disable-dimensions` !)
 $$\big(\frac{h_2}{h_1}\big) =  
  \frac{-1+\sqrt{1+8 F_1^2}}{2}$$
 `Fr` est $F_1$
*/
event init (i = 0)
{
  foreach(){
    zb[]=0.;
    // formule de Bélanger avec h1=1;
    deltah=(((-1 + sqrt(1+8*Fr*Fr))/2)-1);
     h[]=1+(((-1 + sqrt(1+8*Fr*Fr))/2)-1)*(1+tanh((x-x_rs)/.01))/2;
    // vitesse associée par conservation du flux Fr/h + translation 
     u.x[]=Fr/h[]+W;
    }
 
#ifdef gnuX
  printf("\nset grid\n");
#endif
}

/** position du choc $x_c$ */
event findS(t<tmax;t+=0.2) {
   FILE * fp = fopen ("WT.txt", "a");
   double xc=0.;
// si la hauteur dépasse un peu h1=1, c'est le ressaut   
   foreach(){
     if(h[]>1.0001){  xc=xc;} else { xc=x;};
   }
    fprintf (fp, "%g %g   \n", t, xc );
    fflush(fp);
}
#ifdef gnuX
// plot à la volée avec gnuplot
event plot (t<tmax;t+=0.01) {
    printf("set title 'Ressaut en 1D --- t= %.2lf '\n"
   "p[%g:%g][-.5:5]  '-' u 1:($2+$4) t'free surface' w l lt 3,"
   "'' u 1:3 t'velocity' w l lt 4,\\\n"
   "'' u 1:4 t'topo' w l lt -1,\\\n"
     "'' u 1:(sqrt($3*$3/($2*%g))) t'Froude' w l lt 2,\\\n"
     "'' u 1:($2*%g+$3*$3/2) t'Charge' w l lt 5,\\\n"
     "'+' u (%lf):(1+%lf) t'jump theo' \n",
           t,X0,X0+L0,G,G,x_rs+W*t,deltah/2);
    foreach()
    printf ("%g %g %g %g %g\n", x, h[], u.x[], zb[], t);
    printf ("e\n\n");
}
#else
// affichage uniquemeny du dernier temps
event end(t=tmax ) {
    foreach()
    printf ("%g %g %g %g %g\n", x, h[], u.x[], zb[], t);
}
#endif
/**
# Run and Results
## Run
To compile and run:

~~~bash
 qcc -disable-dimensions -DgnuX=1  belanger.c -o belanger -lm; ./belanger | gnuplot
 
 make belanger.tst
 make belanger/plots;
 source  c2html.sh belanger
 
~~~

## Plot of jump
 
Plot of jump comparision 
with analytical result 
$$\big(\frac{h_2}{h_1}\big)   
 = \frac{-1+\sqrt{1+8 F_1^2}}{2}$$
 for `Fr=2.5` and `h1=1` we have 
 `h2=(((-1 + sqrt(1+8*Fr*Fr))/2))`
~~~gnuplot  result, free surface (blue) and bottom (black)
set xlabel "x"
 Fr=2.5
 h2=(((-1 + sqrt(1+8*Fr*Fr))/2))
 p [:][-1:5] 'out' t'free surface' lc 3,'' u 1:3 w p t'speed', '' u 1:4 t'zb' w l lc -1,(x>5?h2:NaN) t'h2'  w l lc 1,(x<5?1:NaN)  t'h1' w l lc 1
 
~~~


Plot de la position du choc en fonction du temps. 
   Dans le cas de Bélanger $W=0$.  
   En jouant avec $W$ dans le code, on voit que les ressauts mobiles ($W\ne 0$) sont par translation en $W$ des ressauts fixes...

~~~gnuplot 
set xlabel "t"
set ylabel "jump position"
p[0:5][0:10] 'WT.txt' u 1:2 w lp
~~~


# Liens pour aller plus loin
 
 Cas dispersif, tension de surface explicite

 [http://basilisk.fr/sandbox/M1EMN/Exemples/belangerdisp.c]()
 
 cas avec frottement turbulent 


 [http://basilisk.fr/sandbox/M1EMN/Exemples/belangerturb.c]()

 Cas dispersif, mascaret à partir des équations de Boussinesq (implicite) avec frottement turbulent 

 [http://basilisk.fr/sandbox/M1EMN/Exemples/ressaut_mascaret.c]()


 


en multicouche laminaire....

 [http://basilisk.fr/sandbox/M1EMN/Exemples/higuera.c]()
 
# Bibliographie
 
* [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
"Equations de Saint Venant et application, Ecoulements en milieux naturels" Cours MSF12, M1 Sorbonne-Université
*/


