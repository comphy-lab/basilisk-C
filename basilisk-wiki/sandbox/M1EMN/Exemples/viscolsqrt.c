/** 
#Viscous Collapse (Viscous Slump): 

## Problème de Effondrement visqueux sur une pente (sploutch)
Effondrement visqueux: il s'agit de l'effondrement d'un tas 1D initialement rectangulaire et disposé sur une pente $Z(x)$.
Comme l'écoulement est très visqueux, l'inertie est vite négligeable:
on a vite un équilibre entre le terme de pente (de gravité) qui fait tomber et le frottement 
visqueux qui freine. 
Cet écoulement pourrait modéliser de la lave s'écoulant le long d'un volcan (il est générique en géophysique, le terme à changer est celui du frottement, nous verrons sur les autres pages le cas de bingham, fluide à seuil).


Les équations de Saint Venant
$$\frac{\partial h} {\partial t} + \frac{\partial (Q)}{\partial x}=0$$
$$\frac{\partial Q } {\partial t} + \frac{\partial }{\partial x}(\Gamma\frac{Q^2}{ h }  +
\frac{g}{2} (h_{\;}^2))=  -gZ' h  - Cf \frac{Q}{h_{\;}^2}$$
 avec $\Gamma=1$ (au lieu de (5/4)) et $C_f=3\nu$,
sont résolues avec une masse initiale 
$$A=\int_{0}^{x_1} h(x,0) dx$$
sur une pente $\alpha$. Elles
sont résolues pour la partie non visqueuse par `saint-venant.h` 


La partie visqueues est résolue 
par split  
 $$\frac{\partial u}{\partial t} = -C_f \frac{u}{h^2} \text {     discretisation semi implicite}
 \frac{  u^{n+1} - u^n }{\Delta t} = -C_f \frac{u^{u+1}}{h^{n2}}$$
 onimpose $u^{n+1}=0$ quand $h^n$ très petit  (paramètre `dry`).
 

Pour les pentes fortes telles que $|Z'|\gg \partial_x h$, on peut simplier les équations car on n'a plus que l'équilibre entre la pente et le frottement: alors  $Q= -gZ' h^3/(3 \nu)$ et l'équation de la masse devient 
$$\frac{\partial h} {\partial t} - \frac{gZ'}{3\nu}  \frac{\partial }{\partial x} h_{\;}^3   =0,$$
l'équation est de la forme ($k=gZ'/\nu$) d'une équation de transport 
$$ \frac{\partial h} {\partial t} +k  h_{\;}^2 \frac{\partial h}{\partial x} =0$$

Huppert en trouve une solution par les caractéristiques, ici, nous nous proposons de trouver la soution sous forme auto semblable.

Par changement d'échelle, l'invariance par dilatation donne $H/T=H^3/X$, la 
masse conservée  donne $HX=1$, ce qui permet de trouver 
$H=T_{\;}^{-1/3}$ et $X=T_{\;}^{1/3}$
donc une solution auto semblable de la forme $t_{\;}^{-1/3}F(x/t_{\;}^{1/3})$ convient,
 par subsitution
   $$ -\eta  F'(\eta )-F(\eta )+3 k F(\eta )^2 F'(\eta )=0, \text{ soit } (\eta F(\eta)-k F(\eta )^3)'=0
$$
avec $F(0)=0$ on a la forme implicite
$\eta=k (F(\eta)^2)$
et donc la solution est 
$F(\eta)=\sqrt{(\eta)/k},$
la solution du problème 
est simplement pour $0<x<x_f$
$$h =t_{\;}^{-1/3} \sqrt{(x)t_{\;}^{-1/3})/k}= \sqrt{\frac{x}{kt}  }$$

 


On s'est donné une masse initiale $A=\int_{0}^{x_1} h(x,0) dx$ 
qui se déplace entre $0$ le front arrière qui ne bouge pas et $x_f $ position du front,  la masse est donc ultérieurement
$$\int_0^{x_f} \sqrt{\frac{x}{kt}} dx = (2/3) \sqrt{x_f^3/(kt)}$$ qui est $A$ par conservation de la masse
soit $x_f = (\frac{9 A_{\;}^2kt}4)^{1/3}$

## Code
*/
#include "grid/cartesian1D.h"
#include "saint-venant.h"

double Cf,h0,alpha,tmax;
#define MAXLEVEL 13
#define MINLEVEL 11
/**
On se donne une masse initiale $A=\int_{0}^{x_1} h(x,0) dx$ sur une pente $\alpha$.
Initialisation d'un tas de masse 1
*/
event init (i = 0)
{
  foreach() {
    zb[] = -(x-X0)*alpha;
     h[] = h0*(x>0)*(x<=1./h0);
  }
}

u.n[right] = neumann(0);
h[right] = neumann(0);

event friction (i++) {
  // Poiseuille friction (implicit scheme)
  foreach()
  {  	 double ff = h[] < dry ? HUGE : (1. + Cf*dt/h[]/h[]);
    u.x[] /= ff;}
  boundary ({u.x});
}

 event plot (t=0.000001;t < tmax;t+=50) {	
 	 FILE * fp  = fopen ("sol.dat", "w");
//a=.584;p[][]'sol.dat'u ($1/$3**(1./3)):($2*$3**(1./3)),'sol13.dat' u ($1/$3**(1./3)):($2*$3**(1./3)),a*sqrt(x)*(x<(3./2./a)**(2./3))
   fprintf (stderr," %g \n", t);
   printf ("a=.584;p[:][0:.25] '-' u 1:2 t'calc' w l,''u 1:(a*( sqrt($1/$3))*($1<((9/4/a/a*$3)**(1./3))))t'selfs' w l \n");
     for(double x=X0;x<L0+X0;x+=L0/N){
       printf (" %g %g %g \n",x , interpolate (h, x,t),t);
       fprintf (fp,"%g %g %g \n",x , interpolate (h, x,t),t);}
       fclose(fp);
       // foreach()
       // printf (" %g %g %g \n",x ,h[],t);
       printf ("e\n");
}
//event adapt (i++) { adapt_wavelet ({h}, (double[]){1e-3}, MAXLEVEL, MINLEVEL); }
//event adapt (i++) {
//  astats s = adapt_wavelet ({h}, (double[]){1e-5}, MAXLEVEL, MINLEVEL);
 // fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
//}

int main() 
{
  L0 = 30;
  X0 = -2.;
  G = 1;
  N = 512;
  DT=0.001;
  h0 = .25;
  Cf = 3.; // valeur de 3 nu
  alpha=1;
  tmax = 3000/2;

  run();
}


/** 
## Run
Ensuite on compile et on lance:

~~~bash
qcc -fopenmp  -g -O2 -DTRASH=1 -Wall  viscolsqrt.c   -o viscolsqrt
./viscolsqrt | gnuplot
~~~

avec make
 
~~~bash
 make viscolsqrt.tst
 make viscolsqrt/plots
 make viscolsqrt.c.html
 
 source c2html.sh viscolsqrt
 
~~~
 
## Plot

Si on veut générer la figure

~~~bash
./viscolsqrt > v.out
echo "set term png;set output 'slope.png';set multiplot; set nokey ; set xrange[-5:25] ; set yrange [:1];" > dmp
cat v.out >> dmp
mv dmp v.out 
cat v.out | gnuplot
~~~


ce qui produit la comparaison entre la solution semblable et le calcul: 
![collapse](/sandbox/M1EMN/Exemples/Img/slope.png)


Les fronts en fonction du temps sortent de `out`

~~~gnuplot heap as a function of time
set xlabel "x"
set ylabel "h(x,t)"
p[][0:.3]'out' u 1:2 t'h(x,t)' w l
~~~

## Comparaison avec la solution analytique autosemblable $t_{\;}^{-1/3} F(\eta)$

Pour les pentes fortes la solution auto semblable est de la forme $t_{\;}^{-1/3}F(x/t_{\;}^{1/3})$.
La solution du problème 
est simplement pour $0<x<x_f$
$$h =t_{\;}^{-1/3} \sqrt{(x)t_{\;}^{-1/3})/k}= \sqrt{\frac{x}{kt}  }$$
avec $A=\int_{0}^{x_1} h(x,0) dx$ 
 $x_f = (\frac{9 A_{\;}^2kt}4)^{1/3}$

 
  
~~~gnuplot selfsimilar solution
a=.584;
a=1;
set key left
p[-0.25:2.25][0:1.5]'sol.dat'u ($1/$3**(1./3)):($2*$3**(1./3)) w l ,a*sqrt(x)*(x<(3./2./a)**(2./3)) t'anal'
~~~

plus on augmente la résolution, moins la partie gauche va bouger: ici on voit que le bloc a glissé, ce qui est une erreur numérique.
Ce glissement est dû à l'étape visqueuse 

# Links

* same example with [Multilayer](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapsesqrt_ML.c)
* same example with [Bingham](http://basilisk.fr/sandbox/M1EMN/Exemples/bingham_collapse_noSV.c)

# Bibliographie
* [Huppert H.](http://www.itg.cam.ac.uk/people/heh/Paper49.pdf)
 "Flow and instability of a viscous current along a slope"
 Nature volume 30 1982 p 427  
* [M2EMN non Newtonian flows](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/mainM2EMN.pdf) Tas visqueux sur fond incliné
 

*/
