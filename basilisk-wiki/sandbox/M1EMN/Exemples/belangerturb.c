/** 
# Saint- Venant (Shallow Water) Bélanger relation with friction and slope

## Scope
 
We test the Bélanger relation for velocity and water height before and after a steady standing jump.
 We study the influence of friction (turbulent $u^2$) and driving angle of topography ($-Z'_b>0$) (gravity is tilted)
 
 
 ~~~gnuplot hydrolic jump configuration
 set arrow from 2,4 to 3,1.5 front
 set label  "g" at 2.5,2 front
 set label  "gravity tilted" at 2.5,4 front
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
 set label  "depth h1" at 3.6,1.5 front
 set label  "depth h2" at 8.2,3.2 front
 set label  "water" at 1.,.9 front
 set label  "air" at 1.4,2.05 front
 set xlabel "x"
 set ylabel "h"
 p [0:10][0:6]0 not,(x<5? 1+x/25 :3) w filledcurves x1 linec 3 t'free surface'
 reset
 ~~~
 
## Equations

 
 Saint-Venant

 $$\frac{\partial }{\partial t}  h \; +\; \frac{\partial }{\partial x} uh=0$$
 $$ \frac{\partial }{\partial t} hu +
 \frac{\partial }{\partial x}  \dfrac{(hu)^2}{h} +
 \frac{\partial }{\partial x}g\dfrac{h^2}{2}
 = - gh \frac{d }{d x} Z_b-C_f {|u|}u
 + \frac{\partial }{\partial x}(\nu_e h \frac{\partial }{\partial x}u)
 $$
 over an inclined plane $Z_b=-Z'_b x$ with the $-Z'_b$ positive angle of the tilted plane.
 
 
  
## Démonstration 
### Bélanger
Conservation de la quantité de mouvement
$$U_1^2h_1+g\frac{h_1^2}{2}= U_2^2h_2+g\frac{h_2^2}{2}
$$ 
 Conservation de la masse:
$$ U_1h_1=U_2h_2 $$
 
à partir de ces relations de conservation nous allons écrire l'expression de $U_2/U_1$ et $h_2/h_1$ en fonction du nombre de Froude $F_1^2=\frac{U_1^2}{gh_1}$.
 
 cette équation du second degré a une racine positive qui est:
 $$\big(\frac{h_2}{h_1}\big) = \frac{-1+\sqrt{1+8 F_1^2}}{2}$$
 On a donc maintenant toutes les quantités avant et après le ressaut fixe:
 $$\big(\frac{h_2}{h_1}\big)= \big(\frac{U_1}{U_2}\big)=
 \big(\frac{F_1}{F_2}\big)^{2/3}
 = \frac{-1+\sqrt{1+8 F_1^2}}{2}$$

### Chézy
 
 the equilibrium velocity
 $$0= - g h Z'_b - C_fu^2$$
 is the Chézy friction velocity $u_{chezy}=\sqrt{-g h Z'_b/C_f}$, which arises after the jump, where we have a long calm river.
 
 As the flow raet is constant $Q=u_{chezy}h$, so
 $h=\left(\frac{C_f Q}{g(-Z'_b)} \right)^{2/3}$
 
 
# Code 
 */
#include "grid/cartesian1D.h"
#include "saint-venant.h"

double tmax,Fr,deltah,x_rs,W,Zpb,Cf,uchezy,hmax,xF1=0,hF1=1;

int main(){
  X0 = 0;
  L0 = 5;
  N = 128*2;
  G = 1;
  tmax = 100;
  Fr = 2.5;    // Nombre de Froude initial
  Zpb = -.1;

  W = 0;  //vitesse du ressaut, ici 0
  x_rs=1.8;    // position nitiale du saut dans [X0;X0+L0]
  DT = 0.002;
    
 double CF[6]={.15, .17 , .18 , 0.2  , .25 , .3};  // loop on Froude
    
 for (int iF=0;iF< sizeof(CF) / sizeof(CF[0]);iF++){
        Cf = CF[iF];
        run();
        fprintf(stderr,"%lf %lf \n",Cf, xF1);
  }

}
/** 
sortie libre
*/
h[left] = dirichlet(1);
h[right] = h[];
u.n[left] =dirichlet(Fr);
u.n[right] = neumann(0) ;

/**
 initialisation: Belanger with $h_1=1$ and $h_2$ computed
 $$\big(\frac{h_2}{h_1}\big) =  
  \frac{-1+\sqrt{1+8 F_1^2}}{2}$$
 
*/
event init (i = 0)
{
  foreach(){
    zb[]=0.;
    // formule de Bélanger avec h1=1;
    deltah=(((-1 + sqrt(1+8*Fr*Fr))/2)-1);
     h[]=1+deltah*(1+tanh((x-x_rs)/.01))/2;
    // vitesse associé par conservation du flux Fr/h + translation 
     u.x[]=Fr/h[]+W;
      
      deltah = pow(sq( Fr*sqrt(G*1)*1)/(G*(-Zpb)/Cf),1./3);
      h[]= ((x-x_rs)<0?  1: deltah) ;
      
       //Q^2 = (hmax^3*G*(-Zpb)/Cf);
    }
  boundary({zb,h,u});
}

/**
### At each timestep
 
 We  use a simple explicit scheme for imposed slope ($-Z'_b$ ) implicit scheme to implement quadratic bottom
 friction with $C_f$. i.e.
 $$
 \frac{\partial \mathbf{u}}{\partial t} = -Z'_b g h  - C_f|\mathbf{u}|\frac{\mathbf{u}}{h}
 $$
  */
event friction (i++) {
    foreach() {
        double a = h[] < dry ? HUGE : 1. + Cf*dt*norm(u)/(h[]);
        foreach_dimension()
        u.x[] = (u.x[] - G*Zpb*dt)/a;
    }
    boundary ((scalar *){u});
}

/**
 we put a NON pysical longitudinal dissipation (as in Razis et al.  paper)
 $$ \frac{\partial h  \mathbf{u}}{\partial t}  = \nu_e \frac{\partial }{\partial x } (h \frac{\partial u}{\partial x } )$$
 */
 event friction_long (i++) {
   foreach_dimension() {
     face vector g[];
     scalar a = u.x;
     double nue=.0;
       foreach_face()
           g.x[] = nue*((h[] + h[-1,0])/2)*(a[] - a[-1,0])/Delta ;
       foreach ()
           if(h[]>dry){ u.x[] += dt/Delta*(g.x[1,0] - g.x[] + g.y[0,1] - g.y[])/((h[] + h[-1,0])/2);}
      boundary ((scalar *){u});
   }}


/**
 compute Chezy velocity associated to the slope
 and position of the jump for plots and analysis*/

event plot (t<tmax;t+=0.05) {
    hmax=0;
    foreach(){
        hmax=(h[]>hmax?h[]:hmax);
        xF1=(sq(u.x[])/(G*h[])<1?xF1:x);
        hF1=(sq(u.x[])/(G*h[])<1?hF1:h[]);
    }
    uchezy=sqrt(hmax*G*(-Zpb)/Cf);
   
#ifdef gnuX
    printf("set title 'Ressaut en 1D --- t= %.2lf '\n"
	 "p[%g:%g][-.5:5]  '-' u 1:($2+$4) t'free surface' w l lt 3,"
	 "'' u 1:3 t'velocity' w l lt 4,\\\n"
	 "'' u 1:4 t'topo' w l lt -1,\\\n"
     "'' u 1:(sqrt($3*$3/($2*%g))) t'Froude' w l lt 2,\\\n"
 //    "'' u 1:($2*%g+$3*$3/2) t'Charge' w l lt 5,\\\n"
     "'' u 1:( $2*(-1 + sqrt(1+8*$3*$3/%g/$2))/2) t'hj' w l lt 6,\\\n"
     "'+' u (%lf):(%lf) t'jump theo' , %lf , %lf \n",
           t,X0,X0+L0,G,G,xF1,hF1,uchezy,Fr*1/hmax);
    foreach()
    printf ("%g %g %g %g %g\n", x, h[], u.x[], zb[], t);
    printf ("e\n\n");
#endif
}
/** final plot
 */
#ifdef gnuX
#else
event end(t=tmax ) {
    foreach()
    fprintf (stdout,"%g %g %g %g %g\n", x, h[], u.x[], zb[], t);
}
#endif
/**
# Run and Results
## Run
To compile and run:

~~~bash
 qcc  -DgnuX=1  belangerturb.c  -lm; ./a.out 2>log | gnuplot
 
 make belangerturb.tst
 make belangerturb/plots;
 
 
 source  ../c2html.sh belangerturb
 
~~~

## Plot of jump
 
Plot of jump comparision 
with analytical result 
$$\big(\frac{h_2}{h_1}\big)   
 = \frac{-1+\sqrt{1+8 F_1^2}}{2}$$
 
 
~~~gnuplot  result, free surface (blue) and bottom (black)
 set xlabel "x"
 p [:1][1:1.2] 'out' t'free surface' w l lc 3,1+ 0.2*x lc 2,1+.25*x lc 1,1+.3*x
 
 
~~~

 
 
 
~~~gnuplot  result, free surface (blue) and bottom (black)
 set xlabel "x"
 p [:][-1:5] 'out' t'free surface' w l lc 3,'' u 1:3 t'speed' w l, \
   '' u 1:4 t'zb' w l lc -1
 
~~~

Position of shock
 
 ~~~gnuplot
 set xlabel "Cf"
 p [:][0:] 'log' t' x Jump' w lp
 
 ~~~

 
 
 
 
 
# Liens
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/belanger.c]() Cas normal
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/belangerdisp.c]() Cas dispersif
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/belangerturb.c]() Cas "turbulent"
 
 
# Bibliographie
 
* [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
"Equations de Saint Venant et application, Ecoulements en milieux naturels" Cours MSF12, M1 SU

 * Razis, Kanellopoulos, van der Weele,
 "Continuous hydraulic jumps in laminar channel flow"
 J. Fluid Mech. (2021), vol. 915, A8,
 
*/


