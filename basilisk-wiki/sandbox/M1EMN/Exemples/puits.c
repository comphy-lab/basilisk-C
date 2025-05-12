/**
#    Ideal fluid "source/ well",


The classical Shallow Water Equations in 1D with no friction but with a source of magnitude $q_s$ ( if $q_s<0$ it is a well).
 The source term is $q_s\delta(x-x_s)$ where $x_s$ is the position of the object and $\delta$ the dirac distribution.
 Her we ar in 1D, so that $\int_{-\infty}^{+\infty}\delta(x-x_s)=1$
$$
 \left\{\begin{array}{l}
 \partial_t h+\partial_x Q=q_s \delta(x-x_s)\\
	\partial_t Q+ \partial_x \left[\dfrac{Q^2}{h}+g\dfrac{h^2}{2}\right]
= 0
        \end{array}\right. 
$$
with at  $t=0$, a lake on the rest  $h(x,t=0)=2$.
As the flow is not viscous, the convective term is important.
*/
#include "grid/cartesian1D.h"
#include "saint-venant.h"

double tmax,q_puits,x_puits,Dx_puits;

u.n[left]   = dirichlet(0);
u.n[right]  = dirichlet(0);

/**
 position of domain `X0`, length `L0`, no dimension `G=1`,  
run with 512 points (less is not enough)
*/
int main() {
  X0 = -5.;
  L0 = 10.;
  G = 1;
  N = 512*2;
  tmax=30;
  
  run();
}
/**
 start by a lake, the soil is flat $z_b=0$
 */
event init (i = 0)
{
    foreach(){
        zb[] = 0;
        h[] =  2;
        u.x[] =0;
    }
    boundary({zb,h,u});
#ifdef gnuX
    printf("\nset grid\n");
#endif
}

/**

 Here the well it self:
 We solve the problem by splitting:
 $$\frac{h^*-h^n}{\Delta t} + \frac{\partial Q^n}{\partial x} =0, \text{ then } \frac{h^{n+1}-h^*}{\Delta t} =  q_s \delta(x-x_s)$$
 The dirac is in practice approximated by a rectangle of hight $1/\Delta_{xpuit}$ and width $\Delta_{xpuit}$,
 we take by inspection  $\Delta_{xpuit}=2 \Delta$ (a multiple of $2 \Delta$ is fine.
 
  $$h^{n+1}=h^* +  \Delta t q_s \delta(x-x_s)$$
 
 and as we add no $Q$
$$\frac{\partial}{\partial t}V =\frac{\partial}{\partial t}\int h dx
 = \int_{-x_0}^{x_s^-} \frac{\partial h }{\partial t} + q_s + \int_{x_s^+}^{x_0} \frac{\partial h }{\partial t}=
  -[[Q]]_{-x_0}^{x_s^-}   +  q_s -[[Q]]_{x_s^+}^{x_0} = q_s$$

*/
event puits(i++){
    q_puits = -.125;
    x_puits=0;
    foreach(){
        Dx_puits=2*Delta;
        h[]+=(fabs((x-x_puits))<=Dx_puits/2 ? q_puits*dt/Dx_puits  :0 );
        h[]=max(h[],0);
    }
}
/**
 monitoring of
  $h(X_0/2,t)$,  Velocity $u(X_0/2,t),  $\partial_t h(X_0/2,t)$,  Volume $V= \int h dx$ and $dV/dt$
*/
event mesure(i++){
    static double hold=1,dhdt=0,dVdt=0,Vold=0,V=0;
    V=0;
    foreach()
       V+=h[]*Delta;
    if (t>.01){
        dhdt=(interpolate(h,X0/2,0)-hold)/dt;
        dVdt=(V-Vold)/dt;
        fprintf(stderr," %g %g %g %g %g %g\n",t,interpolate(u.x,.5*X0,0),interpolate(h,.5*X0,0),dhdt,V,dVdt);}
    hold=interpolate(h,X0/2,0);
    Vold=V;
}
/**
Output in gnuplot if the flag `gnuX` is defined
*/
#ifdef gnuX
event plot (t<tmax;t+=0.01) {
    printf("set title 'Ressaut en 1D --- t= %.2lf '\n"
      "h(x)=(((x-0)<-t)+((x-0)>-t)*(2./3*(1-(x-0)/(2*t)))**2)*(((x-0)<2*t));t= %.2lf  ; "
      "p[%g:%g][-.5:2]  '-' u 1:($2+$4) t'free surface' w l lt 3,"
           "'' u 1:($1<0?$3/($1+5):NaN) t'Q/x' w l lt 4,\\\n"
      "'' u 1:4 t'topo' w l lt -1,\\\n"
      "'+' u (0):(0.296296296296296) t'theo Qmax', h(x) t 'theo h(x)'\n",
           t,t,X0,X0+L0);
    foreach()
    printf ("%g %g %g %g %g\n", x, h[], u.x[], zb[], t);
    printf ("e\n\n");
}
/**
Output at the end if not defined
*/
#else
event end(t=tmax ) {
    //foreach()
    // printf ("%g %g %g %g %g\n", x, h[], u.x[], zb[], t);
}
#endif
/**

## Run
To compile and run with gnuplot:

~~~bash
 qcc -DgnuX=1   -O2   -o puitX  puits.c -lm;  ./puitX | gnuplot
 
 make puits.tst
 make puits/plots
 make puits.c.html
 
 source ../Exemples/c2html.sh puits

~~~

## Plots
 
 
 Total volume as function of time, the slope is indeed `q_puits` and $dV/dt$ is indeed `q_puits` so that we plot $V0+(dV/dt)t$ to control
 
~~~gnuplot
set xlabel "t"
 p[][0:] 'log' u 1:5 w lp t'V' ,'' u 1:($1<30?(20+$6*$1):NaN)  t'V0+(dV/dt)t'w l,20-.125*x
 
~~~
 
Plot of $h(X_0/2,t)$, we notice the stair steps due to the reflection of the signal at the wall.
 Velocity $u(X_0/2,t)$ is positive, it corresponds to the release of water through the well, wave by wave.
 As time increases, the wave velocity is smaller and smaller ($\sqrt{gh}$) that is why the steps are longer and longer...

 
~~~gnuplot
 set xlabel "t"
 p[0:][-.5:2.5]'log' u 1:2 w l t'u(x0/2,t)','' u 1:3 t'h(x0/2,t)' w l,'' u 1:(2+$6*$1/10) w l t'(V/10)',.25/2/2, ''u 1:(2-$3) w l

 
 
~~~

# Remarks
 
 * playing with the width of the source
 
 
 
 * A linearised theory or expansion wave theory can be done to estimate the velocity and the stairs steps....
 I guess that the step is $q_{puits}/2$ (one half right one half left), the velocity is $\sqrt{g h}q_{puits}/2$, the length of the step is $x_0/\sqrt{g h}$...
 
 * In we increase the flow rate, the $dV/dt$ is finaly no more linear?
 
 There is a regime where $\frac{\partial h}{\partial t} \sim q_0$ constant in $x$.
 All teh level decreases uniformly in $x$.  Hence, by mass conservation we estimate:
 
 $Q \simeq - q_0 (x+X_0)$ for $x<0$ and
 $Q = - q_0 (x-X_0)$  for $x>0$.
 
 The flux rate $Q$ is linear in $x$ with a discontinuity in 0: we see a "Z"!!!
 
 
 
 
# Links
 
 
 see non viscous dam break with [standard C](http://basilisk.fr/sandbox/M1EMN/Exemples/svdb.c) 
and with [Basilisk](http://basilisk.fr/sandbox/M1EMN/Exemples/damb.c)

# Bibliographie
 
* [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
"Equations de Saint Venant et application, Ecoulements en milieux naturels" Cours MSF12, M1 UPMC





 Version 1: confine: wash your hands in the well

*/

   

 
 
 
