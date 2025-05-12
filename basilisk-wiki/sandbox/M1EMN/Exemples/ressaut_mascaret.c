/**
# From hydrolic jump to undular bores,  "du ressaut au mascaret"

 
 
## The problem

 
 How an hydrolic jump moves and is deformed by dispersion?
 
 
 This is done  with the Green-Naghdi equations (see switch `SV`) and with Saintt-Venant with Boussinesq correction (and here with turbulent friction).
 

The left level of water is  of elevation $h_1$, velocity is $u_1$, the Froude number is $F_1=u_1^2/(gh_1)$, this jump goes on a level of depth $h_2$ of still water ($u_2=0$).

We want to see how a hydrolic jump changes into an ondular bore.

 
 
## Dispersive terms

 The full system with Boussinesq correction ($\frac{h^3}{3}\frac{\partial }{\partial t}\frac{\partial^2 u}{\partial x^2}$) and friction ($- C_f |Q|\frac{Q}{h^2}$)  reads:
 $$\frac{h^{n+1}-h^n}{\Delta t}+ \frac{\partial Q}{\partial x}=0$$
 $$\frac{Q^{n+1}-Q^n}{\Delta t}+\frac{\partial }{\partial x}(\frac{Q^2}{h}+\rho g \frac{h^2}{2})=\frac{h^3}{3}\frac{\partial }{\partial t}\frac{\partial^2 u}{\partial x^2} - C_f |Q|\frac{Q}{h^2}$$
 
 we do here an extra splitting,
 we estimate first Saint Venant estimate $h^*$ and $Q^*$ without dispersion, next we add the extra dispersive term.
 
 
 We can compare with the Serre Green-Nagdhi, we will compare with multilyer latter.
 
 
## Code

 if `SV` defined, then we use the Saint-Venant model,
 else we use Serre-Green-Naghdi model.
 
We define an extra switch `DISP`, when it is defined, we solve the Saint-Venant model plus an extra term (Boussinesq) whose origin is explained thereafter
*/

#undef SV
#define SV
#define DISP


#include "grid/multigrid1D.h"
#ifdef SV
#include "saint-venant.h"
#else
#include "green-naghdi.h"
#endif
/**
 extrapolation of the Integral of Airy function
*/
double IntgAiry(double x)
{ double ai;
    if((x>-4)&&(x<3)){
    ai=-0.47800749642926166 + (4 + x)*(0.11541826749089747 + (-3 + x)*
     (-0.0013612022054726482 + x*(-0.03341069433229331 + (2 + x)*
    (0.006480518707412483 + (-2 + x)*(0.003011122308321213 - 0.0008866964376086442*(3 + x))))));}
    else{ai=0/0.1;}
  return ai;
}
/**
Neumann conditions
*/
u.n[left] = neumann(0);
h[left] = neumann(0);
u.n[right] = neumann(0);
h[right] = neumann(0);
double h0,h1,h2,W,u1,u2,U1,U2,tmax;
scalar F[],f[];
F[left] = neumann(0);
f[left] = neumann(0);
/**
Main and parameters

 The domain is 75 long, centered on the origin.
 The problem is without dimensions. The problem is
 solved in one dimension with 2048 grid points.
*/
int main()
{
  X0 = -20;
  L0 = 75;
  G = 1.;
  h0 = 1;
  h1 = h0;
  h2 = h0-.1;
  DT = HUGE;
#ifdef SV
   N = 1024;
#else
   gradient = NULL;
   N = 1024;
#endif
    tmax = 45;  // L0=75 512 45
  run();
}
/**
 The initial conditions are 
 $h(x<0,t=0)=h_1>h_2$ and  $u(x<0,t=0)=u_1$ 
 
 
 and for positive $x$:
  $h(x>0,t=0)=h_2<h_1$  and $u(x>0,t=0)=u_2$.
  
  
 we define
 $U_1= -\sqrt{ g (\frac{h_1+h_2}{2})\frac{h_2}{h_1}}$ and $U_2= -\sqrt{ g (\frac{h_1+h_2}{2})\frac{h_1}{h_2}}$.
 They correspond to the standing jump:
 if $u_1=U_1<0$ (subcritical $|u_1/\sqrt{g h_1}|<1$) and $u_2=U_1<0$ (supecritical $|u_2/\sqrt{g h_2}|>2$), then $W=0$.
 
 If we do a moving jump, we impose any value to $u_1$, then the velocity after the shock is $u_2=U_2 + (u_1 -U_1)$ and the shock velocity is
 $W=(u_1 -U_1)$.  Here we take as value of $u_1=U_1  - U_2$ in order to have a flow at rest after the shock $u_2=0$,
 the  shock velocity is $W=-U_2>0$.
 
*/

event init (i = 0){
   foreach(){
     h[] = (h1+h2)/2 + (h2-h1)/2*tanh(x) ; // x < 0. ? h1 : h2;
     U1 = -sqrt(G/2 *(1+h2/h1)*h2);
     U2 = -sqrt(G/2 *(1+h1/h2)*h1);
     u1 = U1 + (- U2) ;
//      u1 = U1;   // uncomment for standing jump
     W = 0 + (u1 -U1);
     u2 = U2 + (u1 -U1);
     u.x[]= (u1+u2)/2 + (u2-u1)/2*tanh(x) ; // x < 0. ? u1 : u2;
     zb[] = 0.;
   }
}
/**
 old value usefull for dispersion
*/
scalar un[];
event init_u (i = 0) {
    foreach() {
        un[] = u.x[];
    }
}

#ifdef DISP
/**
## Poor's man dispersive model
 Let say that pressure is the sum of the hydrostatic pressure plus a perturbation of it:
 $p= \rho g h + \Pi$. This perturbation is linked to the up to now neglected transverse momentum:
 $-\frac{\partial \Pi}{\partial y} \simeq \rho \frac{\partial v}{\partial t}$,
 as we linearise around a mean level of water say $h_0$.
 Hence, as $v \simeq - y \frac{\partial u}{\partial x}$, we subsitute in transverse momentum and integrate $-\frac{\partial \Pi}{\partial y}$
 to obtain  $\Pi \simeq \rho \frac{y^2-h_0^2}{2} \frac{\partial^2 u}{\partial t \partial x }$
 we estimate the extra pressure gradient that will modify the longitudinal momentum with   $\Pi(h)=0$ 
 $$ \int_0^y  \Pi dy \simeq   \rho [(\frac{y^3}{6}  - \frac{y h_0^2}{3})]_0^{h_0}\frac{\partial^2 u}{\partial t \partial x}$$
 
 $$ \int_0^{h_0}  \Pi dy \simeq -  \rho (\frac{h_0^3}{3})\frac{\partial^2 u}{\partial t \partial x}$$
so that 
 $$- \int_0^h\frac{\partial \Pi}{\partial x} dy \simeq \frac{1}{3} \rho h_0^2 \frac{\partial^3 u}{\partial t \partial x^2}$$
 the coefficient nicely allows to reobtain   the dispersion relation
 $\omega^2 =gk \tan k h_0$ at small $kh_0$ (up to order 3).
 
 The non hydrostatic part of the pressure gives an extra contribution to Saint Venant
 solver $$\alpha h^3 \frac{\partial }{\partial t}\frac{\partial^2 u}{\partial x^2}$$
 we do here an extra splitting
 we estimate first Saint Venant estimate $h^*$ and $Q^*$:
 $$\frac{h^{*}-h^n}{\Delta t}+ \frac{\partial Q}{\partial x}=0$$
 $$\frac{Q^{*}-Q^n}{\Delta t}+ \frac{\partial  }{\partial x}(\frac{Q^2}{h}+ \rho g \frac{h^2}{2}) = 0$$
 then as $u^*=Q^{*}/h^*$ and $u^n=Q^{n}/h^n$, the longitudical momentum equation  may be written as  $u^{*}-u^n = F$
 and in the new source term $(\alpha h^3 \frac{\partial }{\partial t}\frac{\partial^2 u}{\partial x^2})\Delta t$
 the variation $\frac{\partial u}{\partial t}\Delta t$ is replaced by the total  increase
 $(u^{n+1}-u^n)$   so that the source  is:
 $\alpha h^2  \frac{\partial^2  }{\partial x^2}(u^{n+1}-u^n)$
 then, the splitting reads:
 $$u^{n+1}-u^n = F + \alpha h^2 \frac{\partial^2}{\partial x^2}(u^{n+1}-u^n)$$
  it allows to compute implicitely the increase in time of $u$
      $$(1 -  \alpha h^2  \frac{\partial^2  }{\partial x^2}) (u^{n+1}-u^n) = F$$
*/

event disperssse(i++){
/**
 Store old value and estimate of  explicit variation of velocity after Saint Venant splitting  $u^{*}-u^n = F$, the corrected
 velocity increase is then $u^{n+1}-u^n =f$ where $f =\left((1-\frac{h^2}{3}\frac{\partial^2  }{\partial x^2})\right)^{-1}F$
 
*/
    scalar uo[];
    foreach(){
      uo[] = un[];
      F[] =  (u.x[] -uo[]);
      f[] = F[];
    }
   
    double dh = change (u.x, un);
    if ( (i > 1 && dh < 1e-7) || i==10000000) {
        foreach()
        fprintf (stderr, "%g %g\n", x, h[]);
        return 1; /* stop */
    }
    /** compute total increase of $u$
    by inversion of
     $(1 -
     \frac{h^2}{3} \frac{\partial^2 }{\partial x^2 } )  f  = F$
     */
    double eps=.0,fn=0.,omega=.25,s;
    do
    { eps=0;
        foreach()
        {   s =   sq(h[])/3/sq(Delta);
            fn = (F[] + s*(f[-1] + f[1]))/(1+ 2*s);
            eps = eps + sq(f[]-fn);
            f[] = omega * fn + (1-omega) * f[] ;
        }
    }while(sqrt(eps)>1.e-15);
 
    /**
     we update the velocity
     $u(x, t+\Delta t) = u (x, t) + f$
     */
    foreach(){
     u.x[] = uo[] + (f[]);
     un[] = u.x[];
    }


}
#endif
/**
 We use a simple implicit scheme to implement quadratic bottom
 friction by splitting i.e.
 $$
 \frac{d\mathbf{u}}{dt} = - C_f|\mathbf{u}|\frac{\mathbf{u}}{h}
 $$
 with $C_f=.05$. */
event friction(i++){
    foreach() {
        double a = h[] < dry ? HUGE : 1. + .03*dt*norm(u)/h[];
        foreach_dimension()
        u.x[] /= a;
    }
}
/**
 output
*/
event output1(t+=.5){
    static FILE * fp = fopen("maxht.OUT","w");
    double hmax;
    hmax=0;
    foreach()
     if(h[]>=hmax)hmax=h[];
    fprintf (fp, " %g %g \n",t,hmax);

}


#ifdef gnuX
event output  (t += .1; t <= tmax) {
    fprintf (stdout, " reset \n set title \"t=%g\"\n",t);
    fprintf (stdout, " set label \"s\" at %g,%g \n",W*t,h0);
    fprintf (stdout, "p[:][:]  '-' u 1:2  not  w l linec 3, '' u 1:3  w l\n");
#else
  event output (t += 1; t <= tmax) {
#endif 	
	
  foreach()
    fprintf (stdout, "%g %g %g %g %g \n", x , h[],u.x[], (x- W*t)/pow(t+dry,1./3),t);
#ifdef gnuX
    fprintf (stdout, "e\n\n");
#else    
    fprintf (stdout, "\n");     
#endif 
}

event outputf (t = end) {
    foreach()
      fprintf (stdout, "%g %g %g %g %g \n", x , h[],u.x[], (x- W*t)/pow(t+dry,1./3),t);
    foreach()
      fprintf (stderr, "%g %g %g %g %g \n", x , h[],u.x[], (x- W*t)/pow(t+dry,1./3),t);
    }
/**
 
## compilation
 
~~~bash
qcc -g -O2 -DTRASH=1 -Wall  -DgnuX=1  -o ressaut_mascaret ressaut_mascaret.c -lm
./ressaut_mascaret 2> log | gnuplot
 
 
make ressaut_mascaret.tst
make ressaut_mascaret/plots
make ressaut_mascaret.c.html


source c2html.sh ressaut_mascaret
 
~~~

## Results
 
The standing jump is dispersed, we see the formation of the ondular bore.

~~~gnuplot Fluid depth profile.
set xlabel 'h'
set ylabel 'z'
 plot [0:][0:1.5]'out' u 1:($5==5?$2:NaN)  w l t 't=05',\
 ''  u 1:($5==10?$2:NaN) w l t 't=10',\
 ''  u 1:($5==20?$2:NaN) w l t 't=20',\
 ''  u 1:($5==30?$2:NaN) w l t 't=30',\
 ''  u 1:($5==40?$2:NaN) w l t 't=40'
~~~

 
~~~gnuplot max wave elevation
 set xlabel 't'
 set ylabel 'hmax'
 plot [0:][0:1.5]'maxht.OUT' u 1:2  w l
 ~~~
 
~~~gnuplot
iai(x)=((x<3)&&(x>-4)?-0.478+(4+x)*(0.115418+(-3+x)*(-0.00136120+x*(-0.033411+(2+x)*(0.006480+(-2+x)*(0.0030111-0.0008867*(3+x)))))) : NaN)
iai(x)=((x<3)&&(x>-4)?-0.47800749642926166 + (4 + x)*(0.11541826749089747+(-3+x)*(-0.0013612022054726482+x*(-0.03341069433229331+(2+x)*(0.006480518707412483 + (-2+x)*(0.003011122308321213-0.0008866964376086442*(3+x)))))): NaN)
p[-3:2]'out' u ($4):($5>10?$2:NaN) w lp #,''u ($4):(-iai($4*4)+.33)*3/2.*.1+.9
~~~
 
~~~gnuplot
 reset
 set pm3d map
 set palette gray negative
 unset colorbox
 set tmargin at screen 0.95
 set bmargin at screen 0.15
 set rmargin at screen 0.95
 set lmargin at screen 0.15
 set xlabel "x"
 set ylabel "t"
 unset key
 splot [][][.85:1.1]'out' u 1:5:2
~~~
 
 
## Links
 
* [http://basilisk.fr/sandbox/M1EMN/BASIC/disperse.c]()
* [http://basilisk.fr/sandbox/M1EMN/Exemples/belanger.c]()
 
 
# Bibliography
 
 * [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
 "Equations de Saint Venant et application, Ecoulements en milieux naturels" Cours MSF12, M1 UPMC

 
 
 
OK v2
*/ 
