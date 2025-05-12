/** 
# Subcritical flow along a slope with a gaussian bump

Subcritical flow over a bump in a Nusselt film. For small bump we recover the Double Deck structure.

`MLSWbump13`: "MultiLayer Shallow Water", with a "bump", the "13" refers to the linearised solution in $(-ik)^{1/3}$
 in Fourier space.

## Code

We use the Multilayer solver of Saint-Venant
*/
#include "grid/cartesian1D.h"
#include "saint-venant.h" 

double tmax,Fr,hin,hout,u0,zbc,uzl;
double alpha,alphabump,alphamax;
scalar tau[],delta1[],delta2[],Qc[];
FILE * fq;
 
/**
position of domain `X0=0`, length `L0`, no dimension `G=1` 
run with *N* points. We define `hout` and `hin` from the analytical solution, 

The choice of  such that scales is $\nu=1$ and gravity balances the viscous term:
$$0 =  g \alpha + \nu \frac{\partial^2 u}{\partial z^2}$$
so that the basic profile is the Nusselt or Half Poiseuille profile:
$$ u(z)=\frac{\alpha g (2 h_{in} - z) z}{2 \nu}$$
The shear rate at the wall is 
v v v v v v v
$$ \frac{\partial u(z)}{\partial z}|_0=\frac{ \alpha g h_{in} }{3 \nu }$$
For one layer
$$0 =  g \alpha h_{in}  + \nu \frac{3 u_0}{h_{in}}$$
so that  we define $u_0 = (\alpha g/3/\nu) h_{in}^2$, the mass flow rate is $u_0h$ and 
shear rate at the wall is
 $$ \frac{\partial u(z)}{\partial z}|_0=\frac{ u0 }{h_{in}}$$

=============
$$ \frac{\partial u(z)}{\partial z}|_0=\frac{ \alpha g h_{in} }{ \nu }$$
The flow rate, as $\int_0^1 (2 - \zeta) \zeta d \zeta/2=1/3$, is say $u_0h$.
$$Q = \frac{\alpha g h_{in}^3}{3 \nu}= u_0h$$
so that   we define $u_0 = (\alpha g/3/\nu) h_{in}^2$ the mean velocity.  
The velocity at the surface $u_e=\frac{\alpha g (h_{in}^2)}{2 \nu} = \frac{3}{2}u_0$
$$ u(z)=\frac{\alpha g h_{in}^2 (2  - (z/h_{in}) (z/h_{in})}{2 \nu} $$
For one layer with Poiseuille closure indeed for a mass flow rate   $u_0h$, shear rate at the wall is
 $$ \frac{\partial u(z)}{\partial z}|_0=\frac{ 3 u_0 }{h_{in}}$$
and equilibrium between slope and friction holds:
$$0 =  g \alpha h_{in}  - \nu \frac{3 u_0}{h_{in}}$$
*************
$$ \frac{\partial u(z)}{\partial z}|_0=\frac{ \alpha g h_{in} }{ \nu }$$
The flow rate, as $\int_0^1 (2 - \zeta) \zeta d \zeta/2=1/3$, is say $u_0h$.
$$Q = \frac{\alpha g h_{in}^3}{3 \nu}= u_0h$$
so that   we define $u_0 = (\alpha g/3/\nu) h_{in}^2$ the mean velocity.  
The velocity at the surface $u_e=\frac{\alpha g (h_{in}^2)}{2 \nu} = \frac{3}{2}u_0$
$$ u(z)=\frac{\alpha g h_{in}^2 (2  - (z/h_{in}) (z/h_{in})}{2 \nu} $$
For one layer with Poiseuille closure indeed for a mass flow rate   $u_0h$, shear rate at the wall is
 $$ \frac{\partial u(z)}{\partial z}|_0=\frac{ 3 u_0 }{h_{in}}$$
and equilibrium between slope and friction holds:
$$0 =  g \alpha h_{in}  - \nu \frac{3 u_0}{h_{in}}$$
^ ^ ^ ^ ^ ^ ^

*/

int main() {
  int web=1;
  X0 = 0;
  L0 = .1-X0;
  G = 1;
  N = 128*2/web;
  tmax=12/web; 

  nu=1;
  alpha=1;   
  hout = 1;
  hin = 1;
  u0 = alpha*G/3/nu*hin*hin;
  Fr=u0/sqrt(G*hin);
/** 
Note that as here $u_0=1/3$ we  have a subcritical flow $Fr=1/3<1$.
First, the small bump case to compare with the "linearised Double Deck"
*/
  nl = 100/web ;
  alphamax=.5; 
  t=0;  
//  run();
//  system("cp xproff.txt xproff_L.txt");
/** 
Second, moderate bump case  with separation"
*/  
  nl = 100/web ;
  alphamax=2.25; 
  t=0; 

  run();
  system("cp xproff.txt xproff_M.txt");
  system("cp uproff.txt uproff_M.txt");
 

 // run();
 // system("cp xproff.txt xproff_M.txt");
 // system("cp xproff.txt xproff_M.txt");
/** third, moderate bump case  with only one layer!"
*/ 
  nl = 1;
  alphamax=2.25;  
  t=0;
  run(); 
  system("cp xproff.txt xproff_1.txt");
}
/** 
Initialisation,  we initialize `h` and `hc`.  
We define a scalar field `hc` to check for convergence on $h$. */
scalar hc[];
/**
 time 0, the bottom is the flat inclined plate */
event init (i = 0)
{
  foreach(){
    zb[] = - alpha*x ;
    delta1[]=0;
    h[] = hc[] =  hin ; 
   }
   boundary({zb,h,hc});

 /**
 initial velocity Poiseuille */
 int l = 0;
 foreach() { 
   if ( nl ==1) {
     u.x[] = u0;
   } else{
   double z = 0;
      l = 0;
      z -= (layer[0]/2)*hin;
      for (vector u in ul) {
       z += (layer[l++])*hin;
       u.x[] = alpha*G/nu*(2*hin - z)*z/2;
      }
     }
    }

/**
Boundary condtion, neumann for $h$ and impose  flux 
   and neumann at the exit, we impose the output position of the free surface.
*/
  
  h[left]   =  neumann(0.);  
  eta[left] =  neumann(0.);
  h[right] =  neumann(0.);
  eta[right] = dirichlet(zb[] + hout);  
/**
 left velocity?  */ 
  zbc=0; 
  l = 0;
   for (vector uz in ul) {
     zbc += (layer[l++]); 
     uzl =(alpha*G/nu*(2*1 - zbc)*zbc/2);
     uz.n[left] =  dirichlet(u0/h[left]);
     //dirichlet(uzl*h[left]*h[left]);//
     // dirichlet(u0/h[left]); //  dirichlet(u0/h[left]);// dirichlet(u0);//  dirichlet(uzl*h[left]*h[left // 
     uz.n[right] = neumann(0.);
   // fprintf (stderr,"coucou %d u=%g z=%g\n",l, uzl,zbc);
  }
  
}
/**
 We update the bump */ 
event updatebump (i++) {
  alphabump=alphamax*(1-exp(-t)) ;
  double eps=.15,xc=0.05;  
  foreach()   
      zb[] = - alpha*x + alphabump*eps*exp(-sq((x-xc)/eps/eps/eps)) ;
  } 
/**
 Compute friction, 
*/ 
event friction (i++) {
//   fprintf (stderr,"coucou %d \n",i);
/**
 Poiseuille friction: implicit scheme (time splitting, here commented)
 $$\frac{\partial u}{\partial t} = - 3 \nu \frac{u}{h^2}$$
 hence, if implicited for $u$, explicited for $h$:
 $$ u^{n+1} = u^n -3 \nu \Delta t \frac{u^{n+1}}{(h^{n})^2}$$
*/ 
 if(nl==1){ 
  foreach()
  {    double ff = h[] < dry ? HUGE : (1. +   3*nu*dt/h[]/h[]);
    u.x[] /= ff;
    }
   foreach()
    tau[]=3*u.x[]/(h[]); 
  }
/**
Multilayer case, the derivative is twice the value in the first layer, 
*/   
else{
  vector u0 = ul[0];
  foreach()
    tau[]=2*u0.x[]/(h[]/nl);
   }
}
/**
We check for convergence. */
event logfile (t += 0.1) {
  double dh = change (h, hc);
  fprintf (stderr,"coucou t= %g  -- var=  %g,  alpha = %lf u0 = %g \n", t, dh, alphabump, u0);
  if (i > 0 && dh < 1e-6){
   fprintf (stderr,"conv\n");  
    return 1;
  }
}
/** 
And save fields, 
*/
event sauv (t<=tmax;t+=.01) { 
/** compute the flux
and displacement thickness  and velocity at the surface

*/
 vector ue;
 foreach() { 
      double z = zb[]*0,psi=0,Qm=0;
      ue=ul[nl-1];
      int l = 0;
      z -= (layer[0]/2)*h[];
      for (vector uz in ul) {
       z += (layer[l])*h[];  // was  [l++]
       psi += (1-uz.x[]/ue.x[])*(layer[l])*h[];
       Qm += uz.x[]*(layer[l])*h[];
       l++;
    if(abs(Qm>10)) {fprintf (stderr,"%g %g %g %g %g\n",x,z,Qm,h[],layer[l]); 
         getchar();} //debug
      }
    if ( nl ==1) psi = 1./3.;    
  delta1[]=  psi; 
  Qc[] = Qm;//(h[]-psi)*ue.x[]; 
}
 /** save profiles along $x$ 
*/
 FILE *fq =  fopen("xproff.txt","w");  
 foreach() { 
      fprintf(fq,"%g %g %g %g %g %g %g %g %g %g %g\n", x, h[], tau[], ue.x[],t,Fr,nu,delta1[],zb[],eta[],Qc[]);}
  fprintf(fq," \n");
  fflush(fq);
  fclose(fq);

 /** save profiles along $x$ and $z$ of velocity  
*/
if (nl>1){
FILE *fp =  fopen("uproff.txt","w"); 
 foreach() { 
   double z = zb[]*0;
      int l = 0;
      z -= (layer[0]/2)*h[];
      for (vector u in ul) {
       z += (layer[l++])*h[]; 
       fprintf(fp,"%g %g %g %g \n", x, z, u.x[], h[]);
      }
      if ( nl ==1)  u.x[] = u0;
   fprintf(fp,"\n");
}
  fflush(fp);
  fclose(fp);
}
}

/**
 
## Run
To compile and run  :

~~~bash
   qcc -O2 -o MLSWbump13 MLSWbump13.c -lm; ./MLSWbump13
~~~
 
using `Makefile`

~~~bash
  make MLSWbump13.tst
  make MLSWbump13/plots
  make MLSWbump13.c.html 
  open MLSWbump13.c.html
~~~  

## Plots

Note no upstream influence (perturbation appear at the bump foot), 
boundary layer separation is after the bump only for multilyer case.
Shape of the bump and the soil $z_b$, depth $h$ (multilayer),  friction   $\tau$ (multilayer), and as well  
depth $h_1$ (one layer),  friction   $\tau_1$ (one layer) 
 
~~~gnuplot
 set xlabel "x"
 alpha=1
 nu=1 
 u0=1./3
 #p[][:7]'xproff.txt' t'h'w l,''u 1:($11/u0) t'Q' w l,''u 1:($4/u0) t'ue' w l,''u 1:3 t'tau' w lp,\
 #'' u 1:($9+$1)t'bump' w l linec -1
 p[][:7]'xproff_M.txt' t'h'w l,''u 1:($11/u0) t'Q' w l,''u 1:($4/u0) t'ue' w l,''u 1:3 t'tau' w lp,\
 '' u 1:($9+$1)t'bump' w l linec -1,'xproff_1.txt' t'h1'w l,''u 1:3 t'tau1' w l
 ~~~

velocity profiles in the entry region and in around the bump + half Poiseuille

~~~gnuplot
 reset
set ylabel "y"
 p(x,y)=abs(x-y)<.01 
 u(x)=x*(2-x)/2
 #p'uproff.txt'u ($3):(p($1,.0025)?$2:NaN) t 'e',''u ($3+.5):(p($1,.05)?$2:NaN)t 'b',\
 #''u (u($2)):2 t'p' w d,''u (u($2)+.5):2 t'p' w d
 p'uproff.txt'u ($3):(p($1,.0025)?$2:NaN) t 'e',''u ($3+.5):(p($1,.05)?$2:NaN)t 'b',\
 ''u (u($2)):2 t'p' w d,''u (u($2)+.5):2 t'p' w d
~~~


The double deck solution is in Fourier space for each mode $k$ of the fourier transform of the bumb $f_k$:
the linearized solution is $\tau = U'_0 +  \alpha TF^{-1} [\tau_k]$ with $\alpha \ll 1$, 
for each mode $k$ in Fourier Space ($f_k =TF[f]$):
$$\tau_k  =  U'_0 (3 Ai(0)) (-i  k U'_0)^{1/3} f_k,$$
 
   $$p_k = U'_0^2 (3 Ai'(0)) (-i k U'_0)^{-1/3} f_k$$

                              
~~~gnuplot
 p[0.0:0.07][:2]'xproff.txt' t'h'w l,''u 1:3 t'tau' w lp,'' u 1:($9+$1)t'bump' w l 
~~~

bla

~~~gnuplot
 
$$p_k = U'_0^2 (3 Ai'(0)) (-i k U'_0)^{-1/3} f_k$$
 
$$p_k = U'_0^2 (3 Ai'(0)) (-i k U'_0)^{-1/3} f_k$$

~~~gnuplot
 set xlabel "x"
 eps=.15
 alphabump=0.5*epx
 p[0.0:0.07][:2]'xproff.txt' t'h'w l,''u 1:3 t'tau' w lp,'' u 1:($9+$1)t'bump' w l linec -1,\
  'img/lin13.txt' u ($1*(eps**3)+0.05):((1+0.5*eps*$2*(1)**(4./3)/eps))w l t'lin'
~~~
 
~~~gnuplot
 set xlabel "x"
 eps=.15
 alphabump=0.5*eps
 p[0.0:0.07][:2]'xproff_L.txt' t'h'w l,''u 1:3 t'tau' w lp,'' u 1:($9+$1)t'bump' w l linec -1,\
  '../Img/lin13.txt' u ($1*(eps**3)+0.05):((1+0.5*eps*$2*(1)**(4./3)/eps))w l t'lin' 
~~~

~~~gnuplot
 set xlabel "x"
 eps=.15
 alphabump=0.5*eps
 p[0.0:0.07][:2]'xproff_L.txt' t'h'w l,''u 1:3 t'tau' w lp,'' u 1:($9+$1)t'bump' w l linec -1,\
  '../Img/lin13.txt' u ($1*(eps**3)+0.05):((1+0.5*eps*$2*(1)**(4./3)/eps))w l t'lin'
~~~
 
bla

                              
=============
Check that $\frac{u_e^2}{2}+ p$ varie peu en entrée (with $p=G(h+Z_b)$) donc $u_e \partial_x u_e + \partial_x p=0$,
 puis ensuite effectivement $\partial_x p =cst$.

~~~gnuplot
 p[][:]'xproff_L.txt' t'h'w l,''u 1:($4*$4/2+($2+$9)) t'Ber' w l
~~~

Check that displacement thickness starts from Blasius and goes to Poiseuille:

~~~gnuplot
 u0=1./3
 p[][:]'xproff_L.txt' u 1:8,.33,1.7*sqrt(x/u0)
~~~


*************
Check that $\frac{u_e^2}{2}+ p$ varie peu en entrée (with $p=G(h+Z_b)$) donc $u_e \partial_x u_e + \partial_x p=0$,
 puis ensuite effectivement $\partial_x p =cst$.

~~~gnuplot
 p[][:]'xproff_L.txt' t'h'w l,''u 1:($4*$4/2+($2+$9)) t'Ber' w l
~~~

Check that displacement thickness starts from Blasius and goes to Poiseuille:

~~~gnuplot
 u0=1./3
 p[][:]'xproff_L.txt' u 1:8,.33,1.7*sqrt(x/u0)
~~~

 
*Grenoble 02/11/2018* 

## Bibliography

* [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/CISM/tuyau_CISM.pdf) Double Deck (pipe flows)

* [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
"Ecoulements en milieux naturels" Cours MSF12, M2 SU

* [F. Chouly Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/choulylagree12.pdf) 
 Comparison of computations of asymptotic flow models in a constricted channel


*/
