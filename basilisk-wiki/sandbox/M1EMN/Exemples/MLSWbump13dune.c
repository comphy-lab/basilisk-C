/** 
# Subcritical flow along a slope with a dune

Subcritical flow over a bump in a Nusselt film. For small bump we recover the Double Deck structure.

`MLSWbump13`: "MultiLayer Shallow Water", with a "bump", the "13" refers to the linearised solution in $(-ik)^{1/3}$
 in Fourier space.
 Furthermore her we add erosion and sedimentation so that the dune moves.
 the initial slope is not erodible, so that we do not digg below 0.
 

 
 
 
~~~gnuplot hydrolic jump configuration
 set arrow from 2,4 to 3,1.5 front
 set label  "g" at 2.5,2 front
 set label  "gravity tilted" at 2.5,4 front
 set arrow from 4,0 to 4,1 front
 set arrow from 4,1 to 4,0 front
 
 set xlabel "x"
 set ylabel "h"
 p [0:4][0:1.2]'xproff.txt' u ($1*10):10 t'dune ' w l
 reset
~~~
 
 
 
 
 
## Code: fluid part

We use the Multilayer solver of Saint-Venant
*/
#include "grid/cartesian1D.h"
#include "saint-venant.h" 

double tmax,Fr,hin,hout,u0,zbc,uzl;
double alpha,alphabump,alphamax,masse,c0;
double lambda=0.0,V=1200.,taus=.9,Qt=.005;
scalar tau[],delta1[],delta2[],Qc[];
scalar qs[],qsx[];
scalar Ar[],Arm[],taud[],zbump[];
FILE *fq;
 
/**
position of domain `X0=0`, length `L0`, no dimension `G=1` 
run with *N* points. We define `hout` and `hin` from the analytical solution, 

The choice of  such that scales is $\nu=1$ and gravity balances the viscous term:
$$0 =  g \alpha + \nu \frac{\partial^2 u}{\partial z^2}$$
so that the basic profile is the Nusselt or Half Poiseuille profile:
$$ u(z)=\frac{\alpha g (2 h_{in} - z) z}{2 \nu}$$
The shear rate at the wall is 
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

*/

int main(int argc , char * argv[]) {
  if (argc>1){fprintf (stderr,"coucou  z=%g\n",atof(argv[1])); V=atof(argv[1]); }
//int web=1;
  X0 = 0;
  L0 = (.5/2-X0)*2;
  G = 1;
  N = 1024/2;
  tmax=5;

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
  nl = 50/2;
  alphamax=1; 
  t=0; 
  system("rm dune.txt;touch dune.txt"); 
  run();
   FILE *fc =  fopen("vm.dta","a");  

  //fprintf (fc,"%g %g\n",V*pow(masse,.75),c0*0.005*V*pow(masse,.25));  ; 
  fprintf (fc,"%g %g %g\n",masse,V,c0);  ; 
  fflush(fc);
  fclose(fc);
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
Boundary condition, neumann for $h$ and impose the flux
   and neumann at the exit, we impose the output position of the free surface.
*/
  
  h[left]   =  neumann(0.);  
  eta[left] =  neumann(0.);
  h[right] =  dirichlet(  hout); // neumann(0.);
  eta[right] =  dirichlet(zb[] + hout);


Arm[left]=dirichlet(0);
Ar[left] =dirichlet(0);
qs[left] = dirichlet(0);


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
 Compute friction, 
*/ 
event friction (i++) {
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
  fprintf (stderr,"coucou t= %g  -- var=  %g,  alpha = %lf masse = %g c0 = %g \n",
    t, dh, alphabump, masse, c0);
  if (i > 0 && dh < 5e-6){
   fprintf (stderr,"conv\n");  
    return 1;
  }
}



/**
## Code: solid part

-------------- Calcul du nouveau fond au pas de temps suivant -------------
 
 Here is the change of the topography due to erosion and sedimentation.
 
 
 
First, we compute the flow rate of granular in the bed load:
 $$\frac{\partial}{\partial x }q_s  + V q_s=  Ar$$
where extraction term
 $$Ar = V (\tau - \Lambda \frac{\partial f}{\partial x}- \tau_s)_+$$
note the threshold corrected by the slope $\tau_s+\Lambda \frac{\partial f}{\partial x}$, if
 shear stress   $\tau$ is too small,
 there is no source if the wall .
 
note that  we have a $\frac{\partial}{\partial x }q_s$ term,
 so that in a homogenous case, or steady case
 $q_s=(\tau - \Lambda \frac{\partial f}{\partial x}- \tau_s)_+$
 
 
 
We second we update the topography (Exner relation)
 $$\frac{\partial z_{bump}}{\partial t }    = - \frac{\partial}{\partial x }q_s$$
 
 
 As the computation uses the hydrodynamic time scale, we write with *Qt*  the ratio of times
 $$\frac{1}{Q_t} \frac{\partial z_{bump}}{\partial t } = - \frac{\partial}{\partial x}q_s$$
 
where $Q_t=$ hydro time/ granular time, so $Q_t \ll 1$
*/
/**
 We update the bump $z_{bump}$, the total topography is $z_b=-\alpha x + z_{bump}$ */
event updatebump (i++) {
    alphabump=alphamax;
    double eps=.15,xc=0.05;
    eps=.2,xc=0.075;
    alphabump=1.;
    foreach(){
        if(t<1){zbump[] = 0.5*alphabump*eps*exp(-sq((x-xc)/eps/eps/eps*.5)) ;}
        else{
            zbump[] =  max(zbump[] + Qt*(dt)*(-qsx[]) , 0);
            if(zbump[]<1e-10)zbump[]=0;
        }
        zb[]= - alpha*x + zbump[];
    }
}
/**
 

*/

event fchangesol (i++) {
    if(t>1){
/**
 First compute the corrected stress due to the slope of the bump: it is favorable  in the lee side
 $$\tau - \Lambda \frac{\partial f}{\partial x}$$
 then compute the source of $\partial_x q_s  + V q_s$:
 $$Ar = V (\tau - \Lambda \frac{\partial f}{\partial x}- \tau_s)_+$$
   */
 
 foreach()
  {  
  taud[]=tau[] - lambda*(zbump[+1]-zbump[-1])/2/Delta;
  Ar[]=0;
  
 if(taud[]>taus){
          Ar[]=V*(taud[]-taus); }
     else{Ar[]=0;}
  };
        
/** valeur maximale admissible de l'arrachement:

 if $A_r$ too large, $z_{bump}$ may be negative, we dig too much. So we find the max value possible $A_{rm}$
 $$\frac{z_{bump}^{n+1}- z_{bump}^{n}}{Qt \Delta t}  = -\frac{\partial}{\partial x }q_s$$
but $$ - \frac{\partial}{\partial x }q_s =
 - A_r +  V q_s$$
 so the max is
 $$A_{rm} =\frac{  z_{bump}^{n}}{Qt \Delta t} +V q_s$$
 
 
         */
 
 foreach()
  {
   Arm[] =  zbump[]/dt/Qt +  V*(qs[ ]+qs[-1])/2;
   if(Ar[]>=Arm[]){ Ar[]=Arm[];}
  }
/** calcul du flux de materiau  $q_s$ by implicit integration of
 
   $$\frac{\partial}{\partial x }q_s  =  - V q_s + Ar $$

discretized as:
 $$\frac{ q_s(i) - q_s(i-1)}{\Delta x}=
 - V \frac{ q_s(i) + q_s(i-1)}{2} + \frac{ Ar_{i} + Ar_{i-1}}{2} $$
 
 we put $q_s=0$ at the entrance and up to $x=.03$ (ad hoc value)
 
 we recompute then the new $\frac{\partial}{\partial x }q_s$
 as $\frac{q_s(i) - q_s(i-1)}{\Delta x}$
 which is  used for the Exner relation
 
 this loop is right in a *for(i=0;i<N;i++)* context
*/  
 foreach()
  {qs[]= (qs[-1] *(1/Delta - 0.5 * V) +  (1*Ar[]+1*Ar[-1])/2)/(1./Delta + 0.5 * V );  
     // qs[] = qs[]*(x>.03);
      //qs[] = qs[]*(x<(L0-.03));
  }
  foreach()
   qsx[] = (qs[]-qs[-1])/Delta;
  
/**
estimate the velocity of the dune $c_0=q_s/z_b$, verify mass conservation.
 */
   masse=0;
   double zbm=0;
   foreach()
   { 
    if(zbump[]>zbm){
      zbm=zbump[]; 
    if((zbump[-1]>0)&&(zbump[1]>0)) c0=(qs[]/zbump[]+qs[]/zbump[])/2; }
    masse = masse + (zbump[]) * Delta; }
}
}
/**
## Code: save/ plot
Finaly we save fields,
*/
event sauv (t<=tmax;t+=.25) {
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
          ue.x[] = 1;  // bug
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
      fprintf(fq,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
        x, h[], tau[], ue.x[],t,Fr,nu,delta1[],zb[],zbump[],eta[],Qc[],qs[],qsx[]);}
  fprintf(fq," \n");
  fflush(fq);
  fclose(fq);

 FILE *fg =  fopen("dune.txt","a");  
 foreach() { 
      fprintf(fg,"%g %g %g %g %g %g\n",
        x, zbump[],t,c0,masse,V);}
  fprintf(fg," \n");
  fflush(fg);
  fclose(fg);
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

event gnuplot (t+=.1) {
   static FILE * fp = popen ("gnuplot 2>/dev/null", "w");
     if (i == 0)
        fprintf (fp, "set term X11\n");
   // fprintf (fp,
    fprintf (fp,
             "set title \"t=%lf \" \n p [0:][:]'-' u 1:2 w l lc -1,'' u 1:3  w l lc 3\n",t);
    foreach()
        fprintf (fp, "%g %g %g\n", x, zb[] ,eta[]);
    
    fprintf (fp, "e\n\n");
    fflush (fp);
}

/**
 
## Run
To compile and run  :

~~~bash
   qcc -O2 -o MLSWbump13dune MLSWbump13dune.c -lm; ./MLSWbump13dune

    ./MLSWbump13dune 1200


~~~
 
using `Makefile`

~~~bash
  make MLSWbump13dune.tst
  make MLSWbump13dune/plots
  make MLSWbump13dune.c.html 
  open MLSWbump13dune.c.html
~~~  

## Plots

bla
bla 



starting from a given bump, it is deformed by the flow, and chanegs the flow, and moves
 finally at constant velocity
 
 
~~~gnuplot
 reset
 p [:][0:1]'dune.txt'u 1:($3>0? $2:NaN) t'' w l
~~~

bla
bla
 
 
Plot of the dune $z_{bump}$, the shear stress on it
$\partial_z u$, and the flux of sediment  $q_s$ at $t_{max}$.
 
~~~gnuplot
p[.15:]'xproff.txt'u 1:10 t'dune ' w l,'' u 1:3 t 'shear' w l,'' u 1:13 t'qs' w l
~~~
 
bla
bla
~~~gnuplot
p 'log' u 3:15 w l
~~~
 
bla
bla
~~~gnuplot
p 'log' u 3:12 w l
~~~
 
 
## Links
 

 
## Bibliography

* [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/CISM/tuyau_CISM.pdf) Double Deck (pipe flows)

* [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
"Ecoulements en milieux naturels" Cours MSF12, M2 SU

* [K.K.J. Kouakou Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/kouakoulagree06.pdf) 
 Evolution of a model dune in a shear flow

* file:///Users/pyl/macintoshHD/ZZZZZZ/Aeq0/help.html

*/
