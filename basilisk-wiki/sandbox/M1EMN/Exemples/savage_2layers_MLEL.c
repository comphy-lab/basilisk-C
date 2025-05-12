/**
#    Collapse of a heap with "Two-layer fluid model" on inclined plane    

 
From the paper by Shah et al. "Two-layer fluid flows on inclined surfaces"

 ![animation  ](savage_2layers_MLEL/animate.gif)


We solve with MLEL [MultiLayer Euler Lagrange solver](hydro.h) in 1D :

 the 2 layers set :

$$\left\{\begin{array}{l}
         \partial_t h_0+\partial_x (hu)_0=0\\
         \partial_t h_1+\partial_x (hu)_1=0\\
  \partial_t (hu)_0+ \partial_x (hu^2)_0= - g \partial_x \eta   - \partial_x (h q)_0 + \left[q \partial_x z \right]_0 
   - \mu g (h_0+h_1) \\
  \partial_t (hu)_1+ \partial_x (hu^2)_1= - g  \partial_x \eta - \partial_x (h q)_1   + \left[q \partial_x z \right]_1   -c_f u_1^2
      \end{array}\right.$$
 




- where the friction is Coulomb friction in the granular part, and turbulent in the upper layer.  This is oversimplified....

- where   $q$ is the (hydrostatic)
pressure deviation $q(z)  = \int_0^z g \frac{\Delta \rho(T)}{\rho}dz$ due to (small) density variations which are described by a scalar field $T$, 
 which is
advected by the flow, and an associated "equation of state"
$\Delta\rho(T)$ :
$$
\begin{aligned}
\partial_t \left( h T \right)_0 + \partial_x \left(
h u T \right)_0 & = 0, \\
\partial_t \left( h T \right)_1 + \partial_x \left(
h u T \right)_1 & = 0
\end{aligned}
$$

- with at  t=0, a given double heap $h_0(|x|<0.5,t=0)=.5$ and $h_1(|x|<.5,t=0)=.5$ so that $\int_{-\infty}^{\infty} h_0(x,0)dx =1$.

- where 
$R=\rho_1/\rho_0$ is relative density of upper fluid
and 
$M=\mu_0/\mu_1$ is relative density of lower  fluid.



The self similar solution and the avalanche cases have been done with the diffusive wave approximation allready (see two first links below). 
The other links correspond to the one fluid case (second Huppert problem).
*/

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#define rho0 1.
double R ; 
#define drho(T) ((R-1.)*T)
#include "layered/dr.h"
// declare parameters
double M;
double tmax;
/**
 position of domain `X0`, length `L0`, no dimension `G=1`
run with not so points (less is not enough)
*/
int main() {
  X0 = -2;
  L0 = 15;
  N = 1024;
  nl = 2;
  G = 1.;
  M = 5 ;
  R = 1.3; 
  tmax=70; 
  DT = 1;
  //DT = (L0/N)*(L0/N)/500;
  run();
}
/** start by a rectangular column  of different densities, slope -0.1 
*/
event init (i = 0)
{
  foreach(){
   zb[] = - 0.1*x;
 
   h[0,0,0]= (0.5*(fabs(x-3)<1))/2;
   h[0,0,1]= (0.5*(fabs(x-3)<1))/2;   

   T[0,0,0]= 1;//(1*(fabs(x)<.5))/2;
   T[0,0,1]= 0;//((fabs(x)<.5))/2;
}   
}
u.x[left] =dirichlet(0);
u.x[right]=dirichlet(0);
eta[right]=neumann(0);
/** 
 
so we can write the split for friction 
$$
 \left\{\begin{array}{l}
  \partial_t (u)_0 &= 
    - \mu g (h_0+h_1)\\
  \partial_t (u)_1 & =
  -c_f \frac{u_1^2}{h_1}
      \end{array}\right. 
$$
so that a semi implicit discretization gives 
 $$\frac{  u_0^{n+1} - u_0^n }{\Delta t} =   ...  \text{  and  }
 \frac{  u_1^{n+1} - u_1^n }{\Delta t} =   -c_f \frac{u_1^{n+1}u_1^{n}}{h_1^n}$$
    inverted by hand to give $u_0^{n+1}$ and $u_1^{n+1}$ as function of the previous ones...
*/
 
event coulomb_friction (i++) {
   fprintf (stdout,"%g \n",t);
  double mu,U,c1;
  
  mu=.3;
  c1=0.05;
  foreach() {
    U=sqrt(u.x[0,0,0]*u.x[0,0,0]);
    if(U>0){
      u.x[0,0,0] = max(U -dt *( mu * G )*(1+h[0,0,1]/h[0,0,0]),0)*u.x[0,0,0]/U;}
    U=sqrt(u.x[0,0,1]*u.x[0,0,1]);  
    double a = h[0,0,1] < dry ? HUGE : 1. + c1*dt*U/h[0,0,1];
    u.x[0,0,1] /= a;
        }
      }
  
/** saving some profiles 
*/



event outputfile (t <= tmax;t+=2) {
  foreach()
  fprintf (stderr,"%g %g %g %g %g \n", x, h[0,0,0], h[0,0,1], t ,zb[]);
  fprintf (stderr,"\n");
}
/** animation 
*/
event animatedplot (t+=.1) {
    static FILE * fp = popen ("gnuplot -persist 2> /dev/null", "w");
    if(t==0) fprintf (fp,"set term gif animate;set output 'animate.gif';set size ratio .333333;set grid\n");
    fprintf (fp,"set title 'Glissement en 1D --- t= %.2lf '\n"
      "t= %.2lf  ; "
      "p[%g:%g][-1.5:1]  '-' u 1:($2+$3+$4) t'free surface' w l lt 3,"
      "'' u 1:($2+$4) t'lower layer' w l lt 4,\\\n"
      "'' u 1:4 t'topo' w l lt -1\\\n",
           t,t,X0,X0+L0);
    foreach()
    fprintf (fp,"%g %g %g %g %g\n", x, h[0,0,0], h[0,0,1], zb[], t);
    fprintf (fp,"e\n\n");
    fflush (fp);  
    if(t==tmax){
        fprintf (fp,"! cp animate.gif a2.gif \n");}

}

 
/**
end of subroutines

 
#  Run

To compile and run  using `make`
 
~~~bash
make savage_2layers_MLEL.tst
make savage_2layers_MLEL/plots
make savage_2layers_MLEL.c.html
source ../c2html.sh savage_2layers_MLEL
open savage_2layers_MLEL.c.html

~~~

#   Results  
 
 An example of collapse along a slope

~~~gnuplot collapse
 set xlabel 'x'
 set ylabel 'h(x,t)'
 set arrow from 0,.47 to 0,0.08
 set label "t" at 0.5,1
 set arrow from 3,.05 to 7,0.05
 set label "t" at 8,0.06
 p [-5:][:1] 'log'  u 1:($3+$2)  t'free surface' w l
~~~
 
the upper and lower layers 
 
~~~gnuplot collapse
 set xlabel 'x'
 set ylabel 'h(x,t)'
 set arrow from 0,.47 to 0,0.08
 set label "t" at 0.2,0.1
 set arrow from 3,.05 to 7,0.05
 set label "t" at 8,0.06
 p [-5:][:1] 'log' u 1:2 t'bot' w l,'' u 1:3 t 'top' w l,'' u 1:($2+$3)  t'free surface' w l
~~~
  
 
 
 
 
# Links

 * [http://basilisk.fr/sandbox/M1EMN/Exemples/savagestaron.c]()
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/damb_dressler.c]() 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/cohesive_savagehutter.c]()
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/front_poul_ed.c]()

 with Bingham
 * [Bingham 1D collapse on a incline](bingham_collapse_noSV.c)

 
# Bibliography
 * [Huppert H.](http://www.itg.cam.ac.uk/people/heh/Paper47.pdf)
 * Kasturi S. Shah, Samuel S. Pegler and Brent M. Minchew
Two-layer fluid flows on inclined surfaces  
J. Fluid Mech. (2021), vol. 917, A54, doi:10.1017/jfm.2021.273
 * Doyle   Hogg  . Mader b,  Sparks 
 A two-layer model for the evolution and propagation of dense and dilute regions of pyroclastic currents
 Journal of Volcanology and Geothermal Research 190 (2010) 365–378
 
 */

 
