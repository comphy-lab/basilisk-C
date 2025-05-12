/**
#    Collapse of a heap with "Two-layer fluid model" on inclined plane    

 
From the paper by Shah et al. "Two-layer fluid flows on inclined surfaces"

 ![animation  ](viscous_collapse_2layers_MLEL/animate.gif)


We solve with [Multilayer Euler Lagrange solver](hydro.h) in 1D,
the 2 layers set :
$$\left\{\begin{array}{l}
         \partial_t h_0+\partial_x (hu)_0=0\\
         \partial_t h_1+\partial_x (hu)_1=0\\
         \partial_t (hu)_0+ \partial_x (hu^2)_0= - g \partial_x \eta   - \partial_x (h q)_0 + \left[q \partial_x z \right]_0 
   +\mu_0\frac{6 (6 h_0 u_0 + 2 h_1 M u_0 - 3 h_0 u_1)}{
   h_0^2 (3 h_0 + 4 h_1 M)}\\
  \partial_t (hu)_1+ \partial_x (hu^2)_1= - g  \partial_x \eta - \partial_x (h q)_1   + \left[q \partial_x z \right]_1
  + \mu_1\frac{(18 u_0 - 12 u_1)}{(3 h_0 h_1 + 4 h_1^2 M)}      \end{array}\right.$$
 




- where the friction is laminar (say $-c_{00} u_0-c_{01} u_1,-c_{10} u_0 -c_{11} u_1$). The coefficients $c_{ij}$ are such that under a pressure gradient 
 the equilibrium is 
$$(hu)_0=-\frac{h_0^2 (h_0 )}{3 \mu_0}(\frac{\partial p_0}{\partial x})-\frac{{h_0}^2 (h_1 )}{2 \mu_0}(\frac{\partial p_1}{\partial x})
\text{  and  }
(hu)_1=-\frac{h_0^2 h_1 }{2 \mu_0}(\frac{\partial p_0}{\partial x})-\frac{{h_1}^2   (3 {h_0}
   {\mu_1}+{h_1} {\mu_0})}{3 {\mu_0} {\mu_1}} \frac{\partial p_1}{\partial x}  \text{ TBC: R?} $$
(as the flow is very viscous, the convective term is finally small, 
and there  is only is a balance between pressure gradient and viscosity so that the flux is the previous lubrication equations ones)

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
#include "hydroNN.h"
#define rho0 1.
double R ; 
#define drho(T) ((R-1.)*T)
//#include "layered/dr.h"
// declare parameters
double M;
double tmax;
double nu_eq(double shear,double pipi){
    double nu_eq;
    nu_eq = 1;
    return nu_eq;}
/**
 position of domain `X0`, length `L0`, no dimension `G=1`
run with 2048*8 points (less is not enough)
*/
int main() {
  X0 = -1;
  L0 = 15;
  N =  2048*8;
  nl = 2;
  G = 1.;
  M = 5 ;
  R = 1.; 
  tmax=200;//400;// 1000; 
  DT=HUGE;
  run();
}
/** start by a rectangular column  of different densities 
*/
event init (i = 0)
{
  foreach(){
   zb[] = - x;
 
   h[0,0,0]= (0.5*(fabs(x)<.5))/2;
   h[0,0,1]= (0.5*(fabs(x)<.5))/2;  // plus haut plante??

 //  T[0,0,0]= 1;//(1*(fabs(x)<.5))/2;
 //  T[0,0,1]= 0;//((fabs(x)<.5))/2;
}   
}
/** 
Remember the equilibrium is 
$$(hu)_0=-\frac{h_0^2 (h_0 )}{3 \mu_0}(\frac{\partial p_0}{\partial x})-\frac{{h_0}^2 (h_1 )}{2 \mu_0}(\frac{\partial p_1}{\partial x})
\text{  and  }
(hu)_1=-\frac{h_0^2 h_1 }{2 \mu_0}(\frac{\partial p_0}{\partial x})-\frac{{h_1}^2   (3 {h_0}
   {\mu_1}+{h_1} {\mu_0})}{3 {\mu_0} {\mu_1}} \frac{\partial p_1}{\partial x}$$

so we can write the split for friction 
$$
 \left\{\begin{array}{l}
  \partial_t (u)_0+= 
   +\mu_0\frac{6 (6 h_0 u_0 + 2 h_1 M u_0 - 3 h_0 u_1)}{
   h_0^2 (3 h_0 + 4 h_1 M)} \text{ TBC}\\
  \partial_t (u)_1=
  + \mu_1\frac{(18 u_0 - 12 u_1)}{(3 h_0 h_1 + 4 h_1^2 M)}   \text{ TBC}
      \end{array}\right. 
$$
as 

 $$\frac{\partial u_0}{\partial t} = -C_{00} u_0 -C_{01} u_1 $$
 $$\frac{\partial u_1}{\partial t} = -C_{10} u_0 -C_{11} u_1 $$
  so that a semi implicit discretization gives 
 $$\frac{  u_0^{n+1} - u_0^n }{\Delta t} =    -C_{00} u_0^{n+1} -C_{01} u_1^{n+1} \text{  and  }
 \frac{  u_1^{n+1} - u_1^n }{\Delta t} =   -C_{10} u_0^{n+1} -C_{11} u_1^{n+1}$$
 where the $C_{ij}$ depend on $M$ $R$ and the $h_0^n$ and $h_1^n$.  
 The matrix is inverted by hand to give $u_0^{n+1}$ and $u_1^{n+1}$ as function of the previous ones...
*/
event friction (i++) {
  foreach()
  {  
     double U0 = u.x[0,0,0];
     double U1 = u.x[0,0,1];
     double h0 = h[0,0,0];
     double h1 = h[0,0,1];
// friction (implicit scheme)
   u.x[0,0,0] = h0*(h0*h1*(3*h0 + 4*h1*M)*U0 + 6*dt*(2*h0*U0 + 3*h1*U1))*
   pow(36*pow(dt,2) + h1*(3*h0 + 4*h1*M)*pow(h0,2) + 12*dt*(3*h0*h1 + pow(h0,2) + M*pow(h1,2)),-1);
   u.x[0,0,1] = (h1*(3*h0 + 4*h1*M)*U1*pow(h0,2) + 6*dt*(6*h0*h1*U1 + 3*U0*pow(h0,2) + 2*M*U1*pow(h1,2)))*
   pow(36*pow(dt,2) + h1*(3*h0 + 4*h1*M)*pow(h0,2) + 12*dt*(3*h0*h1 + pow(h0,2) + M*pow(h1,2)),-1);
    }
}
/** saving some profiles 
*/
event outputfile (t <= tmax;t+=100) {
  foreach()
  fprintf (stderr,"%g %g %g %g %g \n", x, h[0,0,0], h[0,0,1], t ,zb[]);
  fprintf (stderr,"\n");
}
/** animation 
*/
event animatedplot (t+=10) {
    static FILE * fp = popen ("gnuplot -persist 2> /dev/null", "w");
    if(t==0) fprintf (fp,"set term gif animate;set output 'animate.gif';set size ratio .333333\n");
    fprintf (fp,"\nset grid\n");
    fprintf (fp,"set title 'Glissement en 1D --- t= %.2lf '\n"
      "t= %.2lf  ; "
      "p[%g:%g][-1:1]  '-' u 1:($2+$3+$4) t'free surface' w l lt 3,"
      "'' u 1:($2+$4) t'lower layer' w l lt 4,\\\n"
      "'' u 1:4 t'topo/L0' w l lt -1\\\n",
           t,t,X0,X0+L0);
    foreach()
    fprintf (fp,"%g %g %g %g %g\n", x, h[0,0,0], h[0,0,1], zb[]/L0, t);
    fprintf (fp,"e\n\n");
    fflush (fp);  
    if(t==tmax){
        fprintf (fp,"! cp animate.gif a2.gif \n");}

}

#if 0
// inspired by 
// http://basilisk.fr/src/test/overflow.c
// do more beautifull animations
void setup (FILE * fp)
{
  fprintf (fp,
     "set pm3d map\n"
     "# jet colormap\n"
     "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25"
     " 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625"
     " 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
     "unset key\n"
     //    "set cbrange [10:20]\n"
     "set xlabel 'x '\n"
     "set ylabel 'depth (m)'\n"
     "set xrange [-2:20]\n"
     "set yrange [0:1]\n"
     );
}

void plot (FILE * fp)
{
  fprintf (fp,
     "set title 't = %.2f '\n"
     "sp '-' u 1:2:3\n",t);
  foreach (serial) {
    double z = 0*zb[];
    fprintf (fp, "%g %g %g %g\n", x, z, u.x[], T[]);
    foreach_layer() {
      z += h[];
      fprintf (fp, "%g %g %g %g\n", x, z, u.x[], T[]);
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  //  fprintf (fp, "pause 1\n");
  fflush (fp);  
}

event gnuplot (t += 10)
{
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  if (i == 0)
    setup (fp);
  if (getenv ("DISPLAY")) {
    fprintf (fp, "set term x11\n");
    plot (fp);
  }
  fprintf (fp,
     "set term pngcairo font \",10\" size 800,500\n"  
     "set xrange [-2:20]\n"
     "set output 'plot-%04d.png'\n", i);
  plot (fp);
}

event figures (t += 10)
{
  FILE * fp = popen ("gnuplot 2> /dev/null", "w");  
  setup (fp);
  fprintf (fp,
     "set term pngcairo font \",10\" size 800,500\n"
     "set xrange [-1:20]\n"
     "set output 'T-%g.png'\n", t);
  plot (fp);
}

event moviemaker (t = end)
{
  system ("for f in T*.png; do convert $f ppm:- && rm -f $f; done | "
    "ppm2mp4 movie.mp4");
}

#endif
/**
end of subroutines

 
#  Run

To compile and run  using `make`
 
~~~bash
make viscous_collapse_2layers_MLEL.tst
make viscous_collapse_2layers_MLEL/plots
make vviscous_collapse_2layers_MLEL.c.html
source ../c2html.sh viscous_collapse_2layers_MLEL
open viscous_collapse_2layers_MLEL.c.html

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
  
 
~~~gnuplot collapse
set key left
set xlabel "\eta= x/t^{1/3}"
set ylabel "t^{1/3}h_u, t^{1/3}h_l" 
tt=20
M=5
R=1.2
p[-0.5:1.5][] 'log' u ($1/$4**(1./3)):($4>=tt?(($3)*$4**(1./3)) : NaN) t'Comp._{upper}'w l,\
'' u ($1/$4**(1./3)):($4>=tt?($2*$4**(1./3)) : NaN) t'Comp._{lower}'w l,\
'' u ($1/$4**(1./3)):($4>=tt?(($2+$3)*$4**(1./3)) : NaN) t'Comp._{sum}'w l,\
 sqrt(x/R) t'sqrt' w l,sqrt(x/M*R) t'sq....'' w l
 ~~~
 
 
# Links
 * [http://basilisk.fr/sandbox/M1EMN/TEST/shah_flux.c]()
 * [http://basilisk.fr/sandbox/M1EMN/TEST/shah.c]()
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse.c]()
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_noSV.c]()
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_ML.c]()

 with Bingham
 * [Bingham 1D collapse on a incline](bingham_collapse_noSV.c)

 
# Bibliography
 * [Huppert H.](http://www.itg.cam.ac.uk/people/heh/Paper47.pdf)
 * Kasturi S. Shah, Samuel S. Pegler and Brent M. Minchew
Two-layer fluid flows on inclined surfaces  
J. Fluid Mech. (2021), vol. 917, A54, doi:10.1017/jfm.2021.273
 * E.E. Doyle  , A.J. Hogg  H.M. Mader, R.S.J. Sparks
 A two-layer model for the evolution and propagation of dense and dilute regions of pyroclastic currents
 Journal of Volcanology and Geothermal Research 190 (2010) 365–378

 */

 