/**

# Collapse of a rectangular heap on an incline with Bagnold velocity profile

 In hydrodynamics, it is usual to simplify momentum equation in neglecting inertia. Hence, it is usual to consider the friction/ slope/ pressure gradient equilibrium and deduce the flow rate  $Q$ as a function of $h$, $-Z'_b$ and
 $\partial_xh$.
 Exemple:  Huppert sploutch  :
 $Q = -g\frac{\partial h}{\partial x} h^3/(3\nu)$, flood wave, diffusive wave,
 and Liu and Mei/ Balmforth (Bingham flow) etc. 
 
 
 Then one puts this in mass conservation and has only ONE equation:
 $$\frac{\partial h}{\partial t} +  \frac{\partial Q}{\partial x}=0$$
 
 Here we consider that the flow follows the linearised $\mu(I)$ rheology and is always a
 Bagnold profile (in Huppert's case it is always a half Poiseuille etc).
 
 
Linearisation at small $I_\alpha$,around $\tan \alpha-\mu_s$ small with  $\tan \alpha=-Z'_b>0$
 we obtain the flow rate: 
 $$ Q = \frac{2 }{5} ( \frac{-Z_b' - \mu_s -\partial_x h}{\Delta \mu})I_0 \sqrt{g h} \frac{h^2}{d}$$
  that we substitute in
 $$\frac{\partial h}{\partial t} +  \frac{\partial Q}{\partial x}=0$$
 we solve it juste like we do for Huppert viscous collapse
 [http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_noSV.c]() 
and 
 [http://basilisk.fr/sandbox/M1EMN/Exemples/viscolsqrt_noSV.c]()
 Huppert's full sploutch (first and second problems) or Balmforth Bingham collapse
 [http://basilisk.fr/sandbox/M1EMN/Exemples/bingham_collapse_noSV.c]() cases. 
 
 
Note that as $Z_b'<0$ and  $\mu_s>0$ we must have $-Z_b' - \mu_s >0$ to have a flow. 


# Gif animated of the slump :
 
An animation of the result:
 
![a gif animated](bagnold_ondesimple_noSV/slump.gif)
 

We use as previously the Finite Volume method with $ad$ $hoc$ estimation of the flux

~~~gnuplot
set size ratio .3
set samples 9 
set label "h i-1" at 1.5,3.1
set label "h i" at 2.5,3.15
set label "h i+1" at 3.5,2.5
set xtics ("i-2" 0.5, "i-1" 1.5, "i" 2.5,"i+1" 3.5,"i+2" 4.5,"i+3" 5.5)
set arrow from 2,1 to 2.5,1
set arrow from 3,1 to 3.5,1
set label "Q i" at 2.1,1.25
set label "Q i+1" at 3.1,1.25

set label "x i-1/2" at 1.5,0.25
set label "x i" at 2.4,0.25
set label "x i+1/2" at 3.,0.25

set label "x"  at 0.5,2+sin(0) 
set label "x"  at 1.5,2+sin(1)
set label "x"  at 2.5,2+sin(2) 
set label "x"  at 3.5,2+sin(3) 
set label "x"  at 4.5,2+sin(4)  
f(x) = (2+sin(x))*(x<=4)
p[-1:7][0:4] f(x) w steps not,f(x) w impulse not linec 1
~~~

## Code
*/

#include "grid/cartesian1D.h"
#include "run.h"

scalar h[];
scalar Q[];
scalar hp[],hpc[],nu[],Zb[];
scalar dQ[]; 
double n,dt,S,A,tmax,inct=0.0001;
char s[80];

Q[left] = dirichlet(0);
h[left] = neumann(0);
hp[left] = neumann(0);
Q[right] = neumann(0);
h[right] = neumann(0);

int main() {
  L0 = 10;
  X0 = 0;
  S = 1.;
  N = 512;
  DT = .00001;
  tmax = 8; // 8

{ 
  sprintf (s, "x.txt");
  FILE * fp = fopen (s, "w"); 
  fclose(fp);
  sprintf (s, "shape.txt");
  FILE * fs = fopen (s, "w");  
  fclose(fs);
  inct = DT;
  run();
 }
}
/** 
initial heap, a  rectangle 

if we change this in `h[] =1*(fabs(x)<1.5) `we generate a front just  like 
in 
[http://basilisk.fr/sandbox/M1EMN/Exemples/front_poul_ed.c]()
as the Pouliquen's front is the friction/ slope/ pressure gradient  equilibrium (negligible inertia), and mass conservation.


*/
event init (t = 0) {
  foreach(){
    h[] =1*(fabs(x-2)<1.5)  ;
    Q[]=0;
    }
  boundary ({h,Q});
  }
/**
front tracking*/
event output (t += .1; t < tmax) {
  double  xf=0,xe=0;
    inct = inct*2;
  inct = min(1,inct);
  foreach(){
   xf = h[] > 1e-4 ?  max(xf,x) :  xf ;
   xe = h[] > 1e-4 ?  min(xe,x) :  xe ;
  }

  sprintf (s, "front.txt");
  FILE * f = fopen (s, "w");  
  foreach()
    fprintf (f, "%g %g %g \n", fmin((x-xf),0), h[], xe-xf);   
  fclose(f);
  fprintf (stderr, "%g %g %g \n", t, xf, xe);
  sprintf (s, "x.txt");
  FILE * fp = fopen (s, "a");   
   fprintf (fp, "%g %g  \n", t, xf); 
  fclose(fp);
}

event printdata (t += 0.25, t < tmax){
  sprintf (s, "shape.txt");
  FILE * fp = fopen (s, "a");  
  foreach()
    fprintf (fp, "%g %g %g %g \n", x, h[], Q[], t);
  fprintf (fp, "\n");
  fclose(fp);
}

event integration (i++) {
  double dt = DT;

  dt = dtnext (dt);

  foreach()
    hp[] =  ( h[0,0] - h[-1,0] )/Delta;
  boundary ({hp});

  foreach()
    hpc[] =  ( h[1,0] - h[-1,0] )/2/Delta;
  boundary ({hp});

/**
 $-Z_b' - \mu_s -\partial_x h>0$ to have a flow. 
*/  
  foreach(){
      S = 1;
      //S = (x < 4 ? S : 0.00);
      A = S - hp[];
      nu[] = (A >0.0 ? A : 0.00); 
    }
  boundary ({nu});      
/**
 compute
 $$ Q = \frac{2 }{5} ( \frac{S  -\partial_x h}{\Delta \mu}) I_0\sqrt{g h} \frac{h^2}{d}$$
 with say $S=-Z_b' - \mu_s$
 and  by choice of units $\frac{2 }{5} ( \frac{I_0}{\Delta \mu}) \sqrt{g } \frac{1}{d}=1$
 so that we code 
 $Q =   ( S  -\partial_x h)_+  h^{1/2+2}$
 
*/
    foreach()
      Q[] =  nu[] * pow((h[]+ h[-1,0])/2,5./2);
    boundary ({Q});
/**
 that we substitute in
 $$\frac{\partial h}{\partial t} +  \frac{\partial Q}{\partial x}=0$$
 
 */
  foreach()
    dQ[] =  ( Q[1,0] - Q[0,0] )/Delta;
  boundary ({dQ});

  foreach()
    h[] +=  - dt*dQ[];
  boundary ({h});  
}

/**

## Run
Then compile and run, either by hand either with make:

~~~bash
qcc  -g -O2 bagnold_ondesimple.c -o bagnold_ondesimple
./bagnold_ondesimple > out

source ../c2html.sh bagnold_ondesimple
~~~

## Results

 
~~~gnuplot
 reset
 set xlabel "x"
 set ylabel "h(x,t)"
 p[][0:1.5]'shape.txt'  w l
 
~~~

 
 
 Playing here with `gnuplot` and `every` :
 

plot the first and last shapes :

~~~gnuplot
  reset
  set xlabel "x"
  set ylabel "h(x,0),h(x,infinity)"
  stats 'shape.txt' u 1;
  k = (STATS_records/512)-1
  p[][0:1.5]'shape.txt' every :::k::k t'last' w lp,'' every :::0::0 t 'init' w l

~~~


Front as function of time

~~~gnuplot 
  reset
  set xlabel "t"
  set ylabel "x"
  set key bottom
 p[0.001:]'x.txt'  t'cul' w l
~~~

`gnuplot` generates a gif:

~~~gnuplot 
reset
set xlabel "x"
!echo "p[:][0:1]'shape.txt' ev :::k::k title sprintf(\"t=%.2f\",k*.05) w l;k=k+1 ;if(k<30) reread " > mov.gnu
set term gif animate;
set output 'slump.gif'
k=0
l'mov.gnu'
~~~

 

 
# Links

* [http://basilisk.fr/sandbox/M1EMN/Exemples/floodwave.c]() cases.
* [Bingham 1D collapse on a incline](bingham_collapse_noSV.c)
* [Herschel-Bulkley 1D collapse on a incline](herschel-column-noSV.c)
* [http://basilisk.fr/sandbox/M1EMN/Exemples/bingham_collapse_noSV.c]()
* [Bingham RNSP collapse on a incline](bingham_collapse_ML.c)
* [Viscous collapse of a incline](http://basilisk.fr/sandbox/M1EMN/Exemples/viscolsqrt_noSV.c)


# Bibliographie 

* Neil J. Balmforth, Richard V. Craster, Alison C. Rust, Roberto Sassi 
["Viscoplastic flow over an inclined surface"](http://www.math.ubc.ca/~njb/Research/revslump.pdf),
J. Non-Newtonian Fluid Mech. 139 (2006) 103–127
* K.F. Liu, C.C. Mei, 
"Slow spreading of Bingham fluid on an inclined plane", J. Fluid Mech. 207 (1989) 505–529.
* [http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/mainM2EMN.pdf]()
* Borzsonyi et al 2005 Two Scenarios for Avalanche Dynamics in Inclined Granular Layers
* Borzsonyi et al 2008
Avalanche dynamics on a rough inclined plane 

*/


