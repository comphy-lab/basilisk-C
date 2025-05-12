/**

# Dupuit Boussinesq Equation (with and without Rain)

## Problem
 
Filling of a 1D aquifere due to a continuous rain,
there is an output of water in $x=0$ where water goes out from soil. 
 In the aquifere we have Darcy flow.
 
 
 
 
 
 
 Animation shows the surface of the aquifer (indigo line) in the porous media. The green line corresponds to the filling height due to rain if there is no output at left ($x=0$ is the output), at right a Neumann ;  
 
![a gif](/sandbox/M1EMN/HYDGEO/dupuitboussinesq/drough.gif)
 



## Equation 

Explicit resolution of 1D Dupuit Boussinesq Equation ($\phi$ porosity)
$$\phi \frac{\partial h}{\partial t} +  \frac{\partial Q(h)}{\partial x}=R, $$
where the flow rate $Q$ obeys Darcy's law ($K$  hydraulic conductivity):
$$Q= 
 - K h \Bigg(   \frac{\partial   h}{\partial   x}
  \Bigg)$$
domain is an aquifere for $0\le x \le \infty$, initialy, it is empty at first, $h(t=0=0)$ then, for  $t>0$ rain is falling, with 
  $R$ rate of rain.
   
   


## About discretisation 

Finite volumes 
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
scalar hp[],hpc[],nu[];
scalar dQ[]; 
double R,phi,K;
double n,dt,tmax; 
char s[80];

h[left]  = dirichlet(0);
h[right] = neumann(0);
 

int main() {
  L0 = 10;
  X0 = 0;
  R = .1;
  phi = 1.;
  K = 1;
  N = 512; 
  DT = (L0/N)*(L0/N)*.20;
  tmax = 5;

{ 
  sprintf (s, "x.txt");
  FILE * fp = fopen (s, "w"); 
  fclose(fp);
  sprintf (s, "shape.txt");
  FILE * fs = fopen (s, "w");  
  fclose(fs);
  run();
 }
}
event init (t = 0) {
  foreach(){
    h[] = 0 ;
    Q[] = 0;
    }
  boundary ({h,Q});
  }
 

event printdata (t += 0.05, t < tmax){
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

  foreach()
    nu[] =  K  *  ( h[0,0] + h[1,0] )/2  ;
  boundary ({nu});      

  foreach()
    Q[] =  - nu[]  *hp[] ;
    boundary ({Q});

  foreach()
    dQ[] =  ( Q[1,0] - Q[0,0] )/Delta;
  boundary ({dQ});

  foreach(){
    h[] +=   dt*(- dQ[] + R )/phi;
  }
  boundary ({h});  
}

/**

## Run
Then compile and run, either by hand either with make:

~~~bash
qcc  -g -O2 dupuitboussinesq.c -o dupuitboussinesq
./dupuitboussinesq > out

make dupuitboussinesq.tst;make dupuitboussinesq/plots    
make dupuitboussinesq.c.html ; open dupuitboussinesq.c.html
 
source ../c2html.sh dupuitboussinesq
~~~

## Results

 Playing here with `gnuplot` and `every` :
 

 
 `every I:J:K:L:M:N`

`I` Line increment/ `J` Data block increment/ `K` The first line/ `L` The first data block/ `M` The last line/  `N` The last data block

ex: 

~~~bash

every ::3 skip the first 3 lines
every ::3::5  plot from the 4-th to 6-th lines
every ::0::0  plot the first line only
every 2::::6  plot the 1,3,5,7-th lines
every :2  plot every 2 data block
every :::5::8 plot from 5-th to 8-th data blocks
 
 
 p'shape.txt' ev :20::0::100  // plot from  0 to 100th every 20
~~~


plot the first and last computed shapes

~~~gnuplot
  reset
  set key bottom
  set xlabel "x"
  set ylabel "h(x,0), R Tmax"
  stats 'shape.txt' u 1;
  k = (STATS_records/512)-1
  p[][0:]'shape.txt' every :::k::k t'last' w lp,'' every :::0::0 t 'init' w l ,0.5 t 'rain'

~~~



`gnuplot` generates a gif:
filling the aquifere (with the rain if no output):



~~~gnuplot 
reset
set xlabel "x"
set ylabel "h(x,t)    R t "
!echo "p[:][:.5]'shape.txt' ev :::k::k title sprintf(\"t=%.2f\",k*.05) w l,.1*k*.05 t 'rain with no output';k=k+1 ;if(k<100) reread " > mov.gnu
set term gif animate;
set output 'drough.gif'
k=0
l'mov.gnu'
~~~

Self similar solution : solution of the height of the aquifere as a function of time is proportional to $t$ and the self similar variable is $x/t$

~~~gnuplot
  reset
  set key bottom
  set xlabel "x/t"
  set ylabel "h(x,t)/t "
  stats 'shape.txt' u 1;
  k = (STATS_records/512)-1
  p[:1][0:]'shape.txt' u ($1/($4)):(($2)/$4) ,'' every :::0::0 t 'init' w l  

~~~


 
# Links
 
 
 * see [http://basilisk.fr/sandbox/M1EMN/HYDGEO/dupuit2D.c]()
 
 * see [http://basilisk.fr/sandbox/M1EMN/HYDGEO/dupuitboussinesq.c]()
 
# Bibliography

* Adrien Guerin [PhD](http://adrienguerin.fr/papers/manuscrit.pdf)
Dynamique de l’écoulement dans un aquifère non confiné, 2015

* [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv_aquifere.pdf)
 "Equations de Saint Venant et application, Ecoulements en milieux naturels,
écoulements en milieux souterrains" Cours MSF12, M1 Sorbonne Université


*/

