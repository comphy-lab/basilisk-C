/** 
# Solution of the potential flow  
 

## the problem :
...

solution of the "Basic problem" in 2D, potential incompressible flow 
$$u =\frac{\partial \psi}{\partial y}, \;\; v = -\frac{\partial \psi}{\partial y},\;\;\; 
  \frac{\partial u}{\partial x} +  \frac{\partial v}{\partial y}=0$$
  to be solved in the upper half space $1>y>0$ and $-0.5 < x < 0.5$ (the walls are in $x =\pm 1/2$ $y>0$, and $\frac{\partial p}{\partial x}=0$ on the walls)
  filled with a porous media,  pressure is given at the top.
   At the bottom, there is a very small  slit
 ($- \varepsilon < x < \varepsilon$, $y=0$) where pressure is constant: 
 $p=0$. We take $\varepsilon=.1$.  And there is no penetration for ($|x|>.1$, $y=0$), so normal velocity is  $\frac{\partial p}{\partial y}=0$ 
 
 
 

# Code
*/
#include "run.h"
#include "poisson.h"
#include "embed.h"
#define MAXLEVEL 8

scalar psi[], source[];
double eps=.1;  // OK Up to 0.025 for level 8
face vector beta[];
mgstats mgp;

/**

## boundary conditions    
 far away  

...

*/ 
psi[right] = neumann(0); 
psi[left]  = dirichlet(y);   
psi[top]   = dirichlet(y) ;  
psi[bottom] = dirichlet(y) ;  

/**
 domain is unit
*/
int main()
{    
    L0=1.;
    Y0=0;
    X0=-L0/2.;
    init_grid (1 << MAXLEVEL);
    run();
}
/**
 coefficient  is constant 
 */
event init (i = 0) {
    foreach_face() {
        beta.x[] = 1; 
    }  
}
/**
no source
*/
event defaults (i = 0)
{ 
  foreach()
    psi[] = source[] = 0.;
  boundary ({psi});
}
/**
At every timestep, but after all the other events for this timestep
have been processed (the '`last`' keyword), we update the pressure
field $p$ by solving the Poisson equation with  coefficient
$\beta$. */

event streamfunction (i++, last)
{
/** 
solve $$\nabla \cdot (\beta \nabla \psi  )= s$$ 
with [http://basilisk.fr/src/poisson.h](http://basilisk.fr/src/poisson.h)
*/
  mgp = poisson (psi, source, beta);
}
/**
error
*/
event logfile (i++)
{
    stats s = statsf (psi);
    fprintf (stderr, "%d %g %d %g %g %g\n",
             i, t, mgp.i, s.sum, s.min, s.max);
}
/**
Save in a file
*/
event sauve (i++,last)
{
FILE *  fpc = fopen("stream.txt", "w");
output_field ({psi}, fpc, linear = true);
fclose(fpc);
fprintf(stdout," end\n");
}
/**
# Results
## Run
To compile and run:

~~~bash
 qcc -O2 -Wall -o ailepot ailepot.c -lm
./ailepot
~~~

or more clean

~~~bash
 make ailepot.tst; make ailepot/plots ; make ailepot.c.html
~~~


## Plots

Results, figures of iso pressure

~~~gnuplot
set pm3d map
set palette rgbformulae 22,13,-31;
unset colorbox
set xlabel "x  iso \psi"
 set size (.5*1.6),1
splot [][:] 'stream.txt' u 1:2:3   not
reset
~~~

Figure of iso pressure,  we see the transition from Lamb solution to constant pressure gradient solution via the 
infinite sum of sinks from Paterson's solution.
  
~~~gnuplot
L0=1.
reset
set view map
set size (.5*1.6),1
unset key
unset surface
 set contour base
 set cntrparam levels incremental -1.8,.1,2.5
 splot [][:] 'stream.txt' u 1:2:3  w l not
~~~ 
 
Along $x=0$ we have  



~~~gnuplot
set key bottom
set xlabel "y"
set ylabel "psi(0,y)"

 plot [:1][:] 'stream.txt' u 2:(abs($1)<.01?($3):NaN)  t'num.',x t'y'
~~~ 

 
 # Links
 
 
 http://basilisk.fr/src/test/neumann.c
 
 
 

# Bibliography

*  
*  
*/
