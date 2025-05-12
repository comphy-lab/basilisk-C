/** 
# Solution of the potential flow though an orifice

## the problem 

solution of the "Basic problem" in 2D, potential incompressible flow 
$$u =\frac{\partial p}{\partial x}, \;\; v = \frac{\partial p}{\partial y},\;\;\; 
  \frac{\partial u}{\partial x} +  \frac{\partial v}{\partial y}=0$$
  to be solved in the upper half space $y>0$ filled with a porous media so that velocity is 0 at infinity (but pressure is not zero to have a flow). At the bottom, there is a slit
 ($-1<x<1$, $y=0$) where pressure is constant: 
 $p=0$.   And there is no penetration for ($|x|>1$, $y=0$), so normal velocity is  $\frac{\partial p}{\partial y}=0$ 


Solution is proposed in Sneddon and in Lamb.
It is  $x=cosh(p) cos(\psi)$ and  $y=sinh(p) sin (\psi)$

iso lines of pressure are the ellipses
$$ 
\frac{x^2}{cosh^2 p} + \frac{y^2}{sinh^2 p } =1
$$

## Code
*/
#include "run.h"
#include "poisson.h"
#define MAXLEVEL 8

scalar p[], source[];
face vector beta[];
mgstats mgp;

/**

far away, 
$cosh(p)= sinh(p)= e^p/2$, hence iso pressure are $p \simeq  log(2 \sqrt{x^2+y^2})$

Of course, this corresponds to the expected solution of a source far enough.


as $arcsinh(\sqrt{x^2+y^2}) \simeq log(2 \sqrt{x^2+y^2})$ we write the approximate BC on the top right and left (and "mix" for fun). 

On the bottom, the mixed condition.

*/ 
p[right] = dirichlet(log(sqrt(4.*(x*x+y*y)))) ;  
p[left]  = dirichlet(asinh(sqrt(x*x+y*y))) ;   
p[top]   = dirichlet(log(sqrt(4.*(x*x+y*y)))) ;  
p[bottom] = fabs(x)<= 1 ? dirichlet(0): neumann(0);
/**
 domain is large
*/
int main()
{    
    L0=50.;
    Y0=0;
    X0=-L0/2.;
    init_grid (1 << MAXLEVEL);
    run();
}
/**
 coefficient  of porosity is constant 
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
    p[] = source[] = 0.;
  boundary ({p});
}
/**
At every timestep, but after all the other events for this timestep
have been processed (the '`last`' keyword), we update the pressure
field $p$ by solving the Poisson equation with  coefficient
$\beta$. */

event pressure (i++, last)
{
/** 
solve $\nabla \cdot (\beta \nabla p  )= s$ 
with [http://basilisk.fr/src/poisson.h](http://basilisk.fr/src/poisson.h)
*/
  mgp = poisson (p, source, beta);
}
/**
error
*/
event logfile (i++)
{
    stats s = statsf (p);
    fprintf (stderr, "%d %g %d %g %g %g\n",
             i, t, mgp.i, s.sum, s.min, s.max);
}
/**
Save in a file
*/
event sauve (i++,last)
{
FILE *  fpc = fopen("pressure.txt", "w");
output_field ({p}, fpc, linear = true);
fclose(fpc);
fprintf(stdout," end\n");
}
/**
## Run
 Then compile and run:

~~~bash
 qcc -O2 -Wall -o darcyLambSneddon darcyLambSneddon.c -lm
./darcyLambSneddon
~~~

or better

~~~bash
 make darcyLambSneddon.tst;make darcyLambSneddon/plots; make darcyLambSneddon.c.html
~~~


## Results
Results


~~~gnuplot   pressure
set pm3d map
set palette rgbformulae 22,13,-31;
unset colorbox
set xlabel "x  iso p"
 set size (.5*1.6),1
splot [][:] 'pressure.txt' u 1:2:3   not
reset
~~~
 
Figure to compare with Lamb 1932 Art 66 p 73 showing the ellipses of iso pressure, the stream-lines are 
  
~~~gnuplot   pressure
L0=50.
reset
set view map
set size (.5*1.6),1
unset key
unset surface
 set contour base
 set cntrparam levels incremental 0,.1,log(2*L0)
 splot [][:] 'pressure.txt' u 1:2:3  w l not
~~~ 



Plot along $x=0$ showing that analytical solution $p(0,y) = arcsinh(y)$ is close to the log far away (source solution).

~~~gnuplot   pressure
set key bottom
  set logscale x
 plot [:][:] 'pressure.txt' u 2:(abs($1)<.01?($3):NaN)  t'num.',asinh(x),acosh(x),log(2*x)
~~~ 






## bibliography

* Sneddon I.N.  Mixed boundary value problems in potential theory 1966, Wiley

* Lamb Hydrodynamics 1932 

* see [http://basilisk.fr/src/hele-shaw.h](http://basilisk.fr/src/hele-shaw.h)

* see [./darcysilo.c]() to see what happens if we put walls


*/
