/** 
# Solution of the potential flow though a small orifice in the bottom of a silo (by Matched Asymptotic Expansion)
 

## the problem :

solution of the "Basic problem" in 2D, potential incompressible flow 
$$u =\frac{\partial p}{\partial x}, \;\; v = \frac{\partial p}{\partial y},\;\;\; 
  \frac{\partial u}{\partial x} +  \frac{\partial v}{\partial y}=0$$
  to be solved in the upper half space $1>y>0$ and $-0.5 < x < 0.5$ (the walls are in $x =\pm 1/2$ $y>0$, and $\frac{\partial p}{\partial x}=0$ on the walls)
  filled with a porous media,  pressure is given at the top.
   At the bottom, there is a very small  slit
 ($- \varepsilon < x < \varepsilon$, $y=0$) where pressure is constant: 
 $p=0$. We take $\varepsilon=.1$.  And there is no penetration for ($|x|>.1$, $y=0$), so normal velocity is  $\frac{\partial p}{\partial y}=0$ 


 There is here a small parameter $\varepsilon$, we do a Matched Asymptotic Expansion using this small parameter. 

## solution in complex plane for asymptotic approximations of the problem

### Inner Problem: Lamb
First we look at small $\varepsilon$ near the slit itself, so we rescale $x=\varepsilon \tilde x, y=\varepsilon \tilde y$, 
pressure stays at the same scale. The problem is now to 
solve the laplacian of pressure in the upper half plane $\tilde y >0$, with the slit
 ($-c < \tilde {x} < c$, (here $c=1$) $\tilde y=0$) where pressure is constant, say 
 $\tilde p=P_0$.   And there is no penetration on the wall for ($| \tilde x | > c$, $\tilde y=0$), 
 so normal velocity is  $\frac{\partial \tilde p}{\partial \tilde y}=0$ 



Solution for the case of a slit with the upper half plane $\tilde >0$, full of fluid is proposed in Sneddon and in Lamb 
($c$ is the width of the slit, which is one half).
It is  $\tilde x=c cosh(\tilde p) cos (\tilde \psi)$ and  $\tilde y=c sinh(\tilde p) sin (\tilde \psi)$

with  $\tilde z/c  = cosh(F(\tilde z))$ with $\tilde z = \tilde x +  i \tilde y$ and  $F=\tilde p + i \tilde psi$ a complex potential.
iso pressure are ellipses
$(\tilde x/c/cosh(\tilde p))^2 +   (\tilde y/c/sinh(\tilde p))^2=1$.
 Far enough $sinh(\tilde p) \simeq cosh(\tilde p) \simeq e^{\tilde p}/2$, so that 
$$\tilde p = P_0 + log (\tilde x^2 +   \tilde y^2) + log (2)$$


Note that the $P_0$ is an extra constant that we add for convinioance and matching, the strict Lamb solution is with $P_0=0$ 
(this solution has been tested here [darcyLambSneddon.c]() with $P_0=0$) 
 far enough it is a source. 

   
### Outer problem: Paterson

We turn now to the external problem, for large $x$ and $y$ measured in $\varepsilon$.

Hence we guess that we have to look at the problem of a source/sink $Log(x+i y)$  between two walls in $x=\pm 1$. 
 
A king of method of images is used: we superpose an infinite number of sources/sinks  ($a =1$ distance between one and the other). 
Solution of the problem (for a very small slit $c \ll a$) is
$$\Sigma_{-\infty}^{\infty}Log(z - n a)$$
 
The same problem exists for an infinite sum of  vortices: $(i \Sigma Log(z - n a))$. It was solved by Karman (as cited by Péres 1936), and 
is in Paterson (1984). Thanks to   Weierstrass   Borel factorization:
$$ sin ( \pi z)/(\pi  z) = \Pi_1^{\infty}   (1 -z^2/n^2)$$ 
we take this identity 
$$ sin ( \pi z) = \pi  z  (1 -z^2)  (1 -z^2/4 )...$$
develop it 
$$   sin (\pi z) = \pi  z  (1 -z)(1+z)   (1 -z/2 ) (1+z/2)  ... $$
hence 
$$  Log  sin ( \pi z)  =  Log ( \pi ) +   Log ( z) +  Log ( 1 -z) +  Log ( 1+z  ) +  Log ( 1+z/2 ) +    Log ( 1- z/2 ) +...$$
we change the order to recognize our sum
  $$Log  sin ( \pi z)  =  Log ( \pi ) +   Log ( z) +  Log ( 1 -z)  +  Log ( 1+z ) +    Log ( 2- z ) + Log(1/2) +  
    Log ( 2+ z ) + Log(1/2) +...$$
  which is
  $$Log  sin ( \pi z)  =  Log ( \pi ) +   \Sigma_{n}  Log(z - n )  +  \Sigma_{n>0} 2 Log( 1/n)$$
    
There is a  "constant" which is fact infinite  $(Log(\pi)+2\Sigma_{n>0} Log( 1/n))$,  as Paterson says "except a constant which does not affect velocity"
    
Final expression for complex potential 
 $$Log(sin(\pi z))= \Sigma_{n=-\infty}^{n=\infty} Log(z - n )+( 2 \Sigma_{n=1}^{\infty}   Log( 1/n)+ Log( \pi))$$
 

## Matching 


 for small $x$ and $y$ close to the sink, the outer problem:
  $$Log(sin(\pi z))= Log \pi + Log z + ...$$
this gives the outer pressure (a sink):
$$p= Log \pi + Log r + ...$$
but from Lamb inner problem, far from the slit, 
far away, 
$cosh( \tilde p) \simeq sinh( \tilde  p) \simeq  e^{\tilde p}/2$, 
hence iso pressure are $\tilde p \simeq  P_0+ log(2 \sqrt{\tilde x^2+\tilde  y^2})$ 
so $p \simeq P_0 + Log (2r/(c\varepsilon ))$.
Of course, this corresponds to the expected solution of a sink far enough.

Now we write the asymptotic matching 
$$lim_{\tilde y \rightarrow \infty} \tilde p =lim_{y \rightarrow 0} p$$ 
 
this gives :
$P_0 + Log(2) -Log( \varepsilon)= Log(\pi)$ which provides  the value of pressure $P_0$ that we have to put at the small slit


## Velocity

Velocity is constant and is $\pi$ at the top of the silo. 
Remember  $\tilde x=c cosh(\tilde p) cos (\tilde \psi)$ and  $\tilde y=c sinh(\tilde p) sin (\tilde \psi)$ with $c=1$.

At the slit, where $\psi$ goes from $0$ to $\pi$ and $\tilde y=0$ and pressure is 0, we have 
$d \tilde x= sin \psi d \psi$ and velocity is computed with $\frac{\partial \tilde y}{\partial \tilde p}$ so that
$\frac{\partial \tilde p}{\partial \tilde y} = \frac{1}{sin \psi}$.
So that the integral is 
$$ D = \int_{\tilde x =-1}^{\tilde x =1} \frac{\partial \tilde p}{\partial \tilde y} dx =
\pi.$$


Flow rate in this device is $\pi$.

# Code
*/
#include "run.h"
#include "poisson.h"
#define MAXLEVEL 8

scalar p[], source[];
double eps=.1;  // OK Up to 0.025 for level 8
face vector beta[];
mgstats mgp;

/**

## boundary conditions    
 far away $Log ( sin ( \pi z) ) \simeq  Log (exp( i  \pi z)/2) = -\pi  (i x - y) -Log[2]$
 so that the pressure is indeed linear $(\pi  y  -  Log 2 )$.
 We put at the top the  exact value $2.44658$ indeed note that $\pi - Log[2.]= 2.448$



close to the sink we have the overlapping between the two descriptions which is the source so that as $\varepsilon=.1$

$Log[.1] + Log[\pi] - Log[2]= -1.851$  gives the value of pressure that we put at the small slit

*/ 
p[right] = neumann(0); 
p[left]  = neumann(0);   
p[top]   = dirichlet(2.4465800) ;  
p[bottom] = fabs(x)<= eps ? dirichlet(log(eps) + log(pi) - log(2)): neumann(0);

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
solve $$\nabla \cdot (\beta \nabla p  )= s$$ 
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
# Results
## Run
To compile and run:

~~~bash
 qcc -O2 -Wall -o darcysilo darcysilo.c -lm
./darcyLambSneddon
~~~

or more clean

~~~bash
 make darcysilo.tst; make darcysilo/plots ; make darcysilo.c.html
~~~


## Plots

Results, figures of iso pressure

~~~gnuplot
set pm3d map
set palette rgbformulae 22,13,-31;
unset colorbox
set xlabel "x  iso p"
 set size (.5*1.6),1
splot [][:] 'pressure.txt' u 1:2:3   not
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
 splot [][:] 'pressure.txt' u 1:2:3  w l not
~~~ 
 
Along $x=0$ we have for the outer solution $Re[Log[i sinh(\pi y))$ for the complex potential, so that it is for pressure $log(sinh(\pi y))$.
The linear solution is  $\pi y-log(2)$ (in black), approximation for the outer solution.
We plot (in red)  the numerical solution along $x=0$  compared with analytical outer solution $log(sinh(\pi y))$ (in green) for a sink, they diverge near the slit as expected. 
 Lamb inner solution (in blue) rescaled by $\varepsilon$   gives the solution near the slit.
There is an overlap region for $y=O(\varepsilon)$ where the three are superposed 
(as expected by the Matched Asymptotic Expansion theory).



~~~gnuplot
set key bottom
set xlabel "y"
set ylabel "p(0,y)"

# plot [:][:] 'pressure.txt' u 2:(abs($1)<.01?($3):NaN)  t'num.',log(sinh(pi*x)) 
 plot [:1][:] 'pressure.txt' u 2:(abs($1)<.01?($3):NaN)  t'num.',log(sinh(pi*x)),\
 '../darcyLambSneddon/pressure.txt' u ($2*.1):(abs($1)<.01?($3-1.851):NaN)  t'num. Lamb' w l,pi*(x)-log(2) t'lin. approx' w l lc -1
~~~ 

Composite expansion:
$$p_{composite} = p_{inner}+p_{outer} - p_{common\;limit}$$

along the center line $x=0$, we have
$p_{inner}=\tilde p= log(\pi)-log(2)+log(\varepsilon)+asinh(y/\varepsilon)$, and $p_{outer}=log(sh(\pi y))$,
 and the common behavior that we remove is 
the sink $log(y)+log(\pi)$ so that we plot $p_{composite}$ compared to the numerics:

~~~gnuplot
set key bottom
set xlabel "y"
set ylabel "p(0,y)"
eps=.1
plot [:1][:] 'pressure.txt' u 2:(abs($1)<.01?($3):NaN) t'num.'\
,log(pi)-log(2)+log(eps)+asinh(x/eps)+log(sinh(pi*x))-log(x)-log(pi) t'composite'
~~~ 


# Bibliography

* R. Paterson - A First Course in Fluid Dynamics (1984, Cambridge University Press) page 410 

* J. Peres - Cours de Mécanique des fluides, (1936 Gauthier-Villars) page 181 

* see [https://fr.wikipedia.org/wiki/Théorème_de_factorisation_de_Weierstrass](https://fr.wikipedia.org/wiki/Théorème_de_factorisation_de_Weierstrass)

* see [http://basilisk.fr/src/hele-shaw.h](http://basilisk.fr/src/hele-shaw.h)

* see [darcyLambSneddon.c]()   

* Sneddon I.N.  Mixed boundary value problems in potential theory 1966, Wiley

* Lamb Hydrodynamics 1932 

* PY Lagrée  MAE [http://www.lmm.jussieu.fr/~lagree/COURS/M2MHP/MAE.pdf]()


*/
