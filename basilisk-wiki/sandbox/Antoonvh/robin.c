/**
![Victor Gustave Robin could have looked like Peter Gustav Lejeune Dirichlet (or Claude Louis Marie Henri Navier)](https://alchetron.com/cdn/peter-gustav-lejeune-dirichlet-4db41bd3-b41c-4a34-b3b9-80d9a762eff-resize-750.jpeg)

# Mixed boundary conditions

From [Yong Hui's page](../yonghui/smalltest/robin.c) we know that
for a boundary condition of the form,

$$ a \mathtt{s} + b \frac{\partial{\mathtt{s}}}{\partial \mathbf{n}} = c, 
$$
The discretized version for $2b\neq -a\Delta$ reads,
$$
 s [ \mathtt{ghost} ] = \left[ \frac{2 c \Delta}{2b + a\Delta} \right]+ \Bigg \langle \frac{2b - a\Delta }{2b + a\Delta} \Bigg \rangle s[\quad] 
$$

The ansatz is that this may be expressed as a *linear mix* of a
`dirichlet` and a `Neumann` condition.

$$s[\mathtt{ghost}] = A\cdot \mathtt{dirichlet}(D) + B\cdot\mathtt{neumann}(N).$$

Expanding the macros gives:

$$s[\mathtt{ghost}] = 2AD - As[\quad]+ B\Delta N + Bs[
\quad]$$
$$ = \left[2AD + B\Delta N\right]+  \langle B-A\rangle s[\quad].$$

Each term at the rhs forms a seperate equation. Such that we have four
unkowns ($A,B,D,N$) and only two equations. Such underdetermined
system may have many solutions. Lets check if a solution exist for $A
= 1$ and $N = 0$. The equation for the square-bracketed term then
becomes:

$$2D =  \frac{2 c \Delta}{2b + a\Delta},\rightarrow$$
$$ D = \frac{c \Delta}{2b + a\Delta}.$$

Next, for the second, angle-bracketed-term equation,

$$B - 1 = \frac{2b - a\Delta }{2b + a\Delta} \rightarrow$$
$$B  = \frac{2b - a\Delta }{2b + a\Delta}  + 1.$$

Wrapping it up (in many braces):
*/
#define robin(a,b,c) ((dirichlet ((c)*Delta/(2*(b) + (a)*Delta))) + ((neumann (0))*((2*(b) - (a)*Delta)/(2*(b) + (a)*Delta) + 1.)))
/**
Mind the `devide-by-zero` risk.

## A test

[Vyaas
Gururajan](https://groups.google.com/g/basilisk-fr/c/-YSthke6VBo)
coded a neat test case. Special care must be taken, because it is not
always critical.

~~~gnuplot Looks OK
set xlabel 'x'
set ylabel 'u'
set size square
set grid
plot 'out' u 1:3 w l lw 6 t 'Analytical Solution', '' w l lw 3 t 'Approximate Solution'
~~~
*/
#include "grid/multigrid1D.h"
#include "poisson.h"
/** 
Here is the analaytical solution taken from [Professor
 Fitzpatrick's notes](http://farside.ph.utexas.edu/teaching/329/lectures/node66.html):
*/
double alphaL = 2., betaL  = -1., gammaL = -3.;
double alphaR = 1., betaR  =  3., gammaR = -2.;

double analyticalU (double x) {
  double d = alphaL*alphaR + alphaL*betaR - betaL*alphaR;
  double g = (gammaL*(alphaR + betaR) - betaL*(gammaR - (alphaR+betaR)/3.))/d;
  double h = (alphaL*(gammaR - (alphaR + betaR)/3.) - gammaL*alphaR)/d;
  return (g + h*x + sq(x)/2. - sq(sq(x))/6.);
}

scalar u[], rhs[];

u[left]  = robin (alphaL, -betaL, gammaL); //Mind the sign!
u[right] = robin (alphaR,  betaR, gammaR);

int main() {
  init_grid (1 << 7);
  foreach()
    rhs[] = 1. - 2.*sq(x);
  poisson (u, rhs, tolerance = 1e-9);
  foreach() 
    printf("%g %g %g\n",x, u[], analyticalU(x));
}
