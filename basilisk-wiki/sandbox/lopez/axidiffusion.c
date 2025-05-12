/**
# Axisymmetric diffusion 

In a cylindrical domain that goes from radius $r_1$ to $r_2$ are
imposed the following initial and boundary conditions for a concentration $c(r,t)$,
$$ 
c(0,r) = a + b \log \, r + c \, J_o (r) \, , c(t,r_1) = \alpha  \, , c(t,r_2) = \beta \, .
$$ 
$J_o$ is the Bessel function of zeroth order and $a$, $b$ and $c$ are constant given by
$$
a = \frac{\beta \log r_1 - \alpha \log r_2}{\log(r_1/r_2)} \quad b = \frac{\alpha-\beta}{\log(r_1/r_2)} 
$$
This distribution diffuses radially  according to,
$$ 
\frac{\partial c}{\partial t} = \frac{1}{r}
\frac{\partial}{\partial r}\left( r \frac{\partial c}{\partial r} \right).  
$$ 

If $r_1$ and $r_2$ are zeros of the Bessel function, $J_o(r_1)=
J_o(r_2)= 0$, the time evolution of the concentration is
$$ 
c(t,r) = a + b \log \, r + c \, e^{-t}  J_o (r) \,.
$$
*/

#include "axi.h"
#include "run.h"
#include "diffusion.h"

scalar c[];

#define alpha 1
#define beta 0.
#define r1 2.40483
#define r2 5.52008
#define A ((beta*log(r1) - alpha*log(r2))/log(r1/r2))
#define B ((alpha - beta)/log(r1/r2))
#define C 1.

c[bottom]   = dirichlet(alpha);
c[top]  = dirichlet(beta);

/** 
The width of the square computational box is $(r_2 -r_1)$. The origin
is at the center of the bottom boundary of the domain,which is located
at $r_1$.*/

int main() {
  L0 = r2-r1;
  origin(-L0/2.,r1);
  DT = HUGE [0];
  N = 256;
  init_grid (N);
  run();
}

event init (i = 0) {
  foreach() 
    c[] = A + B*log(y) + C *j0(y);

  boundary ((scalar *){c});
}

/** 
The maximal time step is set to $\Delta t = 0.04$. */

event integration (i++) {
  dt = dtnext (0.04);

 /**
  In the present implementation of *diffusion.h*, $\theta$ has to be
  redone in each step. It is a simple cost given the efficiency of the 
  present implementation. Also, since the diffusivity is 1 in
  this example, the field D in *diffusion* is set to the metric
  factor, *fm*. */
  
  scalar cmv[];
  foreach()
    cmv[] = cm[];
  boundary({cmv});
  
  diffusion (c, dt, fm, theta = cmv);  
}

event output (t += 0.5 ; t <= 2) {
  int no = 30;
  static FILE * fp = fopen("cprof", "w");
  for (int j = 0; j  <= no;  j++) {
    double y = r1 + (r2-r1)/(no + 0.00001)*j;
    fprintf (fp, "%g %g %g %g\n", t, y,
	     interpolate (c, 0., y), A + B*log(y) + C*j0(y) *exp(-t));
  }
  fprintf(fp,"\n");
}

/**
## Results

~~~gnuplot Concentration profiles for several instants.
plot 'cprof' u 2:3  t 'Basilisk','cprof' u 2:4 w l t 'Analytical'    
~~~
*/
