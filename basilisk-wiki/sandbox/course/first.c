/**
# Approximation of the first derivative of a function (1st version)

Given a function $f(x)=\cos(2\pi x)$, compute an approximation of its
first derivative.

We will discretise the function at regular intervals on the segment 
$x\in [0:1]$. */

#include "grid/cartesian1D.h"

/**
The *main()* function is the starting point of any C program. */
int main() {
  /** We declare a new (integer) variable which holds the number of grid points. */  
  int n = 10;
  /** We create the grid. */
  init_grid (n);
  /** We allocate a new scalar field *f*... */
  scalar f[];
  /** initialise it... */
  double k = 2.*pi;
  foreach()
    f[] = cos(k*x);
  /** Our discretisation now looks like this
  
  ~~~gnuplot
  set term pngcairo size 600,200
  set xtics 0,0.1,10
  set grid xtics
  f(x)=cos(2.*pi*x)
  unset key
  plot [-0.05:1.05] f(x), "< seq 0.05 0.1 0.95" u 1:(f($1)) pt 7, \
                          "< seq -0.05 1.1 1.05" u 1:(f($1)) pt 7
  ~~~
  
  where the green values (at the center of the discretisation cell) 
  were set by the *foreach()* loop and the blue "ghost" values were
  set by the call to *boundary()*. By default, the boundary conditions
  are "symmetry" i.e. $f'(0) = f'(1) = 0$. 
  
  We will store the first derivative in field *df*. */
  
  scalar df[];
  
  /** To compute the first derivative, we will use a local *stencil*, for example

  ~~~gnuplot
  set label 'f[]' at 0.25,f(0.25)+0.2 center
  set label 'f[1,0]' at 0.37,f(0.35)+0.2 center
  set label 'f[-1,0]' at 0.13,f(0.15)-0.2 center
  plot [-0.05:1.05] f(x), "< seq 0.05 0.1 0.95" u 1:(f($1)) pt 7, \
                          "< seq -0.05 1.1 1.05" u 1:(f($1)) pt 7
  ~~~
  
  We can then write an approximation of the first derivative for each cell as */

  foreach()
    df[] = (f[] - f[-1,0])/Delta;
  
  /** where $\Delta$ is the grid spacing (1/10 for the moment). Of course this means
  that the derivative is defined at the mid-point between `f[]` and `f[-1,0]` (0.2 
  in the example above). Given that the *foreach()* loop only iterates over the green
  points, this also means that the derivative is not defined for $x=1$. 
  
  Note that we could also have defined `df[]` using `f[]` and `f[1,0]`, it is just a 
  matter of choice.
  
  To check whether this works, we can write the values to *standard output* which is
  redirected to a file called 'out'. */
  
  foreach()
    printf ("%g %g\n", x - Delta/2., df[]);
  
  /** Note that the *x* coordinate in the *foreach()* loop is the position of the 
  center of the cell (the green dots), so that we need to substract $\Delta/2$ to
  get the position of the derivative. */
  
}

/** After the program has finished we can plot both the points in 'out' and the 
exact value of the derivative $-2\pi\sin(2\pi x)$.

~~~gnuplot
unset label
plot [0:1]-2.*pi*sin(2*pi*x), 'out' pt 7 lt 4
~~~ 

Although it looks like the error is small, we need to measure it and do a proper
convergence study. This will be the [next exercise](first2.c). */
