/**

<script type="text/javascript"
    src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>

# A transport problem using the wave propagation algorithm

This example demonstrates how to use the wave propagation solver in 
[wpa_transport.h]() for the transport equation, given by 

$$
  \partial_t \mathbf q + \nabla(\mathbf u \cdot \mathbf q) = 0
$$

where $\mathbf u = (u(x,y),v(x,y))$ is a prescribed velocity field and 
$\mathbf q(x,y,t)$ is a vector of tracer quantities. 

The solver is based on the wave propagation algorithms, first described by 
R. J. LeVeque [LeVeque, 2002](#le:2002). 

To test the accuracy of the solver, we use a velocity field that returns an initial field to its 
original position after a fixed time and compute the difference between the initial
and final times. 

For comparison, we also run the Basilisk solver for the same problem. 

This example was adapted from the advection test [/src/test/advection.c]().
*/

#include "grid/quadtree.h"

/**
To use the wave propagation algorithm, set the macro ${\tt USE\_WAVEPROP}$ to 1. 
*/


#define USE_WAVEPROP 1

/**
Include the wave propagation algorithm (WPA) solver, or the Basilisk solver.
*/

#if USE_WAVEPROP
#include "wpa_transport.h"
#else
#include "advection.h"
#endif


/**
Set tracer field to one or more tracers to be transported by the velocity field.
*/ 
scalar f[];
scalar *tracers = {f};

int main() {
    origin (-0.5, -0.5);   

#if !USE_WAVEPROP    
    gradient = minmod;
#endif    

/**
Run four simulations to see how the error behaves. 
*/

  CFL = 0.8;
  for (N = 64; N <= 512; N *= 2)
    run();
}


/**
## Initial conditions

The wave propagation algorithm requires an estimate of the initial time step. In this problem, subsequent time steps will be adjusted using the CFL parameter, defined above.
*/

#define bump(x,y) (exp(-100.*(sq(x + 0.2) + sq(y + .236338))))

event init (i = 0) {
    foreach() {
        f[] = bump(x,y);    
    }
    boundary((scalar*) tracers);
}

/**
## Velocity field
*/


/**
For the wave propagation algorithm, we define a cell centered velocity field.  
*/

#if USE_WAVEPROP
event velocity(i++) {
    trash ({u,v});
    foreach() {
        u[] =  1.5*pi*sin(2.*pi*(t)/5.)*sin(pi*(x + 0.5))*cos(pi*(y + 0.5))/pi;  
        v[] = -1.5*pi*sin(2.*pi*(t)/5.)*cos(pi*(x + 0.5))*sin(pi*(y + 0.5))/pi;  
    }
    boundary ((scalar *){u, v});    
/**
We define an initial time step, estimated based on the velocity field. 
*/
    dt_initial = 2e-2;  
}

#else

/**
For the the Bell-Colella-Glaz solver used in Basilisk, the velocity field is 
defined at cell faces.  In this example, we use a streamfunction $\Psi(x,y)$ and
define velocities as $u = -\Psi_y$ and $v = \Psi_x$. 
*/
event velocity (i++) {
    vertex scalar psi[];
    foreach_vertex() {
        psi[] = - 1.5*sin(2.*pi*t/5.)*sin((x + 0.5)*pi)*sin((y + 0.5)*pi)/pi;        
    }
    trash ({u});
    struct { double x, y; } f = {-1.,1.};
    foreach_face() {
        u.x[] = f.x*(psi[0,1] - psi[])/Delta;        
    }
    boundary ((scalar *){u});
}
#endif

/**
## Logging information at initial and final times

We can use this to test the numerical conservation of the scheme.
*/

event logfile (t = {0,5})  {
    stats s = statsf (f);
    fprintf (stderr, "# %8.4f %24.16e %24.16e %24.16e\n", t, s.sum, s.min, s.max);
}

/**

## Error computation

At $T=5$, the initial tracer field should return to its initial position.  We then compute
the difference between the initial and final times.  Using the ${\tt norm}$ operator, we
can then compute the average (1-norm), root-mean-square (2-norm) and maximum ($\infty$-norm)
of the error field. 

*/

event field (t = 5) {
    scalar e[];
    foreach() {
        e[] = f[] - bump(x,y);       
    }
    norm n = normf (e);
    fprintf (stderr, "%-5d %24.16e %24.16e %24.16e\n", N, n.avg, n.rms, n.max);
}

/**

## Output

*/

event movie (t += 0.1) {
    if (N == 256) {
        output_ppm (f, linear = true, file = "f.mp4", n = N, min = 0, max = 1);
    }
}

/** 

![Animation of tracer in velocity field](transport/f.mp4)(loop)

*/


/**

## References

~~~bib
@book{le:2002,
    Author = {Randall J. LeVeque},
    Date-Added = {2013-12-26 20:56:55 +0000},
    Date-Modified = {2013-12-27 02:43:06 +0000},
    Publisher = {Cambridge University Press},
    DOI = {10.1017/CBO9780511791253.004},
    Title = {Finite volume methods for hyperbolic problems},
    Year = 2002}
~~~
*/
