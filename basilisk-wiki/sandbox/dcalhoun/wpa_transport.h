/**

<script type="text/javascript"
    src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>

# A solver for the transport problem

This solver solves the conservative transport problem given by
$$
  \partial_t \mathbf q + \nabla(\mathbf u \cdot \mathbf q) = 0
$$
where $\mathbf u = (u(x,y),v(x,y))$ is a prescribed velocity field and 
$\mathbf q(x,y,t)$ is a vector of tracer quantities. 

The solver is based on the wave propagation algorithms, first described by 
R. J. LeVeque [LeVeque, 1997](#le:1997). For a more general reference on 
hyperbolic problems, see [LeVeque 2002](le:2002).
*/

#include "waveprop.h"

/**

## Variable velocity field and tracer components

This example requires a cell-centered velocity field, represented here by Basilisk 
scalars $u$ and $v$.  The solver also expects a list of tracers, defined in the 
example file [transport.c]().  
*/

scalar u[];
scalar v[];

extern scalar* tracers;

/**

## Wave propagation parameters

This transport solver assumes that one or more tracers are all transported with the same velocity
field, defined above.  As a result, there is only a single wave propagating at each cell interface.
Here, we hardwire the number of waves for this solver.  For each wave, we define a limiter, also set by the user (see the example [transport.c]()).  
*/

#define MWAVES 1
int limiters[MWAVES];


event defaults(i=0) {

/**
The wave propagation algorithm in [waveprop.h]() requires several parameters which are specific to the problem being solved.  In this Riemann solver, we are solving a conservation law, and so set
*/

    conservation_law = true;

/**
To obtain a second order accurate scheme for the transport problem, we use the F-wave approach, described in [Bale, et al.](#ba-le-mi-ro:2003).  
*/    
    use_fwaves = true;

/**
Set the number of waves used in this solver. 

*/

    mwaves = MWAVES;

/**

We use the option to set a single limiter for all waves, since in this case, we only 
have a single wave. 
*/
    
    wave_limiter = 6;

/**
Although we only have a single wave in this example, 
we still set the pointer ${\tt mthlim}$ to a valid memory location.  Generally, the user 
may wish to override the default behavior set above and specify limiters for each wave. 
More importantly,  the wave 
propagation algorithm expects that ${\tt mthlim}$ is set and points to an integer array of length
${\tt mwaves}$.
*/

    mthlim = limiters;

/**

Set up the  list of state variables and auxiliary variables expected by the wave propagation
algorithm.

*/

    statevars = list_copy(tracers);
    auxvars = list_concat({u},{v});
}

/**

## Riemann solver for the transport equation


The Riemann solver for the transport problem decomposes the flux difference at each 
interface into eigenvalues and eigenvectors of flux Jacobian. For 
this problem, the flux Jacobian in the x direction, for example, is 
given by 
$$
\mathbf f'(q) = u,
$$ 
which is a scalar, so the eigen-decomposition is trivial.  

This Riemann solver is called for each cell face.  The arguments are the left and 
right states (stored in ${\tt ql}$ and ${\tt qr}$) and left and right auxiliary data 
(stored in ${\tt auxl}$ and ${\tt auxr}$). An additional parameter ${\tt dir}$ is used to 
indication the direction of the Riemann solve (x-direction or y-direction).  

The Riemann solver for the wave-propagation algorithm should then return waves, wave speeds and "fluctuations".  These will be used to compute second order corrections (if required) and update the solution.
*/

void wpa_rpn2(int dir, int meqn, int mwaves, 
              double *ql, double *qr, 
              double *auxl, double *auxr, 
              double *waves, double *speeds, 
              double *amdq, double *apdq, double* flux)  {

    if (mwaves != 1) {
        fprintf(stdout,"Error  (wpa_transport.h : mwaves != 1\n");
        exit(0);
    }

/**
The velocity field in the left and right cells are stored in input auxiliary arrays.
*/

    double ur = auxr[dir];  
    double ul = auxl[dir];


/**
The speed at the interface is taken as the average speed between the two adjacent cells.
*/
    double uhat = (ur + ul)/2.;

/**

### Waves

At each interface, the wave (defined here for propagation in the x-direction) is given by
$$
\mathcal W = u_r q_r - u_l q_l
$$
where $q_r$, $q_l$ are the advected quantities in two neighboring cells, and $u_r$ and $u_l$ are the horizontal velocities stored in the auxiliary array.

For the transport problem, the components of the wave vector correspond to the number of 
tracers we are advecting.
*/

    for(int m = 0; m < meqn; m++) {

        waves[m] = ur*qr[m] - ul*ql[m];    

/**

### Fluctuations

Fluctuations $\mathcal A^+ \Delta Q$ and $\mathcal A^- \Delta Q$ 
are defined according to the speed of the f-wave at each interface.  
*/                
        if (uhat >= 0) {
            amdq[m] = 0;
            apdq[m] = waves[m];    
        }
        else {
            amdq[m] = waves[m];
            apdq[m] = 0;    
        }        

/**
To ensure conservation, we also need a flux function defined in the Riemann solver.  The details of how conservation is imposed in the wave propagation algorithm can be found in 
[Berger and LeVeque](#be-le:1998).
*/

        flux[m] = ur*qr[m]; 
    }

/**

### Speeds

The wave speed at the interface is set and will be used by the wave propagation algorithm to compute second order corrections.
*/ 

    speeds[0] = uhat;
}

/**

## References

~~~bib
@article{le:1997,
    Author = {LeVeque, Randall J.},
    Date-Added = {2013-12-28 22:42:46 +0000},
    Date-Modified = {2016-04-28 20:51:11 +0000},
    Doi = {10.1006/jcph.1996.5603},
    Journal = {J. Comput. Phys.},
    Month = {Mar},
    Number = {2},
    Pages = {327-353},
    Publisher = {Elsevier -- Academic Press},
    Title = {Wave Propagation Algorithms for Multidimensional Hyperbolic Systems},
    Volume = {131},
    Year = {1997},
    DOI = {http://dx.doi.org/10.1006/jcph.1996.5603}}


@book{le:2002,
    Author = {Randall J. LeVeque},
    Date-Added = {2013-12-26 20:56:55 +0000},
    Date-Modified = {2013-12-27 02:43:06 +0000},
    Publisher = {Cambridge University Press},
    DOI = {10.1017/CBO9780511791253.004},
    Title = {Finite volume methods for hyperbolic problems},
    Year = 2002}

@article{ba-le-mi-ro:2003,
    Author = {Bale, Derek and LeVeque, Randall J. and Mitran, Sorin and Rossmanith, James},
    Date-Added = {2014-07-20 19:52:47 +0000},
    Date-Modified = {2014-07-20 19:53:57 +0000},
    Doi = {10.1137/S106482750139738X},
    Eprint = {http://dx.doi.org/10.1137/S106482750139738X},
    Journal = {SIAM J. Sci. Comput.},
    Number = {3},
    Pages = {955-978},
    Title = {{A Wave Propagation Method for Conservation Laws and Balance Laws with Spatially Varying Flux Functions}},
    Url = {http://dx.doi.org/10.1137/S106482750139738X},
    Volume = {24},
    Year = {2003}}

@article{be-le:1998,
    Author = {Berger, Marsha J. and Randall J. LeVeque},
    Date-Added = {2013-12-26 20:56:55 +0000},
    Date-Modified = {2018-08-02 15:42:11 +0100},
    Doi = {10.1137/S0036142997315974},
    Journal = {SIAM J. Num. Anal.},
    Number = 6,
    Pages = {2298--2316},
    Title = {Adaptive mesh refinement using wave-propagation algorithms for hyperbolic systems},
    Volume = 35,
    Year = 1998}
~~~
*/
