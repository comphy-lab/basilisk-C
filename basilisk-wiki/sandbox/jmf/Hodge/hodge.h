/** 
# Hodge Decomposition
Let $\Omega^k$ the space of the k-forms on a closed compact n-dimensional manifold, any k-form $\omega$ can be expressed as
$$ \omega = d \alpha + \delta \beta + \gamma$$
for some k-1 form $\alpha$, a k+1 form $\beta$ and $\gamma \in \Omega^k$.
Where in the exterior calculus (see [Exterior Algebra](https://en.wikipedia.org/wiki/Exterior_algebra))

- d is the differential exterior derivative
- $\delta$ is the codifferential defined by $\star d \star$

and $\star$ is the Hodge operator. 

Given a k-form $\omega$ there are unique $\alpha, \beta,\gamma$ (Hodge theorem, see H. Flanders, Differential Forms with Applications to the Physical Sciences (Dover Books on Mathematics)) which can be found by solving

- $\delta \omega = \delta d \alpha$
- $d \omega = d \delta \beta$
- $\gamma = \omega - d \alpha - \delta \beta$

In $R^2$ a vector $X = (u,v)$ is associated to a 1-form $\omega = u dx + v dy$

we solve $\delta \omega = \delta d \alpha$

- $\delta \omega = \star d \star \omega$ $\rightarrow$ $\delta \omega_X = \star d ( -v dx + u dy)$ 
- $\delta \omega = \star (-v_y dy \wedge dx + u_x  dx \wedge dy)$
- $\delta \omega = u_x  + v_y$

As $\omega_X$ is a 1-form that implies $\alpha \in \Omega^0$ a function of real values, we found then 

- $d \alpha = \alpha_x dx + \alpha_y \ dy$
- $\delta d \alpha = \star d \star d \alpha = \alpha_{xx} + \alpha_{yy}$

The first equation is then
$$ \alpha_{xx} + \alpha_{yy} = u_x  + v_y $$
a Poisson equation for $\alpha$.


Structure for the solver : 

- omega, 
- gradb, 
- rotc, 
- h 

are mandatory vectors giving the Hodge decomposition. If  the scalar alpha and beta exist they will return the scalar fields.
*/

#include "poisson.h"

struct Hodge {
  // mandatory
  vector omega;
  vector gradb;
  vector rotc;
  vector h;
  // optional
  scalar alpha;
  scalar beta;
};

/** Call the Hodge decomposition
*/

trace
mgstats hodge (struct Hodge p)
{

 vector omega =  p.omega;
 vector gradb = p.gradb;
 vector rotc = p.rotc;
 vector h = p.h;

  // Local
 scalar b[];
 scalar alpha[];
 scalar beta[];

/** Check the existence of alpha et beta */
  
 if (p.alpha.i) {
   alpha = p.alpha;
 }

 if (p.beta.i) {
   beta = p.beta;
 }

/** paranoia  */

boundary  ((scalar *){omega});

/** Filling the rhs of the Poisson equation 
*/
  
  foreach()
    {
#if order == 2
      b[] = (omega.x[1,0] - omega.x[-1,0])/2./Delta + (omega.y[0,1] - omega.y[0,-1])/2./Delta;
#else
      b[] =   (omega.x[-2,0] - 8.*omega.x[-1,0] + 8.*omega.x[1,0]  -omega.x[2,0])/12./Delta;
      b[] +=  (omega.y[0,-2] - 8.*omega.y[0,-1] + 8.*omega.y[0,1]  -omega.y[0,2])/12./Delta;
    }
#endif


/** 
Solve $\nabla^2 \alpha = b$ using the Poisson solver
*/

 mgstats s = poisson (alpha, b,tolerance = 1e-6);


/** 
Now we apply the same procedure to $d \omega = d \delta \beta$


$d\omega = \left( -\frac{\partial\,u}{\partial y} + \frac{\partial\,v}{\partial x} \right) \mathrm{d} x\wedge \mathrm{d} y$

and

$d \star d \star \beta = \left( \frac{\partial^2\,\beta}{\partial x ^ 2} + \frac{\partial^2\,\beta}{\partial y ^ 2} \right)d x\wedge d y$

we have still a Poisson equation

$$\nabla^2 \beta =  b = \left( -\frac{\partial\,u}{\partial y} + \frac{\partial\,v}{\partial x} \right)$$
*/


  foreach()
    {
#if order == 2
      b[] = (omega.y[1,0] - omega.y[-1,0])/2./Delta - (omega.x[0,1] - omega.x[0,-1])/2./Delta;
#else
      b[] =  (omega.y[-2,0] - 8.*omega.y[-1,0] + 8.*omega.y[1,0]  -omega.y[2,0])/12./Delta;
      b[] -= (omega.x[0,-2] - 8.*omega.x[0,-1] + 8.*omega.x[0,1]  -omega.x[0,2])/12./Delta;
    }
#endif
  boundary ({beta,b});


  /**
Solve Poisson problem
  */

 s = poisson (beta, b,tolerance = 1e-6);


/** 
Finally from the scalar fields $\alpha$ and $\beta$$ we compute the vectors 
$$gradb = d \alpha$$ 
$$ rotc = \delta \beta$$
and $\gamma$ by subtraction.
*/


 /** 
Compute vectors
  */
  

foreach()
    foreach_dimension()
#if order == 2
    gradb.x[] = (alpha[1,0] - alpha[-1,0]) / 2./Delta ;
#else
  gradb.x[] = (alpha[-2,0] - 8*alpha[-1,0]+ 8.* alpha[1,0] - alpha[2,0]) / 12./Delta ;
#endif

  
struct { double x, y; } a = {-1., 1.};
 
foreach()
    foreach_dimension()
#if order == 2
    rotc.x[] =  a.x * (beta[0,1] - beta[0,-1]) / 2./Delta ;
  #else
  rotc.x[] =  a.x * (beta[0,-2] - 8*beta[0,-1]+ 8.* beta[0,1] - beta[0,2]) / 12./Delta ;
#endif


  boundary ( (scalar *){gradb,rotc});


  /**
Compute h
  */

  foreach()
    foreach_dimension()
    h.x[] = omega.x[] - gradb.x[] - rotc.x[];

  boundary  ((scalar *){h});

  
    return s;  // Fix : better return 

}

/**
[2D Test](hodge.c)
*/
