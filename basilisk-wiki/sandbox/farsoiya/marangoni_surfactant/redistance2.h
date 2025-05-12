/**
# Redistancing of a distance field

The original implementation was by [Alexandre
Limare](/sandbox/alimare/README) used in particular in [Limare et al.,
2022](#limare2022).

This file implements redistancing of a distance field with subcell
correction, see the work of [Russo & Smereka, 2000](#russo2000) with
corrections by [Min & Gibou, 2007](#min2007) and by [Min,
2010](#min2010).

Let $\phi$ be a function close to a signed function that has been
perturbed by numerical diffusion (more precisely, a non-zero
tangential velocity). By iterating on the eikonal equation
$$
\left\{\begin{array}{ll}
\phi_t + \text{sign}(\phi^{0}) \left(\left| \nabla \phi\right| - 1 \right) = 0\\ 
\phi(x,0) = \phi^0(x)
\end{array}
\right.
$$
we can correct or redistance $\phi$ to make it a signed function.

We use a Godunov Hamiltonian approximation for
$\left| \nabla \phi\right|$
$$
\left| \nabla \phi \right|_{ij} = H_G(D_x^+\phi_{ij}, D_x^-\phi_{ij}, 
                                      D_y^+\phi_{ij}, D_y^-\phi_{ij})
$$
where $D^\pm\phi_{ij}$ denotes the one-sided ENO difference finite
difference in the x- direction
$$
D_x^+ = \dfrac{\phi_{i+1,j}-\phi_{i,j}}{\Delta} - 
  \dfrac{\Delta}{2}\text{minmod}(D_{xx}\phi_{ij}, D_{xx}\phi_{i+1,j})
$$
$$
D_x^- = \dfrac{\phi_{i,j}-\phi_{i-1,j}}{\Delta} + 
  \dfrac{\Delta}{2}\text{minmod}(D_{xx}\phi_{ij}, D_{xx}\phi_{i+1,j})
$$
here $D_{xx}\phi_{ij} = (\phi_{i-1,j} - 2\phi{ij} + \phi_{i+1,j})/\Delta^2$.

The minmod function is zero when the two arguments have different
signs, and takes the argument with smaller absolute value when the two
have the same sign.

The Godunov Hamiltonian $H_G$ is given as
$$
H_G(a,b,c,d) = \left\{ \begin{array}{ll}
\sqrt{\text{max}((a^-)^2,(b^+)^2 + (c^-)^2,(d^+)^2)} \text { when } \text{sign}(\phi^0_{ij})
\geq 0\\
\sqrt{\text{max}((a^+)^2,(b^-)^2 + (c^+)^2,(d^-)^2)} \text { when } \text{sign}(\phi^0_{ij}) < 0
\end{array}
\right.
$$
with
$$
x^+ = \text{max}(0, x)\\
x^- = \text{min}(0, x)\\
$$

We use a minmod limiter and a second-order WENO approximation for the
derivatives. */

static inline double minmod3 (double a, double b)
{
  if (a == b || a*b <= 0.)
    return 0.;
  return fabs(a) < fabs(b) ? a : b;
}

#define BGHOSTS 2

foreach_dimension()
static inline double weno_diff_x (Point point, scalar s, int i)
{
  double s1 = (s[2*i] + s[] - 2*s[i])/Delta; 
  double s2 = (s[1] + s[-1] - 2*s[])/Delta;
  return i*((s[i] - s[])/Delta - minmod3(s1, s2)/2.);
}

/**
## Precalculation of the inputs of the Hamiltonian

$$
D^+\phi_{ij} = \dfrac{\phi_{i+1,j} - \phi_{ij}}{\Delta x} - \dfrac{\Delta}
{2}\text{minmod}(D_{xx}\phi_{ij},D_{xx}\phi_{i+1,j})
$$
$$
D^-\phi_{ij} = \dfrac{\phi_{i,j} - \phi_{i-1,j}}{\Delta x} - \dfrac{\Delta}
{2}\text{minmod}(D_{xx}\phi_{ij},D_{xx}\phi{i-1,j})
$$
where
$$
D_{xx}\phi_{ij} = \dfrac{\phi_{i-1,j} - 2\phi_{ij} + \phi_{i+1,j}}{\Delta x^2}
$$
*/

static
void prehamiltonian (Point point, coord  * grapl, coord * gramin, scalar s)
{
  foreach_dimension() {
    grapl->x  = weno_diff_x (point, s,  1);
    gramin->x = weno_diff_x (point, s, -1);
  }
}

/**
## Godunov Hamiltonian

$$
H_G(a,b,c,d) = \left\{ \begin{array}{ll}
\sqrt{\text{max}((a^-)^2,(b^+)^2 + (c^-)^2,(d^+)^2)} \text { when } \text{sign}(\phi^0_{ij})
\geq 0\\
\sqrt{\text{max}((a^+)^2,(b^-)^2 + (c^+)^2,(d^-)^2)} \text { when } \text{sign}(\phi^0_{ij}) < 0
\end{array}
\right.
$$
*/

static
double hamiltonian (Point point, scalar s0, coord grapl, coord gramin)
{
  double hamil = 0;
  if (s0[] > 0)
    foreach_dimension() {
      double a = min(0., grapl.x); 
      double b = max(0., gramin.x);
      hamil += max(sq(a), sq(b));
    }
  else
    foreach_dimension() {
      double a = max(0., grapl.x);
      double b = min(0., gramin.x);
      hamil += max(sq(a), sq(b));
    }
  return sqrt(hamil);
}

/**
## Root extraction for the subcell fix near the interface

$$
\Delta x^+ = \left\{ \begin{array}{ll}
\Delta x \cdot \left( \dfrac{\phi^0_{i,j}-\phi^0_{i+1,j} -
\text{sgn}(\phi^0_{i,j}-\phi^0_{i+1,j})\sqrt{D}}{}\right) 
\text{ if } \left| \phi^0_{xx}\right| >\epsilon \\
\Delta x \cdot \dfrac{\phi^0_{ij}}{\phi^0_{i,j}-\phi^0_{i+1,j}} \text{ else.}\\
\end{array}
\right.
$$
with
$$
\phi_{xx}^0 = \text{minmod}(\phi^0_{i-1,j}-2\phi^0_{ij}+\phi^0_{i+1,j}, 
                     \phi^0_{i,j}-2\phi^0_{i+1j}+\phi^0_{i+2,j}) \\
D = \left( \phi^0_{xx}/2  - \phi_{ij}^0 - \phi_{i+1,j} \right)^2  - 4\phi_{ij}^0\phi_{i+1,j}^0
$$
For the $\Delta x^-$ calculation, replace all the $+$ subscript by $-$, this
is dealt with properly with the `dir` variable in the following function. */

foreach_dimension()
static inline double root_x (Point point, scalar s, double eps, int dir)
{
  double phixx = minmod3 (s[2*dir] + s[] - 2*s[dir], s[1] + s[-1] - 2*s[]);
  if (fabs(phixx) > eps) {
    double D = sq(phixx/2. - s[] - s[dir]) - 4.*s[]*s[dir];
    return 1/2. + (s[] - s[dir] - sign2(s[] - s[dir])*sqrt(D))/phixx;
  }
  else
    return s[]/(s[]- s[dir]);
}

/**
## Forward Euler Integration

Simple Euler integration for the redistance() function. */

static
double forward_euler (scalar dist, scalar temp, scalar dist0, double dt)
{
  double res = 0.;
  foreach (reduction(max:res)) {

    coord gra[2];
    prehamiltonian (point, &gra[0], &gra[1], temp);

    bool interfacial = false;
    foreach_dimension()
      if (dist0[-1]*dist0[] < 0 || dist0[1]*dist0[] < 0)
	interfacial = true;
    
    if (!interfacial) {

      /**
      Far from the interface, we simply use the Hamiltonian defined earlier. */
      
      double delt = sign2(dist0[])*(hamiltonian (point, dist0, gra[0], gra[1]) - 1.);
      dist[] -= dt*delt;
      res = max (res, fabs(delt));
    }
    else {

      /**
      Near the interface, *i.e.* for cells where
      $$
      \phi^0_i\phi^0_{i+1} \leq 0 \text{ or } \phi^0_i\phi^0_{i-1} \leq 0
      $$

      The scheme must stay truly upwind, meaning that the movement of the 0
      level-set of the function must be as small as possible. Therefore the upwind
      numerical scheme is modified to
      $$
      D_x^+ = \dfrac{0-\phi_{ij}}{\Delta x^+} - \dfrac{\Delta x^+}{2} \text{minmod}(D_
      {xx}\phi_{ij},D_{xx}\phi_{i+1,j}) \text{ if } \phi_{ij}\phi_{i+1,j} < 0
      $$
      $$
      D_x^- = \dfrac{\phi_{ij}-0}{\Delta x^-} + \dfrac{\Delta x^-}{2} \text{minmod}(D_
      {xx}\phi_{ij},D_{xx}\phi_{i-1,j}) \text{ if } \phi_{ij}\phi_{i+1,j} < 0
      $$
      which is the correction by [Min & Gibou 2007](#min2007). */
      
      double size = HUGE;
      foreach_dimension()
	for (int i = 0; i < 2; i++)
	  if (dist0[]*dist0[1 - 2*i] < 0.) {
	    double dx = Delta*root_x (point, dist0, 1./HUGE, 1 - 2*i);
	    if (dx != 0.) {
	      double sxx1 = temp[2 - 4*i] + temp[] - 2.*temp[1 - 2*i];
	      double sxx2 = temp[1] + temp[-1] - 2.*temp[];
	      gra[i].x = (2*i - 1)*(temp[]/dx + dx*minmod3(sxx1, sxx2)/(2.*sq(Delta)));
	    }
	    else 
	      gra[i].x = 0.;
	    size = min(size, dx);
	  }
      double delt = sign2(dist0[])*(hamiltonian (point, dist0, gra[0], gra[1]) - 1.);
      dist[] -= min(dt, fabs(size)/2.)*delt;
      res = max (res, fabs(delt));
    }
  }
  return res;
}

/**
## The redistance() function

In 2D, if no specific timestep is set up by the user, we take the most
restrictive one with regards to the CFL condition 
$$
\Delta t = \Delta/2.
$$
*/

trace
int redistance (scalar dist, int imax = 1,
		double dt = L0/(1 << grid->maxdepth)/2.)
{

  /**
  Convergence is reached if residual is below $dt\times 10^{-6}$ */

  double eps = dt*1e-6;

  /**
  We create `dist0[]` which will be a copy of the initial level-set function
  before the iterations and `temp[]` which will be $\phi^{n}$ used for the
  iterations. */
  
  scalar dist0[];
  foreach()
    dist0[] = dist[] ;

  /**
  Time integration iteration loop.

  One can choose between Runge Kutta 2 and Forward Euler temporal integration. */
  
  for (int i = 1; i <= imax; i++) {

    /**
    ## RK3

    We use a Runge Kutta 3 compact version taken from [Shu and Osher,
    1988](#shu1988) made of 3 backward Euler steps.
    
    ### Step 1-2

    $$
    \frac{\widetilde{\phi}^{n+1}  - \phi^n}{\Delta t}  = \text{RHS}^n
    $$
    $$
    \dfrac{\widetilde{\phi}^{n+2}  - \widetilde{\phi}^{n+1}}{\Delta t}  
      = \widetilde{\text{RHS}}^{n+1} 
    $$
    with :
    $$
    \text{RHS} =  \text{sgn} (\phi_{ij}^0)\cdot \left[ H_G\left( D_x^+\phi_{ij}^n, 
                D_x^-\phi_{ij}^n, D_y^+\phi_{ij}^n,
		D_y^-\phi_{ij}^n \right)\right]
    $$
    */
    
    scalar temp[], temp1[], temp2[];
    foreach() {
      temp[] = dist[] ;
      temp1[] = dist[] ;
    }
    forward_euler (temp1, temp, dist0, dt);
    foreach()
      temp2[] = temp1[] ;
    forward_euler (temp2, temp1, dist0, dt);
    
    /**
    ### Intermediate value
    $$
    \widetilde{\phi}^{n+1/2} = 
      \dfrac{3}{4}\widetilde{\phi}^{n} + \dfrac{1}{4}\widetilde{\phi}^{n+2}
    $$
    */
    
    foreach() {
      temp1[] = (3.*dist[] + temp2[])/4.;
      temp2[] = temp1[];
    }

    /**
    ### Step 3
    $$
    \widetilde{\phi}^{n+3/2} - \widetilde{\phi}^{n+1/2} = \widetilde{RHS}^{n+1/2}
    $$
    */
    
    forward_euler (temp2, temp1, dist0, dt);
    
    /**
    ### Final Value
    $$
    \widetilde{\phi}^{n+1} = \widetilde{\phi}^{n} + \dfrac{2}{3}\widetilde{\phi}^{n+3/2}
    $$
    */

    double res = 0;
    foreach (reduction(max:res)) {
      res = max(res, 2./3.*fabs(dist[] - temp2[]));
      dist[] = (dist[] + 2.*temp2[])/3.;
    }

    /**
    Iterations are stopped when $L_1 = \text{max}(|\phi_i^{n+1}-\phi_i^n|) < \epsilon$ */

    if (res < eps)
      return i;
  }
  return imax;
}

/**
## References

~~~bib
@Article{shu1988,
  author        = {Chi-Wang Shu and Stanley Osher},
  title         = {Efficient implementation of essentially non-oscillatory shock-capturing schemes},
  year          = {1988},
  volume        = {77},
  pages         = {439-471},
  issn          = {0021-9991},
  doi           = {10.1016/0021-9991(88)90177-5},
}

@article{russo2000,
  title = {A remark on computing distance functions},
  volume = {163},
  number = {1},
  journal = {Journal of Computational Physics},
  author = {Russo, Giovanni and Smereka, Peter},
  year = {2000},
  pages = {51--67}
}

@article{min2007,
  author        = {Chohong Min and FrÃ©dÃ©ric Gibou},
  title         = {A second order accurate level set method on non-graded adaptive cartesian grids},
  year          = {2007},
  volume        = {225},
  pages         = {300-321},
  issn          = {0021-9991},
  doi           = {10.1016/j.jcp.2006.11.034},
}

@article{min2010,
  author        = {Chohong Min},
  title         = {On reinitializing level set functions},
  year          = {2010},
  volume        = {229},
  pages         = {2764-2772},
  issn          = {0021-9991},
  doi           = {10.1016/j.jcp.2009.12.032},
}

@hal{limare2022, hal-03889680}
~~~
*/
