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

We use a minmod limiter. */

static inline double minmod3 (double a, double b)
{
  if (a == b || a*b <= 0.)
    return 0.;
  return fabs(a) < fabs(b) ? a : b;
}

#define BGHOSTS 2

/**
## Time derivative

This function fills *dphi* with $\partial_t\phi$ .*/

static
void dphidt (scalar phi, scalar dphi, scalar phi0, double cfl)
{
  foreach() {
    double dt = cfl*Delta;
    
    /**
    We first calculate the inputs of the Hamiltonian
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
    
    coord gra[2];
    for (int i = 0, j = 1; i < 2; j = 1 - 2*++i)
      foreach_dimension() {
	double s1 = (phi[2*j] + phi[] - 2.*phi[j])/Delta; 
	double s2 = (phi[1] + phi[-1] - 2.*phi[])/Delta;
	gra[i].x = j*((phi[j] - phi[])/Delta - minmod3(s1, s2)/2.);
      }

    /**
    We check for interfacial cells. */
    
    bool interfacial = false;
    foreach_dimension()
      if (phi0[-1]*phi0[] < 0 || phi0[1]*phi0[] < 0)
	interfacial = true;
    
    dphi[] = - sign2(phi0[]);
    
    if (interfacial) {
      
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
	for (int i = 0, j = 1; i < 2; j = 1 - 2*++i)
	  if (phi0[]*phi0[j] < 0.) {

	    /**
	    We compute the subcell fix near the interface.

	    $$
	    \Delta x^+ = \left\{ \begin{array}{ll}
	    \Delta x \cdot \left( \dfrac{\phi^0_{i,j}-\phi^0_{i+1,j} -
	    \text{sgn}(\phi^0_{i,j}-\phi^0_{i+1,j})\sqrt{D}}{}\right) 
	    \text{ if } \left| \phi^0_{xx}\right| >\epsilon \		\
	    \Delta x \cdot \dfrac{\phi^0_{ij}}{\phi^0_{i,j}-\phi^0_{i+1,j}} \text{ else.}\ \
	    \end{array}
	    \right.
	    $$
	    with
	    $$
	    \phi_{xx}^0 = \text{minmod}(\phi^0_{i-1,j}-2\phi^0_{ij}+\phi^0_{i+1,j}, 
	    \phi^0_{i,j}-2\phi^0_{i+1j}+\phi^0_{i+2,j}) \		\
	    D = \left( \phi^0_{xx}/2  - \phi_{ij}^0 - \phi_{i+1,j} \right)^2 
                - 4\phi_{ij}^0\phi_{i+1,j}^0
	    $$
	    For the $\Delta x^-$ calculation, replace all the $+$ subscript by $-$, this
	    is dealt with properly with the `j` parameter below. */

	    double dx = Delta;
	    double phixx = minmod3 (phi0[2*j] + phi0[] - 2.*phi0[j],
				    phi0[1] + phi0[-1] - 2.*phi0[]);
	    if (fabs(phixx) > 1./HUGE) {
	      double D = sq(phixx/2. - phi0[] - phi0[j]) - 4.*phi0[]*phi0[j];
	      dx *= 1/2. + (phi0[] - phi0[j] - sign2(phi0[] - phi0[j])*sqrt(D))/phixx;
	    }
	    else
	      dx *= phi0[]/(phi0[] - phi0[j]);
	    
	    if (dx != 0.) {
	      double sxx1 = phi[2 - 4*i] + phi[] - 2.*phi[1 - 2*i];
	      double sxx2 = phi[1] + phi[-1] - 2.*phi[];
	      gra[i].x = (2*i - 1)*(phi[]/dx + dx*minmod3(sxx1, sxx2)/(2.*sq(Delta)));
	    }
	    else 
	      gra[i].x = 0.;
	    size = min(size, dx);
	  }
      dphi[] *= min(dt, fabs(size)/2.)/dt;
    }

    /**
    The Godunov Hamiltonian is
    $$
    H_G(a,b,c,d) = \left\{ \begin{array}{ll}
    \sqrt{\text{max}((a^-)^2,(b^+)^2 + (c^-)^2,(d^+)^2)} \text { when } \text{sign}(\phi^0_{ij})
    \geq 0\\
    \sqrt{\text{max}((a^+)^2,(b^-)^2 + (c^+)^2,(d^-)^2)} \text { when } \text{sign}(\phi^0_{ij}) < 0
    \end{array}
    \right.
    $$
    */
    
    double H_G = 0;
    if (phi0[] > 0)
      foreach_dimension() {
	double a = min(0., gra[0].x); 
	double b = max(0., gra[1].x);
	H_G += max(sq(a), sq(b));
      }
    else
      foreach_dimension() {
	double a = max(0., gra[0].x);
	double b = min(0., gra[1].x);
	H_G += max(sq(a), sq(b));
      }
    
    dphi[] *= sqrt(H_G) - 1.;
  }
}

/**
## The redistance() function */

trace
int redistance (scalar phi,
		int imax = 1,       // The maximum number of iterations
		double cfl = 0.5,   // The CFL number
		int order = 3,      // The order of time integration
		double eps = 1e-6,  // The maximum error on $|\nabla\phi| - 1$
		double band = HUGE, // The width of the band in which to compute the error
		scalar resf = {-1}) // The residual $|\nabla\phi| - 1$
{

  /**
  We create `phi0[]`, a copy of the initial level-set function before
  the iterations. */
  
  scalar phi0[];
  foreach()
    phi0[] = phi[] ;

  /**
  Time integration iteration loop. */
  
  for (int i = 1; i <= imax; i++) {
    double maxres = 0.;
    
    /**
    We use either a RK2 scheme... */

    if (order == 2) {
      scalar tmp1[];
      dphidt (phi, tmp1, phi0, cfl/2.);
      foreach()
	tmp1[] = phi[] + Delta*cfl/2.*tmp1[];
      scalar tmp2[];
      dphidt (tmp1, tmp2, phi0, cfl);
      foreach (reduction(max:maxres)) {
	double res = tmp2[];
	if (resf.i >= 0) resf[] = res;
	if (fabs (res) > maxres && fabs(phi0[]) < band*Delta) maxres = fabs (res);
	phi[] += Delta*cfl*tmp2[];
      }
    }
    
    /**
    ... or a RK3 compact scheme from [Shu and Osher,
    1988](#shu1988). */
    
    else {
      scalar tmp1[];
      dphidt (phi, tmp1, phi0, cfl);
      foreach()
	tmp1[] = phi[] + Delta*cfl*tmp1[];
      scalar tmp2[];
      dphidt (tmp1, tmp2, phi0, cfl);
      foreach()
	tmp1[] = (3.*phi[] + tmp1[] + Delta*cfl*tmp2[])/4.;
      dphidt (tmp1, tmp2, phi0, cfl);
      foreach (reduction(max:maxres)) {
	double res = 2./3.*((phi[] - tmp1[])/(Delta*cfl) - tmp2[]);
	if (resf.i >= 0) resf[] = res;
	if (fabs (res) > maxres && fabs(phi0[]) < band*Delta) maxres = fabs (res);
	phi[] = (phi[] + 2.*(tmp1[] + Delta*cfl*tmp2[]))/3.;
      }
    }
    
    if (maxres < eps)
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
  journal       = {Journal of Computational Physics},
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
  author        = {Chohong Min and Frédéric Gibou},
  title         = {A second order accurate level set method on non-graded adaptive cartesian grids},
  journal       = {Journal of Computational Physics},
  year          = {2007},
  volume        = {225},
  pages         = {300-321},
  issn          = {0021-9991},
  doi           = {10.1016/j.jcp.2006.11.034},
}

@article{min2010,
  author        = {Chohong Min},
  title         = {On reinitializing level set functions},
  journal       = {Journal of Computational Physics},
  year          = {2010},
  volume        = {229},
  pages         = {2764-2772},
  issn          = {0021-9991},
  doi           = {10.1016/j.jcp.2009.12.032},
}

@hal{limare2022, hal-03889680}
~~~
*/
