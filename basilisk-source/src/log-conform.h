/**
# The log-conformation method for some viscoelastic constitutive models

## Introduction

Viscoelastic fluids exhibit both viscous and elastic behaviour when
subjected to deformation. Therefore these materials are governed by
the Navier--Stokes equations enriched with an extra *elastic* stress
$\mathbf{\tau}_p$
$$
\rho\left[\partial_t\mathbf{u}+\nabla\cdot(\mathbf{u}\otimes\mathbf{u})\right] = 
- \nabla p + \nabla\cdot(2\mu_s\mathbf{D}) + \nabla\cdot\mathbf{\tau}_p
+ \rho\mathbf{a}
$$
where $\mathbf{D}=[\nabla\mathbf{u} + (\nabla\mathbf{u})^T]/2$ is the
deformation tensor and $\mu_s$ is the solvent viscosity of the
viscoelastic fluid.

The *polymeric* stress $\mathbf{\tau}_p$ represents memory effects due to
the polymers. Several constitutive rheological models are available in
the literature where the polymeric stress $\mathbf{\tau}_p$ is typically a 
function $\mathbf{f_s}(\cdot)$ of the conformation tensor $\mathbf{A}$ such as
$$
\mathbf{\tau}_p = \frac{\mu_p \mathbf{f_s}(\mathbf{A})}{\lambda}
$$
where $\lambda$ is the relaxation parameter and $\mu_p$ is the
polymeric viscosity.

The conformation tensor $\mathbf{A}$ is related to the deformation of
the polymer chains. $\mathbf{A}$ is governed by the equation
$$
D_t \mathbf{A} - \mathbf{A} \cdot \nabla \mathbf{u} - \nabla
\mathbf{u}^{T} \cdot \mathbf{A} =
-\frac{\mathbf{f_r}(\mathbf{A})}{\lambda} 
$$
where $D_t$ denotes the material derivative and
$\mathbf{f_r}(\cdot)$ is the relaxation function.

In the case of an Oldroyd-B viscoelastic fluid, $\mathbf{f}_s
(\mathbf{A}) = \mathbf{f}_r (\mathbf{A}) = \mathbf{A} -\mathbf{I}$,
and the above equations can be combined to avoid the use of
$\mathbf{A}$
$$
\mathbf{\tau}_p + \lambda (D_t \mathbf{\tau}_p -
\mathbf{\tau}_p \cdot \nabla \mathbf{u} -
\nabla \mathbf{u}^{T} \cdot \mathbf{\tau}_p)  = 2 \mu_p \mathbf{D}
$$

[Comminal et al. (2015)](#comminal2015) gathered the functions
$\mathbf{f}_s (\mathbf{A})$ and $\mathbf{f}_r (\mathbf{A})$ for
different constitutive models. In the present library we have
implemented the Oldroyd-B model and the related FENE-P model for which
$$
\mathbf{f}_s (\mathbf{A}) = \mathbf{f}_r (\mathbf{A}) =
\frac{\mathbf{A}}{1-Tr(\mathbf{A})/L^2} -\mathbf{I}
$$

## Parameters

The primary parameters are the retardation or relaxation time
$\lambda$ and the polymeric viscosity $\mu_p$. The solvent viscosity
$\mu_s$ is defined in the [Navier-Stokes
solver](navier-stokes/centered.h). */

(const) scalar lambda[] = 1.;
(const) scalar mup[] = 1.;

/**
Constitutive models other than Oldroyd-B (the default) are defined
through the two functions $\mathbf{f}_s (\mathbf{A})$ and
$\mathbf{f}_r (\mathbf{A})$. */

void (* f_s) (double, double *, double *) = NULL;
void (* f_r) (double, double *, double *) = NULL;

/**
## The log conformation approach

The numerical resolution of viscoelastic fluid problems often faces the
[High-Weissenberg Number
Problem](http://www.ma.huji.ac.il/~razk/iWeb/My_Site/Research_files/Visco1.pdf). 
This is a numerical instability appearing when strongly elastic flows
create regions of high stress and fine features. This instability
poses practical limits to the values of the relaxation time of the
viscoelastic fluid, $\lambda$.  [Fattal \& Kupferman (2004,
2005)](#fattal2004) identified the exponential nature of the solution
as the origin of the instability. They proposed to use the logarithm
of the conformation tensor $\Psi = \log \, \mathbf{A}$ rather than the
viscoelastic stress tensor to circumvent the instability.

The constitutive equation for the log of the conformation tensor is
$$ 
D_t \Psi = (\Omega \cdot \Psi -\Psi \cdot \Omega) + 2 \mathbf{B} +
\frac{e^{-\Psi} \mathbf{f}_r (e^{\Psi})}{\lambda}
$$
where $\Omega$ and $\mathbf{B}$ are tensors that result from the
decomposition of the transpose of the tensor gradient of the
velocity
$$ 
(\nabla \mathbf{u})^T = \Omega + \mathbf{B} + N
\mathbf{A}^{-1} 
$$ 

The antisymmetric tensor $\Omega$ requires only the memory of a scalar
in 2D since,
$$ 
\Omega = \left( 
\begin{array}{cc}
0 & \Omega_{12} \\
-\Omega_{12} & 0
\end{array} 
\right)
$$
The log-conformation tensor, $\Psi$, is related to the
polymeric stress tensor $\mathbf{\tau}_p$, by the strain function 
$\mathbf{f}_s (\mathbf{A})$
$$ 
\Psi = \log \, \mathbf{A} \quad \mathrm{and} \quad \mathbf{\tau}_p =
\frac{\mu_p}{\lambda} \mathbf{f}_s (\mathbf{A})
$$
where $Tr$ denotes the trace of the tensor and $L$ is an additional
property of the viscoelastic fluid.

We will use the Bell--Collela--Glaz scheme to advect the log-conformation 
tensor $\Psi$. */

#include "bcg.h"

/**
## Variables

The main variable will be the stress tensor $\mathbf{\tau}_p$. The trace of
the conformation tensor, $\mathbf{A}$, is often necessary for
constitutive viscoelastic models other than Oldroyd-B.  */

symmetric tensor tau_p[];
#if AXI
scalar tau_qq[];
#endif
(const) scalar trA[] = 0.;

event defaults (i = 0) {
  if (is_constant (a.x))
    a = new face vector;
  if (f_s || f_r)
    trA = new scalar;

  foreach() {
    foreach_dimension()
      tau_p.x.x[] = 0.;
    tau_p.x.y[] = 0.;
#if AXI
    tau_qq[] = 0;
#endif
  }

  /**
  ## Boundary conditions

  By default we set a zero Neumann boundary condition for all
  the components except if the bottom is an axis of symmetry. */

  for (scalar s in {tau_p}) {
    s.v.x.i = -1; // just a scalar, not the component of a vector
    foreach_dimension()
      if (s.boundary[left] != periodic_bc) {
	s[left] = neumann(0);
	s[right] = neumann(0);
      }
  }
#if AXI
  scalar s = tau_p.x.y;
  s[bottom] = dirichlet (0.);  
#endif  
}

/**
## Numerical Scheme 

The first step is to implement a routine to calculate the eigenvalues
and eigenvectors of the conformation tensor $\mathbf{A}$.

These structs ressemble Basilisk vectors and tensors but are just
arrays not related to the grid. */

typedef struct { double x, y;}   pseudo_v;
typedef struct { pseudo_v x, y;} pseudo_t;

static void diagonalization_2D (pseudo_v * Lambda, pseudo_t * R, pseudo_t * A)
{
  /**
  The eigenvalues are saved in vector $\Lambda$ computed from the
  trace and the determinant of the symmetric conformation tensor
  $\mathbf{A}$. */

  if (sq(A->x.y) < 1e-15) {
    R->x.x = R->y.y = 1.;
    R->y.x = R->x.y = 0.;
    Lambda->x = A->x.x; Lambda->y = A->y.y;
    return;
  }

  double T = A->x.x + A->y.y; // Trace of the tensor
  double D = A->x.x*A->y.y - sq(A->x.y); // Determinant

  /**
  The eigenvectors, $\mathbf{v}_i$ are saved by columns in tensor
  $\mathbf{R} = (\mathbf{v}_1|\mathbf{v}_2)$. */

  R->x.x = R->x.y = A->x.y;
  R->y.x = R->y.y = -A->x.x;
  double s = 1.;
  for (int i = 0; i < dimension; i++) {
    double * ev = (double *) Lambda;
    ev[i] = T/2 + s*sqrt(sq(T)/4. - D);
    s *= -1;
    double * Rx = (double *) &R->x;
    double * Ry = (double *) &R->y;
    Ry[i] += ev[i];
    double mod = sqrt(sq(Rx[i]) + sq(Ry[i]));
    Rx[i] /= mod;
    Ry[i] /= mod;
  }
}

/**
The stress tensor depends on previous instants and has to be
integrated in time. In the log-conformation scheme the advection of
the stress tensor is circumvented, instead the conformation tensor,
$\mathbf{A}$ (or more precisely the related variable $\Psi$) is
advanced in time.

In what follows we will adopt a scheme similar to that of [Hao \& Pan
(2007)](#hao2007). We use a split scheme, solving successively

a) the upper convective term:
$$
\partial_t \Psi = 2 \mathbf{B} + (\Omega \cdot \Psi -\Psi \cdot \Omega)
$$
b) the advection term:
$$ 
\partial_t \Psi + \nabla \cdot (\Psi \mathbf{u}) = 0
$$
c) the model term (but set in terms of the conformation 
tensor $\mathbf{A}$). In an Oldroyd-B viscoelastic fluid, the model is
$$ 
\partial_t \mathbf{A} = -\frac{\mathbf{f}_r (\mathbf{A})}{\lambda}
$$

The implementation below assumes that the values of $\Psi$ and
$\tau_p$ are never needed simultaneously. This means that $\tau_p$ can
be used to store (temporarily) the values of $\Psi$ (i.e. $\Psi$ is
just an alias for $\tau_p$). */

event tracer_advection (i++)
{
  tensor Psi = tau_p;
#if AXI
  scalar Psiqq = tau_qq;
#endif

  /**
  ### Computation of $\Psi = \log \mathbf{A}$ and upper convective term */

  foreach() {
    if (lambda[] == 0.) {
      foreach_dimension()
	Psi.x.x[] = 0.;
      Psi.x.y[] = 0.;
#if AXI
      Psiqq[] = 0.;
#endif
    }
    else { // lambda[] != 0.

      /**
      We assume that the stress tensor $\mathbf{\tau}_p$ depends on the
      conformation tensor $\mathbf{A}$ as follows
      $$
      \mathbf{\tau}_p = \frac{\mu_p}{\lambda} f_s (\mathbf{A}) = 
      \frac{\mu_p}{\lambda} \eta (\nu \mathbf{A} - I)
      $$
      In most of the viscoelastic models, $\nu$ and $\eta$ are 
      nonlinear parameters that depend on the trace of the conformation tensor,
      $\mathbf{A}$.*/

      double eta = 1., nu = 1.;
      if (f_s)
	f_s (trA[], &nu, &eta);

      double fa = (mup[] != 0 ? lambda[]/(mup[]*eta) : 0.);

      pseudo_t A;
      A.x.y = fa*tau_p.x.y[]/nu;
      foreach_dimension()
	A.x.x = (fa*tau_p.x.x[] + 1.)/nu;

      /**
      In the axisymmetric case, $\Psi_{\theta \theta} = \log A_{\theta
      \theta}$. Therefore $\Psi_{\theta \theta} = \log [ ( 1 + fa 
      \tau_{p_{\theta \theta}})/\nu]$. */

#if AXI
      double Aqq = (1. + fa*tau_qq[])/nu;
      Psiqq[] = log (Aqq); 
#endif

      /**
      The conformation tensor is diagonalized through the
      eigenvector tensor $\mathbf{R}$ and the eigenvalues diagonal
      tensor, $\Lambda$. */

      pseudo_v Lambda;
      pseudo_t R;
      diagonalization_2D (&Lambda, &R, &A);
      
      /**
      $\Psi = \log \mathbf{A}$ is easily obtained after diagonalization, 
      $\Psi = R \cdot \log(\Lambda) \cdot R^T$. */
      
      Psi.x.y[] = R.x.x*R.y.x*log(Lambda.x) + R.y.y*R.x.y*log(Lambda.y);
      foreach_dimension()
	Psi.x.x[] = sq(R.x.x)*log(Lambda.x) + sq(R.x.y)*log(Lambda.y);
      
      /**
      We now compute the upper convective term $2 \mathbf{B} +
      (\Omega \cdot \Psi -\Psi \cdot \Omega)$.
	
      The diagonalization will be applied to the velocity gradient
      $(\nabla u)^T$ to obtain the antisymmetric tensor $\Omega$ and
      the traceless, symmetric tensor, $\mathbf{B}$. If the conformation
      tensor is $\mathbf{I}$, $\Omega = 0$ and $\mathbf{B}= \mathbf{D}$.  */

      pseudo_t B;
      double OM = 0.;
      if (fabs(Lambda.x - Lambda.y) <= 1e-20) {
	B.x.y = (u.y[1,0] - u.y[-1,0] +
		 u.x[0,1] - u.x[0,-1])/(4.*Delta); 
	foreach_dimension() 
	  B.x.x = (u.x[1,0] - u.x[-1,0])/(2.*Delta);
      }
      else {
	pseudo_t M;
	foreach_dimension() {
	  M.x.x = (sq(R.x.x)*(u.x[1] - u.x[-1]) +
		   sq(R.y.x)*(u.y[0,1] - u.y[0,-1]) +
		   R.x.x*R.y.x*(u.x[0,1] - u.x[0,-1] + 
				u.y[1] - u.y[-1]))/(2.*Delta);
	  M.x.y = (R.x.x*R.x.y*(u.x[1] - u.x[-1]) + 
		   R.x.y*R.y.x*(u.y[1] - u.y[-1]) +
		   R.x.x*R.y.y*(u.x[0,1] - u.x[0,-1]) +
		   R.y.x*R.y.y*(u.y[0,1] - u.y[0,-1]))/(2.*Delta);
	}
	double omega = (Lambda.y*M.x.y + Lambda.x*M.y.x)/(Lambda.y - Lambda.x);
	OM = (R.x.x*R.y.y - R.x.y*R.y.x)*omega;
	
	B.x.y = M.x.x*R.x.x*R.y.x + M.y.y*R.y.y*R.x.y;
	foreach_dimension()
	  B.x.x = M.x.x*sq(R.x.x)+M.y.y*sq(R.x.y);	
      }

      /**
      We now advance $\Psi$ in time, adding the upper convective
      contribution. */

      double s = - Psi.x.y[];
      Psi.x.y[] += dt*(2.*B.x.y + OM*(Psi.y.y[] - Psi.x.x[]));
      foreach_dimension() {
	s *= -1;
	Psi.x.x[] += dt*2.*(B.x.x + s*OM);
      }

      /**
      In the axisymmetric case, the governing equation for $\Psi_{\theta
      \theta}$ only involves that component, 
      $$ 
      \Psi_{\theta \theta}|_t - 2 L_{\theta \theta} = 
      \frac{\mathbf{f}_r(e^{-\Psi_{\theta \theta}})}{\lambda} 
      $$
      with $L_{\theta \theta} = u_y/y$. Therefore step (a) for
      $\Psi_{\theta \theta}$ is */

#if AXI
      Psiqq[] += dt*2.*u.y[]/y;
#endif
    }
  }
  
  /**
  ### Advection of $\Psi$
  
  We proceed with step (b), the advection of the log of the
  conformation tensor $\Psi$. */

#if AXI
  advection ({Psi.x.x, Psi.x.y, Psi.y.y, Psiqq}, uf, dt);
#else
  advection ({Psi.x.x, Psi.x.y, Psi.y.y}, uf, dt);
#endif

  /**
  ### Model term */
  
  foreach() {
    if (lambda[] == 0.) {

      /**
      If $\lambda = 0$ the stress tensor for the polymeric part
      reduces to that of a Newtonian fluid $\mathbf{\tau}_p = 2 \mu_p
      \mathbf{D}$ with $\mathbf{D}$ the rate-of-strain
      tensor. Note that $\mathbf{\tau}_p$ is in this case independent of
      time. */

      foreach_dimension()
	tau_p.x.x[] = mup[]*(u.x[1,0] - u.x[-1,0])/Delta; // 2*mu*dxu;
      tau_p.x.y[] = mup[]*(u.y[1,0] - u.y[-1,0] +
			   u.x[0,1] - u.x[0,-1])/(2.*Delta); // mu*(dxv+dyu)
#if AXI
      tau_qq[] = 2.*mup[]*u.y[]/y;
#endif
    }
    else { // lambda != 0.
      
      /**
      It is time to undo the log-conformation, again by
      diagonalization, to recover the conformation tensor $\mathbf{A}$
      and to perform step (c).*/

      pseudo_t A = {{Psi.x.x[], Psi.x.y[]}, {Psi.y.x[], Psi.y.y[]}}, R;
      pseudo_v Lambda;
      diagonalization_2D (&Lambda, &R, &A);
      Lambda.x = exp(Lambda.x), Lambda.y = exp(Lambda.y);
      
      A.x.y = R.x.x*R.y.x*Lambda.x + R.y.y*R.x.y*Lambda.y;
      foreach_dimension()
	A.x.x = sq(R.x.x)*Lambda.x + sq(R.x.y)*Lambda.y;
#if AXI
      double Aqq = exp(Psiqq[]);
#endif

      /**
      We perform now step (c) by integrating 
      $\mathbf{A}_t = -\mathbf{f}_r (\mathbf{A})/\lambda$ to obtain
      $\mathbf{A}^{n+1}$. This step is analytic,
      $$
      \int_{t^n}^{t^{n+1}}\frac{d \mathbf{A}}{\mathbf{I}- \nu \mathbf{A}} = 
      \frac{\eta \, \Delta t}{\lambda}
      $$
      */

      double eta = 1., nu = 1.;
      if (f_r) {
#if 0 // Set to one if the midstep trace is to be used.
	scalar t = trA;
	t[] = A.x.x + A.y.y;
#if AXI
	t[] += Aqq;
#endif
#endif
	f_r (trA[], &nu, &eta);
      }

      double fa = exp(-nu*eta*dt/lambda[]);

#if AXI
      Aqq = (1. - fa)/nu + fa*exp(Psiqq[]);
      Psiqq[] = log (Aqq);
#endif

      A.x.y *= fa;
      foreach_dimension()
	A.x.x = (1. - fa)/nu + A.x.x*fa;

      /**
      The trace at time $n+1$ is also needed for some models. */
      
      if (f_s || f_r) {
	scalar t = trA;
	t[] = A.x.x + A.y.y;
#if AXI
	t[] += Aqq;
#endif
      }

      /**
      Then the stress tensor $\mathbf{\tau}_p^{n+1}$ is computed from
      $\mathbf{A}^{n+1}$ according to the constitutive model,
      $\mathbf{f}_s(\mathbf{A})$.  */

      nu = 1; eta = 1.;
      if (f_s)
	f_s (trA[], &nu, &eta);

      fa = mup[]/lambda[]*eta;
      
      tau_p.x.y[] = fa*nu*A.x.y;
#if AXI
      tau_qq[] = fa*(nu*Aqq - 1.);
#endif
      foreach_dimension()
	tau_p.x.x[] = fa*(nu*A.x.x - 1.);
    }
  }
}

/**
### Divergence of the viscoelastic stress tensor

The viscoelastic stress tensor $\mathbf{\tau}_p$ is defined at cell centers
while the corresponding force (acceleration) will be defined at cell
faces. Two terms contribute to each component of the momentum
equation. For example the $x$-component in Cartesian coordinates has
the following terms: $\partial_x \mathbf{\tau}_{p_{xx}} + \partial_y
\mathbf{\tau}_{p_{xy}}$. The first term is easy to compute since it can be
calculated directly from center values of cells sharing the face. The
other one is harder. It will be computed from vertex values. The
vertex values are obtained by averaging centered values.  Note that as
a result of the vertex averaging cells `[]` and `[-1,0]` are not
involved in the computation of shear. */

event acceleration (i++)
{
  face vector av = a;
  foreach_face()
    if (fm.x[] > 1e-20) {
      double shear = (tau_p.x.y[0,1]*cm[0,1] + tau_p.x.y[-1,1]*cm[-1,1] -
		      tau_p.x.y[0,-1]*cm[0,-1] - tau_p.x.y[-1,-1]*cm[-1,-1])/4.;
      av.x[] += (shear + cm[]*tau_p.x.x[] - cm[-1]*tau_p.x.x[-1])*
	alpha.x[]/(sq(fm.x[])*Delta);
    }
#if AXI
  foreach_face(y)
    if (y > 0.)
      av.y[] -= (tau_qq[] + tau_qq[0,-1])*alpha.y[]/sq(y)/2.;
#endif
}

/**
## References

~~~bib
@article{fattal2004,
  title={Constitutive laws for the matrix-logarithm of the conformation tensor},
  author={Fattal, Raanan and Kupferman, Raz},
  journal={Journal of Non-Newtonian Fluid Mechanics},
  volume={123},
  number={2-3},
  pages={281--285},
  year={2004},
  publisher={Elsevier}
}

@article{fattal2005,
  title={Time-dependent simulation of viscoelastic flows at 
         high {W}eissenberg number using the log-conformation representation},
  author={Fattal, Raanan and Kupferman, Raz},
  journal={Journal of Non-Newtonian Fluid Mechanics},
  volume={126},
  number={1},
  pages={23--37},
  year={2005},
  publisher={Elsevier}
}

@article{hao2007,
  title={Simulation for high {W}eissenberg number: viscoelastic 
         flow by a finite element method},
  author={Hao, Jian and Pan, Tsorng-Whay},
  journal={Applied mathematics letters},
  volume={20},
  number={9},
  pages={988--993},
  year={2007},
  publisher={Elsevier}
}

@article{comminal2015,
  title={Robust simulations of viscoelastic flows at high {W}eissenberg 
         numbers with the streamfunction/log-conformation formulation},
  author={Comminal, Rapha{\"e}l and Spangenberg, Jon and Hattel, Jesper Henri},
  journal={Journal of Non-Newtonian Fluid Mechanics},
  volume={223},
  pages={37--61},
  year={2015},
  publisher={Elsevier}
}
~~~

## See also

* [Functions $f_s$ and $f_r$ for the FENE-P model](fene-p.h)
*/
