/** # Introduction */
/** ## Preamble */
/**
Here is detailed how the algorithm of BCG is working, and how it is applied 
to self-similar space.
This code is specifically dedicated to the problem of 
[Lamb-Oseen](http://basilisk.fr/sandbox/cailler/lamb_oseen/lamb.c).
*/
/** ## Self-similar analytical solutions of the *Lamb-Oseen* vortex */
/**
Given the following definitions:

$$
\mathbf{\Xi}:=  
\begin{pmatrix}
\xi\\ 
\eta
\end{pmatrix}
\quad
,
\quad
\boldsymbol{\nabla}_{(\xi, \eta)} :=
\begin{pmatrix}
\partial_{\xi}\\ 
\partial_{\eta}
\end{pmatrix}
\quad
,
\quad
\mathbf{\Psi} :=
\begin{pmatrix}
\varphi\\ 
\psi
\end{pmatrix}
$$

where $\xi = x / \sqrt{\nu t}$, $\eta = y / \sqrt{\nu t}$ and $\tau = \ln (t)$. 
Then, expressions of velocities and pressure field in self-similar cartesian
coordinates are (by applying the *scale invariance method*):

$$
\left\{\begin{array}{lcl}
u(x,y,t) & = & \dfrac{\Gamma}{\sqrt{\nu t}} \varphi \left[ \xi(x,t), \eta(y,t), 
\tau(t) \right ] = \dfrac{\Gamma}{\sqrt{\nu t}}
\left[
-\dfrac{1}{2 \pi} \dfrac{\eta}{\xi^2 + \eta^2}  
\left(1 - \mathrm{e}^{-\dfrac{\xi^2 + \eta^2}{4}} \right )
\right ] \\ 
& & \\
v(x,y,t) & = & \dfrac{\Gamma}{\sqrt{\nu t}} \psi \left[ \xi(x,t), \eta(y,t), 
\tau(t) \right ] = \dfrac{\Gamma}{\sqrt{\nu t}}
\left[
\dfrac{1}{2 \pi} \dfrac{\xi}{\xi^2 + \eta^2}  
\left(1 - \mathrm{e}^{-\dfrac{\xi^2 + \eta^2}{4}} \right )
\right ] \\ 
& & \\
p_{\rho} := p/\rho & = & \dfrac{\Gamma^2}{\nu t} \gamma \left[ \xi(x,t), 
\eta(y,t), \tau(t) \right ] 
& & \\
& = & \dfrac{\Gamma ^2}{32 \pi^2 \text{$\nu $t}}
\left[
-2 \text{Ei}\left(-\dfrac{\eta ^2+\xi ^2}{2 }\right)
+ 2 \text{Ei}\left(-\dfrac{\eta ^2+\xi ^2}{4 }\right)
-\dfrac{4 \left(1-\mathrm{e}^{-\dfrac{\eta ^2 
+\xi ^2 }{4 }}\right)^2}{\eta ^2 + \xi ^2 }
\right]
\end{array}\right. 
$$

The self-similar expression of the pressure field was derived thanks to the
[German Wiki page](https://de.wikipedia.org/wiki/Hamel-Oseenscher-Wirbel#Druck).
It introduces the *exponential integral* function $\text{Ei}$.
*/
/** ## General Equation */
/**
The *Navier-Stokes* equations for the *Lamb-Oseen* problem, when we place 
ourselves in its relative self-similar space, read as:

$$
\left\{\begin{array}{rcl}
\partial_\tau \mathbf{\Psi}
+ 
\left[
\text{Re} \, (\mathbf{\Lambda} \cdot \boldsymbol{\nabla}) \mathbf{\Psi}
- 
\mathbf{\Psi}/2
 \right ] & = & - \text{Re} \, \boldsymbol{\nabla} \gamma
+
\boldsymbol{\nabla}^2 \mathbf{\Psi} \\
\boldsymbol{\nabla} \cdot \mathbf{\Psi} & = & 0 
\end{array}\right.
$$

where we have introduced:

* the *Reynolds number* defined by $\text{Re} = \Gamma / \nu$;
* the substitution $\mathbf{\Lambda} := \mathbf{\Psi} - \text{Re}^{-1}
\, \mathbf{\Xi} /2$;
* the notation simplification $\boldsymbol{\nabla} \equiv 
\boldsymbol{\nabla}_{(\xi, \eta)}$.

*/

/** # Self-Similar Bell-Collela-Glaz algorithm */
/** ## Splitting Strategy */
/**
### Splitting
The momentum equation can be, in a first step, approximated by:

$$
\dfrac{\mathbf{\Psi}^{n+1} - \mathbf{\Psi}^{n}}{\Delta \tau}
+
\left[
\text{Re} \, (\mathbf{\Lambda} \cdot \boldsymbol{\nabla}) \mathbf{\Psi}
- 
\mathbf{\Psi}/2
\right]^{n+1/2}
=
- \text{Re} \, \boldsymbol{\nabla} \gamma^n + \boldsymbol{\nabla}^2 \mathbf{\Psi}^{n+1} 
$$

Using linear stability upon the temporal term in regards to $\mathbf{\Psi}$, 
**splitting** can be applied:

$$
\left\{\begin{matrix}
& \dfrac{\mathbf{\Psi}^{*} - \mathbf{\Psi}^{n}}{\Delta \tau}
&+
\left[
\text{Re} \, (\mathbf{\Lambda} \cdot \boldsymbol{\nabla}) \mathbf{\Psi}
- 
\mathbf{\Psi}/2
\right]^{n+1/2} & = & 0\\ 
+ \\
& \dfrac{\mathbf{\Psi}^{**} - \mathbf{\Psi}^{*}}{\Delta \tau} &  & = 
& - \text{Re} \, \boldsymbol{\nabla} \gamma^n\\ 
+ \\
& \dfrac{\mathbf{\Psi}^{n+1} - \mathbf{\Psi}^{**}}{\Delta \tau} &  & = 
&\boldsymbol{\nabla}^2 \mathbf{\Psi}^{n+1}  
\end{matrix}\right.
$$

### Why evaluating the advection term at $n + 1/2$ 
We are looking for a second-order time scheme, in order to increase the 
precision of our computations. 

* For a given function $f$, a time-integration scheme between steps $n$ and 
$n + 1$ would have given, at first order:

$$
I = \displaystyle\int_{t}^{t + \Delta t} 
\left(
f(t) + (u - t)f'(t) + \mathcal{O}(\Delta t^2)
\right )
\,\, \mathrm{d}u
= I_n + \dfrac{\Delta t^2}{2} \, f'(t) 
+ \mathcal{O}(\Delta t^3)
$$

with $I_n = \Delta t \, f^n$. As it can be seen, the additional term 
$\Delta t^2 \, f'(t) / 2$ is a **1st-order error**, that is cumulated
for each timestep of the simulation. 

For $N = (t_f - t_0) / \Delta t$ the total number of timesteps,
then the *total error* at first-order will be:

$$
E_{tot}^{(1)}
=
(t_f - t_0) \dfrac{\Delta t}{2} \, f'(t)
$$

* *However*, by discretizing between $n$ and a provisional step $n+1/2$, that is 
to say, using a "*middle-point technique*" of discretization, the time 
integration becomes then, at first order:

$$
I = \displaystyle\int_{t}^{t + \Delta t} 
\left(
f(t + \Delta t /2) + (u - (t + \Delta t /2))f'(t + \Delta t /2) 
+ (u - (t + \Delta t /2))^2 \,\dfrac{f''(t + \Delta t /2)}{2}
+ \mathcal{O}(\Delta t^3)
\right )
\,\, \mathrm{d}u
$$

$$
I = I_{n+1/2} + 0 \times f'(t + \Delta t /2)
  + \dfrac{\Delta t^3}{24} \, f''(t)
  + \mathcal{O}(\Delta t^4)
=
I_{n+1/2} + \mathcal{O}(\Delta t^3)
$$

with $I_{n+1/2} = \Delta t \, f^{n+1/2}$.
Therefore, with such a time-scheme method, an order of magnitude in terms of 
error has been saved, meaning that the numerical error of discretization 
has been significantly improved.
Indeed, the *total error* at **second**-order will be:

$$
E_{tot}^{(2)}
=
(t_f - t_0) \dfrac{\Delta t^2}{24} \, f''(t)
$$

*/
/** ## Rewriting the advection term*/
/** 
Let $\mathbf{A}^{n+1/2} = \left[ \text{Re} \, (\mathbf{\Lambda} \cdot 
\boldsymbol{\nabla}) \mathbf{\Psi} - \mathbf{\Psi}/2 \right ]^{n+1/2}$ be the 
advection quantity. 

Since by definition $\lambda_j = \psi_j - {\text{Re}}^{-1} \, \xi_j / 2$, the fact that $\psi_{j,j} = 0$, 
and also that $\xi_{j,j} = 2$ in 2D, then:

$$
(\mathbf{\Lambda} \cdot \boldsymbol{\nabla}) \mathbf{\Psi} 
=
\lambda_j \psi_{i,j} 
= 
(\lambda_j \, \psi_i)_{,j} 
- 
\psi_i \, \lambda_{j,j}
=
(\lambda_j \, \psi_i)_{,j} 
- 
\psi_i \,
(\psi_{j,j} - {\text{Re}}^{-1} \, \xi_{j,j} / 2)
=
(\lambda_j \, \psi_i)_{,j} 
+ 
\psi_i/\text{Re}
$$

$$
\Rightarrow
\text{Re} \, (\mathbf{\Lambda} \cdot \boldsymbol{\nabla}) \mathbf{\Psi} 
=
\text{Re} \, (\lambda_j \, \psi_i)_{,j} 
+ 
\psi_i
$$

$$
\Rightarrow
\mathbf{A}
=
\text{Re} \, (\mathbf{\Lambda} \cdot \boldsymbol{\nabla}) \mathbf{\Psi} 
- \mathbf{\Psi} / 2
=
\text{Re} \, (\lambda_j \, \psi_i)_{,j} 
+ 
\psi_i / 2
$$

By the use of the *Green-Ostrogradsky theorem*, the integration over the 
volume of a cell $\mathcal{C}$ of $\mathbf{A}^{n+1/2}$ is:

$$
\displaystyle\int_{\mathcal{C}}{a_i}^{n+1/2} \,\, \mathrm{d}V
=
\text{Re} 
\displaystyle\int_{\mathcal{\partial C}} \psi_i (\lambda_j n_j) \,\, \mathrm{d}S 
+ \dfrac{1}{2} \displaystyle\int_{\mathcal{C}}\psi_i \,\, \mathrm{d}V
$$

that is:

$$
\displaystyle\int_{\mathcal{C}}{a_i}^{n+1/2} \,\, \mathrm{d}V
=
\text{Re} 
\displaystyle\int_{\mathcal{\partial C}} \psi_i (\psi_j n_j) \,\, \mathrm{d}S 
- \dfrac{1}{2} 
\displaystyle\int_{\mathcal{\partial C}} \psi_i (\xi_j n_j) \,\, \mathrm{d}S 
+ \dfrac{1}{2} 
\displaystyle\int_{\mathcal{C}}\psi_i \,\, \mathrm{d}V
$$
*/
/** ## Approximation in 2D of the integrals*/
/**
Now, for each square-cell of size $h$, and for each direction $d$ (north, 
east, south, west) we sum their value to approximate numerically the integral
terms over cell boundaries; for instance:

$$
\displaystyle\int_{\mathcal{\partial C}} \psi_i (\psi_j n_j) \,\, \mathrm{d}S
\approx 
h \, \sum_d \left.{\psi_i}^{n+1/2}\right|_d 
\left.\left({\psi_j}^{n+1/2} \, n_j\right)\right|_d
$$

As for the volumic term:

$$
\displaystyle\int_{\mathcal{C}} \psi_i \,\, \mathrm{d}V
\approx 
h^2 \, {\psi_i}^n
$$
*/
/** ## Taylor-expansion in space and time of $\left.{\psi_i}^{n+1/2}\right|_d$ */






/** # Basilisk-C Code */
/**
## Self-Similar Bell-Collela-Glaz advection scheme

The function below implements the 2nd-order, unsplit, upwind scheme of
[Bell-Collela-Glaz, 1989](references.bib#bell89), adapted for the 
self-similar case of the *Lamb-Oseen* vortex.
Given a centered scalar field *Ψ*, a face vector field *λf* 
(possibly weighted by a face metric), a timestep *dτ* and a source term field 
*src*, it fills the face vector field *flux* with the components of the 
advection fluxes of *Ψ*. 

*Old*      *New*  *Old*      *New*  *Old*      *New*  
-------  -------  -------  -------  -------  ------- 
x              ξ  un            λn  f              Ψ 
y              η  vn            λt  f2            Ψ2 
z              ζ  wn            λb  fyy          Ψηη 
dt            dτ  uf            λf  fzz          Ψζζ 

Table: **NEW NOTATIONS**
*/

#include "selfsim-coord.h"


void tracer_fluxes_selfsim (scalar Ψ,
		    face vector λf,
		    face vector flux,
                     face vector Rey,
		    double dτ,
		    (const) scalar src
        )
{

  /**
  We first compute the cell-centered gradient of *Ψ* in a locally-allocated
  vector field. */
  
  vector grad_Ψ[];
  gradients ({Ψ}, {grad_Ψ});

  /**
  For each face, the flux is composed of two parts... */

  foreach_face() {

    /**
    A normal component... (Note that we cheat a bit here, `λn` should
    strictly be `dτ*(λf.ξ[i] + λf.ξ[i+1])/((fm.ξ[] +
    fm.ξ[i+1])*Delta)` but this causes trouble with boundary
    conditions (when using narrow '1 ghost cell' stencils)). */

    double λn = dτ*λf.ξ[]/(fm.ξ[]*Delta + SEPS);
    double s = sign(λn);
    int i = -(s + 1.)/2.;
    double Ψ2 = Ψ[i] + (src[] + src[-1])*dτ/4. 
      + s*(1. - s*λn)*grad_Ψ.ξ[i]*Delta/2.;

    /**
    and tangential components... */

    #if dimension > 1
    if (fm.η[i] && fm.η[i,1]) {
      double λt = (λf.η[i] + λf.η[i,1])/(fm.η[i] + fm.η[i,1]);
      double Ψηη = λt < 0. ? Ψ[i,1] - Ψ[i] : Ψ[i] - Ψ[i,-1];
      Ψ2 -= dτ*λt*Ψηη/(2.*Delta);
    }
    #endif
    #if dimension > 2
    if (fm.ζ[i] && fm.ζ[i,0,1]) {
      double λb = (λf.ζ[i] + λf.ζ[i,0,1])/(fm.ζ[i] + fm.ζ[i,0,1]);
      double Ψζζ = λb < 0. ? Ψ[i,0,1] - Ψ[i] : Ψ[i] - Ψ[i,0,-1];
      Ψ2 -= dτ*λb*Ψζζ/(2.*Delta);
    }
    #endif

    flux.ξ[] = Rey.ξ[]*(Ψ2*λf.ξ[]);
  }
}

/**
The function below uses the *tracer_fluxes_selfsim* function to integrate the
advection equation, using an explicit scheme with timestep *dτ*, for
each tracer in the list. */

struct Advection {
  scalar * tracers;
  face vector λf;
  face vector Rey;
  double dτ;
  scalar * src; // mandatory here
};

void advection_selfsim (struct Advection p)
{

  /**
  If *src* is not provided we set all the source terms to zero. */
  
  scalar * lsrc = p.src;
  if (!lsrc)
    for (scalar s in p.tracers)
      lsrc = list_append (lsrc, zeroc);
  assert (list_len(p.tracers) == list_len(lsrc));

  scalar Ψ, src;
  for (Ψ,src in p.tracers,lsrc) {
    face vector flux[];
    tracer_fluxes_selfsim (Ψ, p.λf, flux, p.Rey, p.dτ, src);
    foreach()
      foreach_dimension()
        Ψ[] += p.dτ*(flux.ξ[] - flux.ξ[1])/(Delta*cm[]);
  }

  if (!p.src)
    free (lsrc);
}