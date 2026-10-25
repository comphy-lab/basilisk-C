**THE CONTAINS OF THIS FILE WAS GENERATED WITH AN AI (ChatGPT)** 
**TO DERIVE ALL THE NEEDED EQUATIONS TO ESTABLISH THE SELF-SIMILAR** 
**VISCOUS OPERATOR FOR THE 1st HUPPERT PROBLEM.**

**ALL THE CALCULATIONS WERE THEN MANUALLY CHECKED AND VERIFIED CAREFULLY**

Our goal is to prove the derivation of the self-similar viscous operator 
that is needed to implement to solve the Huppert's heap problem directly into 
the self-similar space, that is, how to obtain:

$$ \boxed{ \mathcal D_U =\frac{1}{R}\left[ A_x\partial_\xi(2MU_\xi) +A_y\partial_\eta(MU_\eta) +A_x\partial_\eta(MV_\xi) \right]. } $$

$$ \boxed{ \mathcal D_V =\frac{1}{R}\left[ A_x\partial_\xi(MV_\xi) +A_y\partial_\xi(MU_\eta) +A_y\partial_\eta(2MV_\eta) \right]. } $$

with:
$$ 
\boxed{ 
  \begin{aligned} 
    A_x(\tau)
    &=
    \dfrac{\nu_1 T_0}{L_0^2} 
    \mathrm{e}^{3\tau/5},
    \\ 
    A_y(\tau)
    &=
    \dfrac{\nu_1 T_0}{H_0^2}
    \mathrm{e}^{7\tau/5}. 
  \end{aligned} 
} 
$$
and here with the dimensionless density and viscosity fields 
$$ \rho=\rho_1 R,\qquad \mu=\mu_1 M, $$
where $\rho_1$ and $\mu_1$ are reference values.

--------------------------------------------------------------------------------
# Establishment of the self-similar viscous operator of the Huppert's problem


The key is to track each term in the Newtonian viscous stress divergence separately. 
Because the self-similar coordinates contract in one direction and expand in the other, 
the derivatives in $\xi$ and $\eta$ acquire different time-dependent prefactors. 
We must also distinguish the two contributions to the shear stress: 
$\partial_y(\mu \, \partial_y u)$ and $\partial_y(\mu \, \partial_x v)$, 
for example, do not have the same prefactor after transformation.

We will derive both components from the physical equations, keeping every scale 
factor visible until the final step. 

## 1. Start from the physical viscous operator

Consider a two-dimensional incompressible Newtonian fluid with spatially varying 
density $\rho(x,y,t)$ and dynamic viscosity $\mu(x,y,t)$. 
The viscous terms in the physical momentum equations are

$$ 
\begin{aligned} 
\mathcal V_x 
&=
\dfrac{1}{\rho} \left[ 
  \partial_x (2\mu \, \partial_x u) 
  + \partial_y \left(\mu(\partial_y u + \partial_x v)\right) 
\right],
\\ 
\mathcal V_y &=
\dfrac{1}{\rho} \left[ 
  \partial_x \left(\mu(\partial_y u + \partial_x v)\right) 
  + \partial_y (2\mu \, \partial_y v) 
\right]. 
\end{aligned} 
$$

Here $u$ and $v$ are the physical horizontal and vertical velocities. 
These expressions come from the divergence of the Newtonian viscous stress tensor,

$$ \boldsymbol{\tau} =\mu\left(\nabla\mathbf u+(\nabla\mathbf u)^T\right), $$

so, in particular, the shear stress is

$$ \tau_{xy}=\tau_{yx} =\mu(\partial_y u+\partial_x v). $$

Important: we retain the viscosity inside every derivative. 
Thus, the derivation includes spatial gradients of $\mu$; 
it does not assume constant viscosity.

For convenience, we define the dimensionless density and viscosity fields by

$$ \rho=\rho_1 R,\qquad \mu=\mu_1 M, $$

where $\rho_1$ and $\mu_1$ are reference values, and let

$$ \nu_1=\dfrac{\mu_1}{\rho_1}. $$

The task is to express $\mathcal V_x$ and $\mathcal V_y$ using the 
self-similar variables.

## 2. Write down the self-similar transformation

We use the transformation established for the planar Huppert heap:

$$ 
\begin{aligned} 
\tau 
&= 
\ln(t_{\rm phys}/T_0), 
\\ 
x &= L_0 \mathrm{e}^{\tau/5}\xi, 
& y &= H_0\mathrm{e}^{-\tau/5}\eta, 
\\ 
u 
&= \dfrac{L_0}{T_0} \mathrm{e}^{-4\tau/5}U, 
& v &= \dfrac{H_0}{T_0} \mathrm{e}^{-6\tau/5}V. 
\end{aligned} 
$$

Here $L_0,H_0,T_0$ are reference scales, while $U,V,\xi,\eta$ are the 
transformed variables.

At fixed physical time, the spatial derivatives transform as

$$ \boxed{ \partial_x=\frac{\mathrm{e}^{-\tau/5}}{L_0}\partial_\xi, \qquad \partial_y=\frac{\mathrm{e}^{\tau/5}}{H_0}\partial_\eta. } $$

These two relations are the heart of the calculation. Notice the opposite signs 
in the exponential factors: the $x$-scale expands like $\mathrm{e}^{\tau/5}$, 
whereas the $y$-scale contracts like $\mathrm{e}^{-\tau/5}$.

Substituting the velocity scalings gives four useful identities:

$$ 
\boxed{ 
  \begin{aligned} 
    \partial_x u 
    &=
    \dfrac{1}{T_0} \mathrm{e}^{-\tau} U_\xi, 
    \\ 
    \partial_y u 
    &= 
    \dfrac{L_0}{H_0T_0} \mathrm{e}^{-3\tau/5} U_\eta, 
    \\ 
    \partial_x v 
    &= 
    \dfrac{H_0}{L_0T_0} \mathrm{e}^{-7\tau/5} V_\xi, 
    \\ 
    \partial_y v 
    &= 
    \dfrac{1}{T_0} \mathrm{e}^{-\tau} V_\eta. 
  \end{aligned} 
} 
$$

Subscripts $\xi,\eta$ denote derivatives with respect to the transformed 
spatial coordinates.

For example,

$$ 
\begin{aligned} 
  \partial_y u 
  &= 
  \dfrac{\mathrm{e}^{\tau/5}}{H_0} 
  \partial_\eta \left( 
    \dfrac{L_0}{T_0} \mathrm{e}^{-4\tau/5} U 
  \right)
  \\ 
  &= 
  \dfrac{L_0}{H_0 T_0} \mathrm{e}^{-3\tau/5} U_\eta, 
\end{aligned} 
$$

because the spatial derivative is taken at fixed $\tau$.

Notice the different time exponents of $U_\eta$ and $V_\xi$. 
This is exactly why we must keep the two parts of the shear stress separate 
during the derivation.

## 3. Establish the acceleration scales

Before transforming the viscous terms in the momentum equations, 
we need the scale of each equation's acceleration.

Since $t_{\rm phys}=T_0\mathrm{e}^\tau$, a time derivative at fixed physical position 
scales as $1/t_{\rm phys}$. The transformed time derivative also contains 
the self-similar coordinate drift, but that does not change the overall 
acceleration scale.

For horizontal momentum, the characteristic acceleration is

$$ a_x=\dfrac{L_0}{T_0^2}\mathrm{e}^{-9\tau/5}. $$

For vertical momentum, it is

$$ a_y=\dfrac{H_0}{T_0^2}\mathrm{e}^{-11\tau/5}. $$

These follow from the velocity scales divided by $t_{\rm phys}$:

$$ \begin{aligned} \frac{u}{t_{\rm phys}} &\sim \frac{L_0}{T_0^2}\mathrm{e}^{-9\tau/5},\\\ \frac{v}{t_{\rm phys}} &\sim \frac{H_0}{T_0^2}\mathrm{e}^{-11\tau/5}. \end{aligned} $$

We will divide the physical viscous terms by $a_x$ in the horizontal equation 
and by $a_y$ in the vertical equation.

## 4. Derive the horizontal viscous operator, term by term

Recall

$$ \mathcal V_x =\dfrac{1}{\rho}\left[ \underbrace{\partial_x(2\mu u_x)}_{\text{normal stress}} +\underbrace{\partial_y(\mu u_y)}_{\text{shear contribution 1}} +\underbrace{\partial_y(\mu v_x)}_{\text{shear contribution 2}} \right]. $$

We now transform each of these three terms independently.

### 4.1 The normal-stress term: $\partial_x(2\mu u_x)$

From the derivative identity above,

$$ u_x=\dfrac{1}{T_0}\mathrm{e}^{-\tau}U_\xi. $$

Therefore,

$$ \begin{aligned} \partial_x(2\mu u_x) &=\frac{\mathrm{e}^{-\tau/5}}{L_0} \partial_\xi\left( 2\mu_1 M\frac{\mathrm{e}^{-\tau}}{T_0}U_\xi \right)\\ &=\frac{\mu_1}{L_0T_0}\mathrm{e}^{-6\tau/5} \partial_\xi(2MU_\xi). \end{aligned} $$

Dividing by the density and the horizontal acceleration scale gives

$$ 
\begin{aligned} 
  \dfrac{\partial_x(2\mu u_x)}{\rho \, a_x} 
  &=
  \dfrac{\mu_1}{\rho_1 R} \dfrac{\mathrm{e}^{-6\tau/5}}{L_0T_0} 
  \dfrac{T_0^2 \mathrm{e}^{9\tau/5}}{L_0} \partial_\xi (2 M U_\xi)
  \\ 
  &=
  \dfrac{1}{R} \underbrace{\left(\dfrac{\nu_1 T_0}{L_0^2}
  \mathrm{e}^{3\tau/5}\right)}_{A_x} \partial_\xi (2 M U_\xi). 
\end{aligned} $$

This identifies the first coefficient:

$$ \boxed{A_x=\frac{\nu_1T_0}{L_0^2}\mathrm{e}^{3\tau/5}.} $$

### 4.2 The first shear contribution: $\partial_y(\mu u_y)$

We have

$$ u_y=\dfrac{L_0}{H_0T_0}\mathrm{e}^{-3\tau/5}U_\eta. $$

Applying the $y$-derivative,

$$ \begin{aligned} \partial_y(\mu u_y) &=\frac{\mathrm{e}^{\tau/5}}{H_0} \partial_\eta\left( \mu_1M\frac{L_0}{H_0T_0} \mathrm{e}^{-3\tau/5}U_\eta \right) \\ &=\frac{\mu_1L_0}{H_0^2T_0} \mathrm{e}^{-2\tau/5}\partial_\eta(MU_\eta). \end{aligned} $$

Now we divide by $\rho a_x$:

$$ 
\begin{aligned} 
  \dfrac{\partial_y(\mu u_y)}{\rho \, a_x} 
  &= 
  \dfrac{1}{R} \dfrac{\mu_1}{\rho_1} \dfrac{L_0}{H_0^2 T_0} \dfrac{T_0^2}{L_0} 
  \mathrm{e}^{-2\tau/5+9\tau/5} \partial_\eta (M U_\eta)
  \\ 
  &=
  \dfrac{1}{R} 
  \underbrace{\left(\dfrac{\nu_1 T_0}{H_0^2} 
  \mathrm{e}^{7 \tau / 5}\right)}_{A_y} 
  \partial_\eta (M U_\eta). 
\end{aligned} 
$$

Thus,

$$ \boxed{A_y=\frac{\nu_1T_0}{H_0^2}\mathrm{e}^{7\tau/5}.} $$

The crucial observation is that this contribution differentiates in $y$ a 
stress built from $u_y$. Its geometric prefactor involves $H_0^{-2}$, 
not $L_0^{-2}$.

### 4.3 The second shear contribution: $\partial_y(\mu v_x)$

Now we consider the other part of the shear stress:

$$ v_x=\dfrac{H_0}{L_0T_0}\mathrm{e}^{-7\tau/5}V_\xi. $$

Consequently,

$$ \begin{aligned} \partial_y(\mu v_x) &=\frac{\mathrm{e}^{\tau/5}}{H_0} \partial_\eta\left( \mu_1M\frac{H_0}{L_0T_0} \mathrm{e}^{-7\tau/5}V_\xi \right) \\ &=\frac{\mu_1}{L_0T_0} \mathrm{e}^{-6\tau/5}\partial_\eta(MV_\xi). \end{aligned} $$

Dividing by $\rho a_x$,

$$ 
\begin{aligned} 
  \dfrac{\partial_y(\mu v_x)}{\rho \, a_x} 
  &= 
  \dfrac{1}{R} 
  \dfrac{\mu_1}{\rho_1} \dfrac{1}{L_0 T_0} 
  \dfrac{T_0^2}{L_0} 
  \mathrm{e}^{-6\tau/5+9\tau/5} \partial_\eta (M V_\xi)
  \\ 
  &= 
  \dfrac{1}{R} \underbrace{\left(\dfrac{\nu_1 T_0}{L_0^2} 
  \mathrm{e}^{3 \tau / 5}\right)}_{A_x} \partial_\eta (M V_\xi). 
\end{aligned} 
$$

Therefore this cross-derivative term has coefficient $A_x$, not $A_y$:

$$ \boxed{ \frac{\partial_y(\mu v_x)}{\rho \, a_x} =\frac{A_x}{R}\partial_\eta(MV_\xi). } $$

This is one of the most important checks in the derivation. 
Although this term is differentiated in $y$, the derivative acts on $v_x$, 
whose prefactor already contains $L_0^{-1}$. 
The factors $H_0$ cancel, leaving the coefficient based on $L_0^{-2}$.

### 4.4 Assemble the horizontal operator

Adding the three transformed contributions yields

$$ \boxed{ \mathcal D_U =\frac{1}{R}\left[ A_x\partial_\xi(2MU_\xi) +A_y\partial_\eta(MU_\eta) +A_x\partial_\eta(MV_\xi) \right]. } $$

This is the horizontal self-similar viscous operator.

## 5. Derive the vertical viscous operator, term by term

The physical vertical viscous term is

$$ \mathcal V_y =\dfrac{1}{\rho}\left[ \underbrace{\partial_x(\mu u_y)}_{\text{shear contribution 1}} +\underbrace{\partial_x(\mu v_x)}_{\text{shear contribution 2}} +\underbrace{\partial_y(2\mu v_y)}_{\text{normal stress}} \right]. $$

The important difference from the preceding calculation is that we now divide by 
the vertical acceleration scale,

$$ a_y=\dfrac{H_0}{T_0^2}\mathrm{e}^{-11\tau/5}. $$

### 5.1 The first shear contribution: $\partial_x(\mu u_y)$

Using

$$ u_y=\dfrac{L_0}{H_0T_0}\mathrm{e}^{-3\tau/5}U_\eta, $$

we obtain

$$ \begin{aligned} \partial_x(\mu u_y) &=\frac{\mathrm{e}^{-\tau/5}}{L_0} \partial_\xi\left( \mu_1M\frac{L_0}{H_0T_0} \mathrm{e}^{-3\tau/5}U_\eta \right) \\ &=\frac{\mu_1}{H_0T_0} \mathrm{e}^{-4\tau/5}\partial_\xi(MU_\eta). \end{aligned} $$

We divide by $\rho a_y$:

$$ 
\begin{aligned} 
  \dfrac{\partial_x(\mu u_y)}{\rho \, a_y} 
  &=
  \dfrac{1}{R} \dfrac{\mu_1}{\rho_1} \dfrac{1}{H_0 T_0} 
  \dfrac{T_0^2}{H_0} \mathrm{e}^{-4\tau/5+11\tau/5} 
  \partial_\xi (M U_\eta)
  \\ 
  &=
  \dfrac{A_y}{R}\partial_\xi (M U_\eta). 
\end{aligned} 
$$

Hence

$$ \boxed{ \frac{\partial_x(\mu u_y)}{\rho \, a_y} =\frac{A_y}{R}\partial_\xi(MU_\eta). } $$

Although the outer derivative is $\partial_x$, this term uses $A_y$. 
The inner gradient $u_y$ introduces $H_0^{-1}$, and after the outer 
derivative and normalization the coefficient is proportional to $H_0^{-2}$.

### 5.2 The second shear contribution: $\partial_x(\mu v_x)$

Recall that

$$ v_x=\dfrac{H_0}{L_0 T_0}\mathrm{e}^{-7\tau/5}V_\xi. $$

Thus,

$$ 
\begin{aligned} 
  \partial_x(\mu v_x) 
  &= \dfrac{\mathrm{e}^{-\tau/5}}{L_0} 
  \partial_\xi \left( \mu_1 M \dfrac{H_0}{L_0 T_0} 
  \mathrm{e}^{-7\tau/5} V_\xi \right)
  \\ 
  &=
  \dfrac{\mu_1 H_0}{L_0^2 T_0} 
  \mathrm{e}^{-8\tau/5} \partial_\xi (M V_\xi). 
\end{aligned} $$

Normalizing by $\rho a_y$,

$$ 
\begin{aligned} 
  \dfrac{\partial_x(\mu v_x)}{\rho \, a_y} 
  &=
  \dfrac{1}{R} 
  \dfrac{\mu_1}{\rho_1} \dfrac{H_0}{L_0^2 T_0} \dfrac{T_0^2}{H_0} 
  \mathrm{e}^{-8\tau/5+11\tau/5} 
  \partial_\xi (M V_\xi)
  \\ 
  &=
  \dfrac{A_x}{R} \partial_\xi (M V_\xi). 
\end{aligned} 
$$

Therefore,

$$ \boxed{ \frac{\partial_x(\mu v_x)}{\rho \, a_y} =\frac{A_x}{R}\partial_\xi(MV_\xi). } $$

Once again, this is the term with $v_x$ and therefore the coefficient 
based on $L_0^{-2}$.

### 5.3 The normal-stress contribution: $\partial_y(2\mu v_y)$

The transformed vertical normal gradient is

$$ v_y=\dfrac{1}{T_0}\mathrm{e}^{-\tau}V_\eta. $$

It follows that

$$ \begin{aligned} \partial_y(2\mu v_y) &=\frac{\mathrm{e}^{\tau/5}}{H_0} \partial_\eta\left( 2\mu_1M\frac{\mathrm{e}^{-\tau}}{T_0}V_\eta \right)\\ &=\frac{2\mu_1}{H_0T_0} \mathrm{e}^{-4\tau/5}\partial_\eta(MV_\eta). \end{aligned} $$

We divide by $\rho a_y$:

$$ 
\begin{aligned} 
  \dfrac{\partial_y(2\mu v_y)}{\rho \, a_y} 
  &=
  \dfrac{1}{R} \dfrac{\mu_1}{\rho_1} \dfrac{T_0^2}{H_0^2 T_0} 
  \mathrm{e}^{-4\tau/5+11\tau/5} 
  \partial_\eta (2 M V_\eta)
  \\ 
  &=
  \dfrac{A_y}{R} \partial_\eta (2 M V_\eta). 
\end{aligned} 
$$

Thus,

$$ \boxed{ \frac{\partial_y(2\mu v_y)}{\rho \, a_y} =\frac{A_y}{R}\partial_\eta(2MV_\eta). } $$

### 5.4 Assemble the vertical operator

Adding the three contributions gives

$$ \boxed{ \mathcal D_V =\frac{1}{R}\left[ A_x\partial_\xi(MV_\xi) +A_y\partial_\xi(MU_\eta) +A_y\partial_\eta(2MV_\eta) \right]. } $$

This is the second equation.

## 6. A compact check of all six terms

Here is the complete result in one place. 
The column headed “coefficient” shows the factor remaining after division by 
density and the appropriate acceleration scale.

| Momentum term before normalization | Resulting transformed term      | Coefficient |
| ---------------------------------- | ------------------------------- | ----------- |
| $\partial_x(2\mu u_x)$         | $\partial_\xi(2MU_\xi)$   | $A_x$   |
| $\partial_y(\mu u_y)$          | $\partial_\eta(MU_\eta)$  | $A_y$   |
| $\partial_y(\mu v_x)$          | $\partial_\eta(MV_\xi)$   | $A_x$   |
| $\partial_x(\mu v_x)$          | $\partial_\xi(MV_\xi)$    | $A_x$   |
| $\partial_x(\mu u_y)$          | $\partial_\xi(MU_\eta)$   | $A_y$   |
| $\partial_y(2\mu v_y)$         | $\partial_\eta(2MV_\eta)$ | $A_y$   |

The coefficients are

$$ 
\boxed{ 
  \begin{aligned} 
    A_x(\tau)
    &=
    \dfrac{\nu_1 T_0}{L_0^2} 
    \mathrm{e}^{3\tau/5},
    \\ 
    A_y(\tau)
    &=
    \dfrac{\nu_1 T_0}{H_0^2}
    \mathrm{e}^{7\tau/5}. 
  \end{aligned} 
} 
$$

The most useful pattern to remember is:

- Terms involving $v_x$ carry $A_x$.
- Terms involving $u_y$ carry $A_y$.
- The $u_x$ normal-stress term carries $A_x$.
- The $v_y$ normal-stress term carries $A_y$.

This pattern is specific to the chosen anisotropic scaling. 
It should not be inferred merely from whether the outer derivative is 
$\partial_x$ or $\partial_y$.

## 7. Why the final operators remain in divergence form

One last detail matters for checking the algebra. Since

$$ \mu=\mu_1M(\xi,\eta,\tau), $$

the transformed expressions retain $M$ inside each spatial derivative, 
for example

$$ \partial_\xi(MV_\xi), $$

rather than being simplified to $M V_{\xi\xi}$. 
If $M$ varies spatially, then

$$ \partial_\xi(MV_\xi) =M_\xi V_\xi+M V_{\xi\xi}. $$

Both contributions are part of the viscous stress divergence. 
The same reasoning applies to all six terms. 
The factor $R^{-1}$ sits outside the stress-divergence expressions because 
the original momentum equation has $\rho^{-1}\nabla\cdot\boldsymbol{\tau}$; 
we have not moved density inside those derivatives.

Final result: the operators follow directly from the physical Newtonian stress 
divergence, the anisotropic spatial transformation, and the horizontal and vertical 
acceleration scales. No additional assumption of constant viscosity is needed.