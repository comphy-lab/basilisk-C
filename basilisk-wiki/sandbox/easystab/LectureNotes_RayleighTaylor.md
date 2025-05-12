
**Rayleigh-Benard instability**

(Lecture notes for the M2R-DET course, D. Fabre, nov. 2018-janv. 2022)

This documents gives mathematical support for chapter 5 of the M2-DET course "introduction to hydrodynamical instabilities".

# Introduction

The situation considered in this chapter is the flow in a rectangular box of height $H$ and width $L$. The flow is driven by difference of temperature between the top and bottom plates: the bottom plate is maintained to temperature $T_0$, the top plate to temperature $T_1$ with $T_0>T_1$. 
Since the hot fluid is lighter, it wants to raise and the cold fluid want to fall. This is hindered by two diffusive effect: 
the fluid viscosity slows down the motion, and the thermal diffusivity will smear out the temperature.

The important parameters are the distance $H$ betwen the two plates, the temperature difference $T_0-T_1$ between the two plates, the mean fluid density $\rho_0$, the gravity $g$, the fluid dynamic viscosity $\mu$, the thermal diffusivity $\kappa$ (both assumed constant), and the thermal dilatation coefficient $\alpha$
(how much the fluid becomes lighter when it is heated). 


The state equation is a *dilatable, incompressible fluid* defined as follows:
$$
\rho(T)= \rho_0 [ 1 - \alpha (T-T_0) ].
$$


## Starting equations

We start with the Navier-Stokes equations for the fluid, coupled to the
heat equation for the temperature and the state equation for an
incompresssible, but dilatable fluid.
$$
\begin{array}{l}
\rho(u_t+uu_x+vu_y)=-p_x+\mu\Delta u\\
\rho(v_t+uv_x+vv_y)=-p_y-\rho g +\mu\Delta v\\
\rho_t + (u \rho_x + v \rho_y) + \rho (u_x+v_y) =0\\
T_t+u T_x+v T_y= \kappa \Delta T\\
\rho = \rho_0 \left[ 1  - \alpha (T-T_0) \right]
\end{array}
$$

We see that there is the volume force in the $v$ equation, dependent on the density, 
this will be the term through which the dilatation will affect the flow 
(the "buoyancy" or "floattability" term, telling how much the fluid "floats" when it is heated). 

# Linear analysis 

## Equations

### Base Flow

The "base flow" solution is the solution of the energy equation in the purely conductive regime ($\Delta T = \partial_{yy} T = 0$) :

$$ 
\overline{T}(y) = T_0 + \frac{(T_1-T_0) y}{H}
$$

The associated density field and pressure field (hydrostatic) are given by :

$$
\overline{\rho}(y)
= \rho_0 \left( 1 + \frac{\alpha (T_0-T_1) y}{H} \right)
$$

$$
\overline{P}(y) =  P_0 - \int \overline{\rho}(y) g dy
$$      

Here $\delta T = T_0 - T_1$ is the difference of temperature between the
lower and upper plates, and assumed positive, and $H$ is the height of the cell.



### Small-perturbation hypothesis and simplifications

We assume that the flow is a small-amplitude deviation to the "base state" previously described :

$$
\left[ \begin{array}{c} u \\ v \\ T \\ \rho \\ p  \end{array} \right]
= \left[ \begin{array}{c} 0 \\ 0 \\ \overline{T}(y) \\ \overline{\rho}(y) \\\overline{P}(y) \end{array} \right]  
+
\left[ \begin{array}{c} u(x,y,t) \\ v(x,y,t) \\ \theta(x,y,t) \\ \rho'(x,y,t) \\ p'(x,y,t) \\   \end{array} \right] 
$$

Under the hypothesis of small perturbations we can do the following simplifications :



- $(i)$ In the NS equations the nonlinear advection terms can be dropped, and $\rho$ can be replaced by $\rho_0$, except in the buoyancy term. 

- $(ii)$ The buoyancy and vertical pressure gradient terms lead to the Boussinesq approximation :

$$
- p_y - \rho g = - p'_y + \rho_0 g \alpha \theta
$$

- $(iii)$ The mass equation can be simplified to $div( u ) = 0$.

- $(iv)$ In the energy equation, the convective derivative is approximated a 
$$
\frac{d T}{dt} = \theta_t+ v \overline{T}_y 
$$



### Resulting set of linearized equations 
With these simplifications we end up with the following system of equations :

$$
\begin{array}{l}
\rho_0 u_t=-p'_x + \mu\Delta u\\
\rho_0 v_t=-p'_y + \mu\Delta v + \rho_0 \alpha g \theta \\
u_x+v_y=0\\
\theta_t = \frac{(T_0-T_1)}{H} v + \kappa \Delta \theta 
\end{array}
$$

## Analytical solution for the case of a vertical cell

If the cell is elongated in the vertical direction ($H \gg L$), it sounds reasonable to assume a monodimensional flow (invariance upon the $y$ direction, and $u \approx 0$) and only retain the vertical velocity $v$. Under the hypothesis of monodimensional flow the pressure can be dropped (same justification as for boundary layer theory).

Hence we are lead to:
$$
\begin{array}{l}
v_t= \nu v_{xx} +  \alpha g \theta \\
\theta_t = \frac{(T_0-T_1)}{H} v + \kappa \theta_{xx} 
\end{array}
$$
where $\nu = \mu/\rho_0$ is the kinematic viscosity.

Of course, the one-dimensional assumption is not justified in the top ($y \approx H$) and bottom ($y \approx 0$) parts of the cell where a flow reversal must occur. But it is a reasonable approximation in the central part of the cell ($0 \ll y \ll H$).

The $x$-dependency of $v$ and $\theta$ can be assumed sinusoidal. Hence :
$$
\begin{array}{l}
v(x,t) &=& V(t) \sin ( k x ) , \\
\theta(x,t) &=& \Theta(t)  \sin ( k x ).
\end{array}
$$

With the choice $k = 2 \pi / L$, this choice satisfies boundary conditions for $v$ at $x=0$ and $x=L$, and further corresponds to the assumed shape of the flow for a single convection cell ($v$ has opposite signs in the left and right halves of the domain).

Injecting the expression in the equations leads to a simple system which can be set in matricial form:
$$
\frac{d}{d t}
{\left[
\begin{array}{c} V \\ \Theta 
\end{array} 
\right]} 
\,= \, 
{\left[
\begin{array}{cc} - \nu k^2 & g \alpha \\ 
\frac{T_0-T_1}{H} & - \kappa k^2
\end{array} 
\right]}
\cdot {\left[
\begin{array}{c} V \\ \Theta 
\end{array} 
\right]} 
$$

Looking for eigenmodes solution with temporal dependency $e^\lambda t$, the
characteristic equation of this two-dimensional linear system is:
$$
det( {\mathcal A} - \lambda {\mathcal I})
= \lambda^2 + \underbrace{( \nu + \kappa) k^2}_{-tr(A)} \lambda 
+ \underbrace{\left[ \nu \kappa k^4 - \frac{T_0-T_1}{H} g \alpha \right]}_{det(A)} = 0
$$

Since the sum of eigenvalues $tr(A)$ is negative, the condition for instability is that the product of eigenvalues $det(A)<0$ (two real eigenvalues, one of them positive). This condition leads to:

$$T_0-T_1 > (\delta T)_c = \frac{ 16 \pi^4 \nu \kappa H}{\alpha g L^4}$$

Below this threshold, the problem is stable and no convection cells can sustain. Above the threshold, small-amplitude perturbations are amplified.
We can thus expect the appearance of a *bifurcation* at this value of $\delta T$.


## Case of a horizontal cell of large dimension

In the case of a horizontal cell of large dimensions ($L \gg H$) it is no longer possible to assume a one-dimensional flow, and both $u$ and $v$ (as well as $p$) have to be retained. On the other hand, we can assume a periodic behaviour in the $x$ direction with dependency $e^{i k z}$. Here $\Lambda = 2 \pi / k$ is the wavelength of the pattern which is still unknown (and not related to the dimension of the cell). 

Hence we assume : 

$$
\left[ \begin{array}{c} u(x,y,t) \\ v(x,y,t) \\ \theta(x,y,t) \\ p'(x,y,t) \\   \end{array} \right] 
=
\left[ \begin{array}{c} \hat{u}(y) \\ \hat{v}(y) \\ \hat{\theta}(y) \\ \hat{p}(y) \\   
\end{array} \right] e^{i k x } e^{\lambda t}
$$

Injecting in the equations leads to a generalized eigenvalue problem which can be set into the following block-matrix form:

$$
\lambda 
\underbrace{\left[ \begin{array}{cccc} 
\rho_0 & 0 & 0 & 0 \\ 
0 & \rho_0 & 0 & 0 \\ 
0 & 0 & 0 & 0 \\ 
0 & 0 & 0 & 1 
\end{array} \right]}_{\mathcal B}
\left[ \begin{array}{c} \hat{u}(y) \\ \hat{v}(y)  \\ \hat{p}(y) \\ \hat{\theta}(y) \\   
\end{array} \right]
=
\underbrace{\left[ \begin{array}{cccc} 
\mu (\partial_y^2 - k^2) & 0 & - ik & 0 \\ 
0 & \mu (\partial_y^2 - k^2) & - \partial_y & \alpha g \\ 
i k & \partial_y & 0 & 0 \\ 
0 & \frac{T_0-T_1}{H} & 0 & \kappa (\partial_y^2 - k^2)
\end{array} \right]}_{\mathcal A} 
\left[ \begin{array}{c} \hat{u}(y) \\ \hat{v}(y) \\ \hat{p}(y) \\ \hat{\theta}(y)  \\   
\end{array} \right]
$$



### Nondimensionalization
After convenient nondimensionalization, the problem can be reexpressed as:

$$
\lambda 
\underbrace{\left[ \begin{array}{cccc} 
1 & 0 & 0 & 0 \\ 
0 & 1 & 0 & 0 \\ 
0 & 0 & 0 & 0 \\ 
0 & 0 & 0 & 1 
\end{array} \right]}_{\mathcal B}
\left[ \begin{array}{c} \hat{u}(y) \\ \hat{v}(y)  \\ \hat{p}(y) \\ \hat{\theta}(y) \\   
\end{array} \right]
=
\underbrace{\left[ \begin{array}{cccc} 
Pr (\partial_y^2 - k^2) & 0 & - ik & 0 \\ 
0 & Pr (\partial_y^2 - k^2) & - \partial_y & Pr \\ 
i k & \partial_y & 0 & 0 \\ 
0 & Ra & 0 & (\partial_y^2 - k^2)
\end{array} \right]}_{\mathcal A} 
\left[ \begin{array}{c} \hat{u}(y) \\ \hat{v}(y) \\ \hat{p}(y) \\ \hat{\theta}(y)  \\   
\end{array} \right]
$$

where $Ra = \frac{\alpha g (T_0-T_1)H^3}{\nu \kappa}$ is the *Rayleigh number*
and $Pr = \nu / \kappa$ is the *Prandtl number*.

<span style="color:green"> 
Note on nondimensionalization : we have to write
$$ \hat{u} = U^* \tilde{\hat{u}} ; \quad 
\hat{v} = U^* \tilde{\hat{v}} ; \quad   
\hat{p} = P^* \tilde{\hat{p}} ; \quad 
\theta = \theta^*  \tilde{\hat{\theta}};
$$
$$
\lambda = {t^*}^{-1} \tilde{\lambda}; \quad 
y = H \tilde{y}; \quad k = H^{-1} \tilde{k}. 
$$
where $U^*,P^*$ and $\theta^*$ are the velocity, pressure and temperature scales, $t^*$ is the time scale.
The length scale is naturally taken as $H$. We chose the other scales as follows:
$$t^* = H^2/\kappa; \quad 
P^* = \rho_0 \kappa U^* /H; \quad 
\theta^* = (T_0-T_1) H U^* / \kappa
$$
We finally drop the tildes to get the final result.
</span>


### Numerical resolution for no-slip boundary conditions

Considering no-slip conditions at the top and bottom plates (i.e. $\hat{u} = \hat{v} = 0$ for $y=0$ and $y=H$), the problem cannot be solved analytically. On the other hand, the problem is easily solved numerically by using discretization in the $y$-direction, leading directly to a matrical problem with the form $\lambda B \hat{X} = A \hat{X}$ where $\hat{X}$ is the discretized version of 
$\hat{u}(y);\hat{v}(y);\hat{p}(y);\hat{\theta}(y)]$.

See [associated program.](http://basilisk.fr/sandbox/easystab/RayleighBenard.m)

The main conclusions are :

- For $Ra < Ra_c = 1708$, the problem is *stable* (i.e. all eigenvalues are negative) for all values of $k$.

- For $Ra > Ra_c = 1708$, the problem is *unstable* (i.e. one eigenvalue is positive) in a range or $k$ centered around $k_c \approx \pi$.

- Interestingly, the value of $Ra_c$ is independent upon the prandtl number $Pr$.

This means that above $Ra_c$, convection rolls are amplified, the most amplified ones corresponding to a dimensional wavelength $\Lambda = 2 \pi H / k \approx 2H$ (approximatively square convection patterns), whatever $Pr$. 


### analytical solution in case of slip boundary conditions 

The case of "slip" boundary conditions consists of replacing the boundary conditions by
$\hat{v} = 0 ; \partial_y \hat{u} =0$ at $y=0$ and $y=H$. This case is not physically justified, but allows to solve the problem in analytical way, and hence was favoured in the earliest studies. See Chandrasekhar or Drazin & Reid for details. This case is treated in [exercice 5.1](http://basilisk.fr/sandbox/easystab/LectureNotes_RayleighTaylor.md#exercice-5.1-rayleigh-benard-convection-with-free-free-boundaries.).


The resolution leads to the same conclusions as in the previous case, but the value of $Ra_c$ has an exact expression : $Ra_c= 27 \pi^4/4 \approx 657.5$. The  associated wavenumber is $k = \pi/\sqrt{2}$ which corresponds to a dimensional wavelength $\Lambda = \sqrt{2} H$, somewhat shorter than with no-slip conditions.





# Nonlinear analysis

To study what happens when convection starts, one has to consider nonlinear effects. In this section we consider a modelisation of the nonlinear effects called the *Lorenz model*. This model is not the most precise modelization of the actual dynamics of convection cells but has a strong historical significance, hence we restrict to this model in this lecture.

Physically, what happens in the nonlinear regime is that the term  $v \partial T/\partial y$, which is a source term in the temperature equation, is modified. 
The effect of the rolls is to "mix" the temperature, and as a result the temperature gradient $\partial T/\partial y$ is decreased compared to the value $(T_0-T_1)/H$ corresponding to the conduction state.

This effect is expected to lead to saturation of the amplitude of the rolls... but may lead to more exotic behaviour ! 


## Derivation of the Lorenz model

### Starting point
The starting point of the model is as follows :

$$ 
\begin{array}{l}
u(x,y,t) \approx X(t) \left( \frac{\pi}{k H}\right) \cos ( \pi y/H ) \cos( k x ),  
\\ 
v(x,y,t) \approx X(t) \sin ( \pi y/H ) \sin( k x ),
\\
\theta(x,y,t) \approx \bar{T}(y) + Y(t) \sin ( \pi y/H ) \sin( k x ) -Z(t) \sin( 2 \pi y /H).
\end{array}
$$

Explanations :

- The chosen expressions for $u$ and $v$ fit with the analytical solutions of the problem considering slip conditions (see exercice 1). The geometrical factor $\left( \frac{\pi}{k H}\right)$ is introduced so that the incompressibility condition is automatically satisfied.


- In the expression for $T$, the second term (of amplitude $Y$) corresponds to the temperature modification associated to the rolls, with the same modal dependance as $v$ (just as in the simplified case of the vertical cell considered above and in the solution with slip conditions treated in exercice 1).

- The additional term $Z(t)$ in the expression of $T$, which is independent of $x$, corresponds to the "flattening" of the temperature profile.
In effect, the averaged (in the x-direction) temperature profile is 
$$< T(y) > = \bar{T}(y) + Z(t) \sin( 2 \pi y /H)$$
Hence if $Z>0$ the temperature gradient is larger close to the boundaries ($y=0$ and $H$) and lower in the bulk, so the temperature profile is "flatter".

### Model
After a series of manipulations, the problem leads to 


$$
\frac{ dX}{dt} = Pr (Y-X)
$$
$$
\frac{ dY}{dt} = r X - Y - X Z;
$$
$$
\frac{ dZ}{dt} = X Y - bZ.
$$

Here $r$ is the *reduced Rayleigh number* : $r=Ra/Ra_c$ ; $Pr$ is again the Prandtl number, and $b$ is a geometrical factor linked to the aspect ratio of the convection cells (hence to the wavenumber $k$). The classical value initially chosen by Lorenz is $b=8/3$.

Considering the model, one can make the following remarks:

- When discarding the $Z$ term, the problem reduces to the same form as the simplified model considered above in the case of the verctical cell. One recovers the essential ingredients of the linear instability mechanism : the "motor terms" corresponding to buoyancy force and enhancing of temperature inhomogeneities by convection and the "opposing terms" corresponding to viscous and thermal dissipations.

- The second equation contains a new term. The physical interpretation is that the enhancing of temperature inhomogeneities by convection is reduced because the temperature gradient is reduced in the central zone of the cell. Hence the temperature gradient is represented by $(r-Z)$ instead of $r$.

- In the third equation, the source term $XY$ means that the "flattening" of the main temperature profile is increased if the vertical velocity $X$ and the temperature perturbations $Y$ have same sign (which is the case in the linear stages), but the effect may become inverse if $X$ and $Y$ happen to have opposite signs. 

- The last term $bZ$ corresponds to the simple effect of thermal diffusion which tends to damp the flattening of the main temperature profile $Z$.


## Dynamics of the Lorenz system.

Dynamics of the Lorenz system is best understood by playing with [the program](http://basilisk.fr/sandbox/easystab/lorenz_convection.m). On can also predict the first bifurcations by calculus (see exercice 5.2)

- The first system undergoes a *supercritical Pitchform bifurcation* for $r = r_{c,1} = 1$. Below this value, the only equilibrium solution is the trivial state $X=Y=Z=0$ which is stable. Above this value, this state becomes instable and new stable equilibrium solutions appear representing steady convection rolls.

- One can further prove that the next bifurcation is a *subcritical Hopf bifurcation* occurring for $r = r_{c,2} =  \frac{ Pr ^2 + b Pr + 3P}{Pr - b - 1} \approx 24.73$ (for $Pr= 10,b=8/3$). 


- Above $r_{c,2}$, the solution is *chaotic*. This word means that two initially arbitrarily close trajectories will always diverge, and after some finite time evolve in a complete different way. In other words, there is an extreme sensitivity to initial conditions, which means that the system, while perfectly deterministic, is actually unpredictible in the long-term.
In this range, trajectories in the phase space converge towards a *strange attractor* with a fractal stucture.

- (NB one can actually prove that the strange attractor arises for $r\gtrsim 24$. In the range $24<r<24.73$, depending on the initial conditions, the dynamics can converge either towards the fixed-point solution or towards the strange attractor).

- A number of other transitions occur at much larger values of $r$, for instance restabilizations, period-doubling events, etc... Try by yourself ! yu may for instance look at the range $r\in [140-170]$.

Note that in such ranges of parameters, the Lorenz model is no longer an accurate representation of what happens in the physical system of the convection cell. However, the Lorenz model played a huge large in the emergence of new concepts. In particular, the extreme sensitivity to initial condition of chaotic systems leading to long-term unpredictibility has became a paradigm of chaotic system, popularized by Lorenz himself as the "butterfly effect". Namely, when applying the idea to meteorology, which is iteslf a chaotic problem, we are led to the conclusion that the flapping of a butterfly's wing in Brasil can result, one month later, to a huge storm in Texas ! 


# Exercices for lecture 5


## Exercice 5.1 : Rayleigh-Benard convection with slip boundary conditions.

(sources: Drazin & Reid, sections 2.8 and 2.10 ; Rieutort ; Chandrasekar)

1. Starting from the stability equations in primitive form for $[\hat u, \hat v, \hat, p, \hat T]$, show that 
the problem can be reduced to a single differential equation of order 6 with the following form (Drazin & Reid, Eq. 8.39, p. 43) :
$$
(\partial_y^2 - k^2) (\partial_y^2 - k^2 - \lambda) (\partial_y^2 - k^2 - \lambda/Pr)\hat v + k^2 Ra \hat{v}= 0
$$ 
2. precise the set of boundary conditions to be used with this equation, in the respective cases of (a) rigid boundaries, and (b) "free-slip" boundaries.

3. In the case of "free-slip" boundaries, show that the least-stable eigenvalue is (D&R, eq. 10.5):
$$
\lambda = -\frac{1}{2}(1+Pr)(\pi^2 + k^2) + \frac{1}{2} \sqrt{ (Pr-1)^2(\pi^2+k^2)^2+  4 k^2 \, Ra \,  Pr/(\pi^2+k^2)}
$$  


4. Show that for a given $k$, the flow is unstable for $Ra >Ra_j(k)$, with 

$$
Ra_j(k) = (\pi^2+k^2)^3/k^2
$$ 

(D&R, Eq. 10.6)

Draw this curve in the $(Ra,k)$ plane.

5. Show that the instability threshold is $Ra_c = min [Ra_j(k)] =  27 \pi^4/4 \approx 657.5$ and give the corresponding value of $k$.

## Exercice 5.2 Study of the Lorenz system.

(Ref : Charru) 

1. Find all fixed-point solutions of the Lorenz system. Show that the system undergoes a supercritical Pitchfork bifurcation for $r>r_{c,1}=1$.

2. Show that if $Pr>b+1$, the stable solution existing for $r>1$ becomes unstable through a Hopf bifurcation 
for $r>r_{c,2}$ with 
$$
r_{c,2} = \frac{ Pr ^2 + b Pr + 3Pr}{Pr - b - 1}.
$$ 
Give the numerical value for $b=1$ (Charru) and for $b=8/3$ (initial choice of Lorenz).


Tip : Study the stability of the fixed-point solution. Write the characteristic polynomial governing its stability under the form $P(\lambda) =  \lambda^3 + A \lambda^2 +B \lambda + C$. Remarking that if a Hopf bifurcation occurs, this polynomial has two purely imaginary roots and one real root, write a condition on the coefficients $A$,$B$,$C$ of the polynomial.

3. Using the program, check that the bifurcation is subcritical.

4. Using the program, verify the existence of a "strange attractor" for $r>r_s$ with $r_s\approx 24$ (try to find the best approximation of $r_s$ !) 

5. Verify that for $r_c<r<r_s$ both stationnary and "strange" solutions are possible depending on initial conditions.

6. Observe the behaviour of the Lorenz system for much larger values or $r$. Observe in particular the behaviour in the range $r \in [145, 170]$. 

 

 



