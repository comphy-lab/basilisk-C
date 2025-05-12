**(This document belongs to the lecture notes for the [M2-DET](http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md) course, D. Fabre, nov. 2018-dec. 2023)**



# Inviscid temporal stability analysis of parallel flows - Kelvin-Helmholtz instability (chap. 6)

In this chapter and the three following ones, we consider the instability properties of a *parallel* flow (i.e independent upon the spanwise direction $x$ and unbounded in this direction) defined as ${\bf u} = \bar{U}(y) {\bf e_x}$ for $y \in  [a,b]$ (the boundaries $a$, $b$ may be finite or infinite).

## Illustration
(Charru, sec. 4.1)

[https://www.youtube.com/watch?v=UbAfvcaYr00]()

[https://fr.wikipedia.org/wiki/Instabilité_de_Kelvin-Helmholtz]()

## Physical mechanism

The instability mechanism can be understood with a Bernoulli argument 
(cf. Charru,  sec. 4.3.1)

## Analytical solution for the zero-thickness shear layer
 
Consider the simplest case $\bar{U} = U_1$ (for $y<0$) and $\bar{U} = U_2$ (for $y<0$). See Exercice 6.0.

In this case the linear stability problem can be solved analytically. The solution is detailed [here](http://basilisk.fr/sandbox/easystab/Correction_Exercices.md#exercice-6.0)

This results in the *Dispersion relation* :
$$
(U_1-c)^2  + (U_2-c)^2  = 0
$$
whose solution can also be written as
$$
c = U_m \pm i \Delta U
$$
Where $U_m = (U_2+U_1)/2$ is the mean velocity and $\Delta U = |U_2-U_1|/2$ is half the velocity jump. One can thus deduce that:

*  There exists a unstable mode (as well as a associated complex-conjugate stable mode) for all $k$.
*  The phase velocity of the perturbation corresponds to the mean velocity, namely  $c_r =  U_m$ 
*  The growth rate $\omega_i$ of the unstable mode is proportional to the wavenumber and given by $\omega_i = k \Delta U$.

This last conclusion is somewhat unphysical, as small wavelengths (large $k$) perturbations can have an arbitrarily high amplification rate.
In practise, a "cutoff" is expected to be reached when the wavelength  $2 \pi k^{-1}$ becomes comparable to the thickness $\delta$. 
Apart from very special cases (such as the piecewise profile treated in Charru, 4.3.2) which have an analytical solution, the problem has to be solved numerically.

  

## Case of a continuous velocity profile

In the general case, one cannot make the hypothesis of a potential flow. One has to consider the linearized equations, either in primitive form (for $u,v,p)$ or in a reduced from called the Rayleigh equation (introducing the streamfunction $\psi$).

### Linearized equations

We will derive here linearized equations for 2D perturbations $(u',v',p')$ . Note that a more general derivation, considering 3D perturbations and nonzero viscosity, can be found in the document [General introduction to chap. 6-7-8-9](/sandbox/easystab/LectureNotes_StabilityOfParallelFlows.md) 

The expansion is:

$$ 
{\left[
\begin{array}{c} u \\ v  \\ p 
\end{array} 
\right]} \,= \, {\left[
\begin{array}{c} \bar{U}(y) \\ 0  \\ \bar{P} 
\end{array} 
\right]} \quad + \quad  {\left[
\begin{array}{c} u'(x,y,t) \\ v'(x,y,t)  \\ p'(x,y,t) 
\end{array} 
\right]}
$$

The linearized equations are:
$$
\left\{ \begin{array}{lcl}
\frac{\partial u'}{\partial t} + \bar{U} \frac{\partial u'}{\partial x} + v' \frac{\partial \bar{U}}{\partial y} &=&\frac{-1}{\rho}  \frac{\partial p'}{\partial x}, \\
\frac{\partial v'}{\partial t} + \bar{U} \frac{\partial v'}{\partial x}  &=&\frac{-1}{\rho}  \frac{\partial p'}{\partial y}, \\
\frac{\partial u'}{\partial x}  + \frac{\partial v'}{\partial y} &=& 0
\end{array} \right.
\qquad \mathrm{ ( Eqs. 1a, 1b, 1c)}
$$

### Eigenmode decomposition

In this chapter we consider perturbations under the form of eigenmodes with the following form:

$$[ u' ; v' ; p' ]  = [ \hat{u} ; \hat{v} ; \hat{p} ]  e^{ i k x} e^{- i \omega t} \equiv [ \hat{u} ; \hat{v} ; \hat{p} ]  e^{ i k (x - c t)}.$$ 



Here $c = \omega/k$ is an alternative notation for the eigenvalue, which will be useful in this chapter. We consider a *temporal formalism* in which $k$ is real and $\omega$ (and $c$) can be *complex*. The instability is again linked to the existence of a mode with positive growth rate $\omega_i$.

Note that $c_r = \omega_r / k$ can be interpreted as the *phase velocity* of the perturbation.

Assuming $\rho = 1$ for simplicity, the linear equations in primitive form are as follows:

$$
\begin{array}{lcl}
i k (\bar{U}-c) \hat{u} + \bar{U}'(y) \hat{v} &=& - i k \hat{p},  \\
i k (\bar{U}-c) \hat{v}   &=& - \partial_y  \hat{p}, \\
i k  \hat{u} + \partial_y\hat{v} &=& 0.
\end{array}
$$

This problem can also be written simply in matricial form as follows:
$$
\left[
\begin{array}{ccc} 
i k (\bar{U}-c) &  \partial_y \bar{U} &  i k \\ 
0 & i k (\bar{U}-c) &  \partial_y \\
i k  & \partial_y  & 0 
\end{array} 
\right] \hat{q} = 0
$$

## The Rayleigh equation

The divergence condition allows to introduce a streamfunction $\psi$. This streamfunction is again searched in eigenmode form:
$\psi(x,y,t) = \hat{\psi}(y) e^{i k (x -ct)}$.
The velocity components are related through :
$$ \hat{u} = \partial_y \hat{\psi} , \quad \hat{v} = - i k \hat{\psi}$$

The two first lines of the matricial problem can be combined to remove the pressure component, leading to a simple equation. Known as the *Rayleigh Equation*:


$$(\bar{U} - c) (\partial_y^2 - k^2) \hat \psi - \partial_y^2 \bar{U} \hat \psi = 0$$

### Mathematical analysis

<span style="color:green">
Along with suitable boundary conditions (either $\hat\psi \rightarrow 0$ as $|y| \rightarrow \infty$ for an unbounded domain, or $\hat{v}=-i k \hat\psi = 0$ at $y=a,b$ for a bounded domain), the Rayleigh equation is a continuous eigenvalue problem for the eigenvalue $c$ (or equivalently $\omega$). However, this equation does not belong to the class of problems for which theorems predicting the existence of solutions exist (such as Sturm-Liouville problems).
See [this document](http://basilisk.fr/sandbox/easystab/LinearSystems.md#linear-partial-difference-problems) for more details on these mathematical aspects.


<span style="color:green">
Experience shows the existence of two kind of solutions:

<span style="color:green">
* First, there can be *regular eigenmode* solutions for discrete values of $\omega$, but only a finite number (in most cases of interest, only 2). These eigenmodes always come in pairs or complex conjugate $\omega$, hence an unstable mode with $\omega_i>0$ is always associated with a stable one with $\omega_i<0$. 

<span style="color:green">
*(NB this property is characteristic of conservative problems where dissipative effects (such as viscosity) are neglected and is associated to time reversibility of the equations)*

<span style="color:green">
* Secondly, for continuous profiles, there exist *generalized eigenmode* solutions for values of $c = \omega/k$  in the *continuous interval*  $[min(\bar{U}), max(\bar{U})]$. Such generalized eigenmodes are *discontinuous* at the position $y_c$ such that $\bar{U}(y_c) = \omega/k$. (as mentionned in chap. 3 these modes cannot lead to instability).
 

### The Rayleigh and Fjortoft theorems

*Rayleigh Theorem:* 
The existence of an inflection point (i.e. a location $y_c$ such as $\bar{U}''(y_c)=0$ within the interval $y$ occupied by the flow is a *necessary* (but not sufficient) condition for instability.

*Demonstration:*
The demonstration consists of multiplying the Rayleigh equation by $\hat{\psi}^* / (\bar{U}(y)-c)$ where $^*$ means complex conjugate, and integrating over the interval $[a,b]$ filled by the flow.
The imaginary part leads to 

$$ 
c_i \int_{a}^{b} \frac{\partial_y^2 \bar{U}}{(\bar{U}-c_r)^2+c_i^2} |\hat{\psi}|^2 d y = 0.
$$  

If there exist an unstable mode ($c_i >0$), the integral must be zero, hence the numerator $\partial_y^2 \bar{U}$ must change sign along the interval.


*Fjortoft theorem :* 
This refinement of the Rayleigh theorem precises that the inflection point must correspond to a *maximum* of the absolute vorticity $|\partial_y \bar{U}|$

*Note:* 
Although it can be only demonstrated as a *necessary criterion*, in practise instability is always met when the Fjortoft criterion is verified. So it can actually be considered as a necessary and sufficient criterion.    




### Numerical solution for a continuous profile:

Let us consider the continuous profile corresponding to a $tanh$ profile:

$$ \hat{U}(y) = U_m + \Delta U \tanh( y/\delta)$$

Note that in the framework of temporal stability theory we can set $U_m=0$ without loss of generality because of Gallilean invariance. 

Numerical resolution of the Rayleigh equayion for this profile is done by the program [http://basilisk.fr/sandbox/easystab/KH_temporal_inviscid.m](KH_temporal_inviscid.m).

This program shows that :

* There exists an unstable mode for $k<k_c$ where $k_c \delta \approx 0.92$.
* This mode is most amplified for $k_c \delta \approx 0.5$ leading to amplification rate $\omega_i \approx 0.18$.
* For all values of $k$ there is a collection of spurious modes with eigenvalues filling the real interval $c \in [U_m-\Delta U_m+\Delta U]$. These modes represent the discretized version of the continous spectrum.
 
# Exercices
 
## **Exercice 0:** Kelvin-helmholtz instability of a zero-thickness shear layer

We consider a discontinous base-flow: $\bar{U}(u) = U_1$ for $y<0$ and  $\bar{U}(u) = U_2$ for $y>0$.

Show that the eigenvalues $c$ are given by 
$$
c = \frac{U_1+U_2}{2} \pm i \frac{U_1-U_2}{2}
$$

[Solution](http://basilisk.fr/sandbox/easystab/Correction_Exercices.md#exercice-6.0) 
 

## **Exercice 1:** Stability of a 2D swirling flow

[Solution](http://basilisk.fr/sandbox/easystab/Correction_Exercices.md#exercice-6.1) 

(cf. Drazin & Reid, section 3.15.3 & exercice 3.2)

We consider a swirling flow with mean velocity 
$\vec{u} = \bar{V}(r) \vec{e}_\varphi$ defined in an annular region $r \in [r_1,r_2]$.



Here $(r, \varphi,z)$ are cylindrical coordinates.



1. Starting from the Euler equations in cylindrical coordinates, write the primitive equations for a perturbation searched under eigenmode form as 
$q' = [u',v,'p'] = \hat{q} e^{i m \varphi - i \omega t}$.

2. Show that the equation can be reduced to a single equation which is the equivalent of the Rayleigh Equation for swirling flows:

$$
(m \Omega(r) - \omega) 
\left( \partial_r^2 + r^{-1} \partial_r -m^2/r^2 \right) 
\hat{\psi} - m/r \partial_r \Xi(r) \hat{\psi} = 0.
$$

Where $\Omega(r) = \bar{V}(r)/r$ is the angular velocity of the swirling flow and $\Xi(r) = \bar{V}(r) + \partial_r \bar{V}(r) \equiv 2 \Omega(r) + r \partial_r \Omega(r)$ is the corresponding vorticity.

Tip : you should introduce a streamfunction $\psi(r,\theta)$ such as $u = \frac{1}{r}\frac{\partial \psi}{\partial \theta}$ and $v = -\frac{\partial \psi}{\partial r}$.



3. Starting from this equation, prove the following criterion : *a necessary condition for instability is the existence of an extremum of the vorticity* $\Xi(r)$ *for some* $r \in [r_1,r_2]$.


## **Exercice 2:** Kelvin-Helmholtz instability of a piecewise-linear shear layer 

(Charru, section 4.3.2; Drazin & Reid, p. 146; Schmid & Henningson). 

[Solution](/sandbox/easystab/Correction_Exercices.md##exercice-6.2) 


Consider a base flow defined as follows :

$$
U(y) = \left\{
\begin{array}{ll} U_1 & \quad (y<-\delta) \\
U_m + \Delta U y/\delta & \quad (-\delta <y<\delta) \\
U_2 & (y>\delta) \\
\end{array}
\right.
$$
with $U_m = (U_1+U_2)/2$ and $\Delta U = (U_2-U_1)/2$.

1. Starting from the solution of the Rayleigh Equation in the three zones and using suitable matching conditions, show that modal perturbations are governed by the dispersion relation with the following form :
$$
4 k^2 \delta^2 (c-U_m)^2 - \Delta U^2 \left[ (2 k \delta -1)^2- e^{-4 k \delta} \right]= 0.
$$

2. Justify that the flow is unstable for $k\delta < k_c \delta$ with $k_c \approx  0.629 \delta^{-1}$ and that for $k \delta \ll 1$ the dispersion relation reduces to the classical one for a shear layer of zero thickness.

3. Plot the amplification rate $\omega_i = k c_i$ as function of $k$. Compare with the case of a contious profile (tanh profile).



## **Exercice 3:** Combined Rayleigh-Taylor-Kelvin-Helmoltz instability 

(Charru, exercice 4.5.1 ; exam jannuary 2018)

Consider two superposed fluid layers with different densities and velocities : 

$$
y<0 : \quad \rho = \rho_1 ; \, u = - U 
$$

$$
y>0 : \quad \rho = \rho_2 ; \, u = U
$$

Show that the dispersion relation governing modal perturbations is given by :
$$
\rho_1 (U+c)^2 + \rho_2 (U-c)^2 - \frac{(\rho_1-\rho_2)g}{k} - \gamma k = 0
$$

Show that this dispersion relation generalizes the cases previously considered.

Discuss the stability conditions as function of the parameters (see Charru, p. 132)

 [Back to main page](http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md)

