

**Spatial and spatio-temporal stability analysis (chap. 9)**
   

# Introduction

   Up to now, we have considered a *temporal stability framework* which consists of considering modal perturbation proportional to 
$e^{i k x}$ and considering their possible growth/decay in time. 
Such perturbations represents wave-like perturbations occupying the whole domain $x \in [ -\infty,\infty]$ and it is somewhat unphysical to consider a single wave of this kind as an initial condition. 
    
   In this chapter we will modify this starting point in order to address two related questions :
    
- what happens if a shear flow is perturbed by an initial condition *localised in space* ? (**Initial value problem**) 
- what happens if the shear flow is *continuously forced by spatially localized activator*  immposing a sinusoidal forcing at a (real) angular frequency $\omega$ ? (**Signaling problem**) 

**Note :** As in chapters 6-7 we restrict to 2D perturbations ($\beta = 0$) and still  assume the base flow to be parallel.

# The signaling problem.


## Starting point
We will first address the second question. Considering that the forcing (applied at $x=0$) is periodic in time $e^{-\omega t}$, it seems reasonable to assume that the perturbation $q'$ of the flow will also be periodic in time. On the other hand, the perturbation can be expected to grow in space. Thus a convenient starting point is :

$$ q'  = \hat{q} e^{i k x} e^{- i \omega t} ; \quad \omega \in \mathbb{R} ; k \in \mathbb{C} 
$$


Or, written differently : $q'  = \hat{q} e^{i (k_r x - \omega t)} 
e^{-k_i x}$.

 Remarks :
 
- In this formalism $-k_i$ represents the *spatial growth* rate in the direction $x \rightarrow + \infty$. A flow is thus *Spatially unstable* if there exists some real $\omega$ such that $k_i <0$. 
- In addition, $c_r = \omega / k_r$ represents the *phase velocity* of the disturbance.
- Note that the condition $||q'||<<1$ necessary for linear stability analysis cannot be expected to hold for all $x$ as soon as the flow is spatially unstable. Necessarily, nonlinearities are dominant after a finite distance $x$. It is however still possible to use the linear equations locally, in some finite interval.

## Mathematical analysis
 
In the spatial framework, the linearized equations 
 (either in primitive form for $\hat{ \bf q}$, or in reduced form (Rayleigh or Orr-Sommerfeld)  for $\hat \psi$) still constitute an eigenvalue problem,
  but the resolution now consists of solving for the eigenvalue $k$ as function of the (real) frequency $\omega$ and the base-flow properties.
 
For temporal stability results, the growth rate does not depend upon the mean velocity $U_m$ because of gallilean invariance. This is not the case for spatial stability results.
 
### Example

For instance, for the piecewise shear layer, the solutions $\omega/k = U_m \pm i \Delta U$ can be inverted to 

$$
k = 2 \omega \frac{ U_m \mp i \Delta U}{ U_m^2 + \Delta U^2}
$$ 
 
## Numerical resolution 

When considering $k$ as an eigenvalue, the linear stability equations (in primitive form) can be set under the matrical form:
$$
k {\mathcal B} \hat{q} = {\mathcal A} \hat{q}
$$

In the inviscid case, the equations can be directly reduced to this form with 
$\hat q = [ \hat{u}; \hat{v}; \hat{p}]$. It is then straightforward to build a numerical program to resolve this eigenvalue problem. 

See the program [KH_spatial_inviscid.m](/sandbox/easystab/KH_spatial_inviscid.m) doing the resolution for the case of the tanh shear layer.

In the viscous case, as the wavenumber $k$ appears quadratically in the equations, it is not possible to write directly a eigenvalue problem for $\hat {\bf q} = [ \hat{u}; \hat{v}; \hat{p}]$. On the other hand, this can be done for $\hat q = [ \hat{u}; \hat{v}; \hat{p} ; \hat{u}_1 ; \hat{v}_1]$ 
where $\hat{u}_1$ and $\hat{v}_1$ are auxiliary unknown fields defined by
$\hat{u}_1 = k \hat{u}$;  $\hat{v}_1 = k \hat{v}$.

See the program [KH_spatial_viscous.m](/sandbox/easystab/KH_spatial_viscous.m) to see how this this method can be implemented to solve the spatial stability problem in the viscous case.

## Typical results for shear flows

We consider as a prototype shear flow the tanh shear layer defined as 
$$\bar{U}(y) = U_m + \Delta U \tanh(y)$$
We note $R = U_m / \Delta U$. If $0<R<1$ the velocity is always in the positive $x$-direction , while if $R>1$ the lower layer goes in the negative direction.

Results obtained from the numerical resolution of the eigenvalue problem ([KH_spatial_inviscid.m](/sandbox/easystab/KH_spatial_inviscid.m)), or from model equations ([Such as the one considered in thie program](/sandbox/easystab/SpatioTemporal_ModelEquations.m))
observe that : 

1. If $R$ is large, one obtains an eigenvalue noted $k = k^+ (\omega,R)$, 
which is unstable in a limited range of $\omega$ (say $Im(k^+)<0$ for 
$\omega_1<\omega<\omega_2$) and is associated to an eigenmode qualitatively 
similar to the one predicted by the temporal theory. This mode corresponds to the same instability as in temporal theory ((Kelvin-helmholtz instability, here), except that it grows with space.  
 

<span style="color:green">
NB in the limit $R \gg 1$ (or $U_m \gg \Delta U$) one can demonstrate that the spatial eigenvalie $k = k^+ (\omega)$ can be deduced from the temporal growth rate $\omega = \Omega(k)$ using the Gaster transformation :

$$
k^+ (\omega) =  \omega / U_m  - \Delta U / U_m \Omega( \omega / U_m)
$$ 
 
 </span>
  
2. As $R$ increases, the amplification rate of the $k^+$ solution increases ; moreover one observes the appearance of one (or possibly several) other branch(es) called $k^-$ for which $Im(k^-)<0$ for all values of $\omega$.

The fact that $-Im(k^-)>0$ means that the amplitude of this mode grows as $x$ increases, which seems to predict an other spatially developing instability.  
However, thanks to a mathematical analysis (explained in next section), it can be shown that this second kind of solutions actually corresponds to perturbations *propagating in the upstream direction*   which are *spatially damped* as $x \rightarrow - \infty$.



3. As $R$ increases further, one reaches a particular value $R=R_c$ for which the spatial branches coincide ; i.e.  $k^+(\omega_0) = k^-(\omega_0)$ for some *real* value $\omega_0$ of $\omega$.


4. For $R<R_c$, one still oberves two branches of solution $k(\omega)$ of the spatial problem but they exchange their identity. As soon as it occurs it is no longer possible to identify them as $k^+$ and $k^-$ and associate them to the developement of perturbations propagating upstream and downstream, respectively.
 
NB : the actual value of $R_c$ is 1.319 for the tanh shear layer (Huerre & Monkewitz, 1988), and $R_c = 1$ for the simple model considered in exercice 1. 

# The initial value problem.

## Examining the failure of the spatial problem

We have seen that for $R<R_c$ the spatial "signaling" problem admits a well-defined solution both for $x>0$ (corresponding to $k^+$ solutions) and $x<0$ (corresponding to $k^-$ solutions), while for $R>R_c$ it leads to results for which we are not able to give a physical interpretation. What happens in such cases ???

To understand this, reconsider the starting point of the spatial stability analysis. We have supposed that both the localised forcing at $x=0$ and the perturbation of flow are proportional to $e^{i \omega t}$. This corresponds to the  *permanent response*, namely the response to a forcing exists for all instants starting from $t=- \infty$. In reality, the forcing starts at a given instant (say $t=0$). Then the perturbation is the sum of the *permanent response* and a *transient response*. Thus the spatial stability results are only relevant to the cases where the transient regime decays.

## Mathematical formulation of the initial value problem

Let us go back to the initial value problem. If the initial condition is square-integrable, we can always perform a Fourier transform in $x$ and write the solution as follows:
$$
 q'(x,y,t) = \int_{-\infty}^{+\infty}  \tilde{q}(k,y,t) e^{i k x } \, dk
$$

for each value of $k$, the evolution is given by
$$
\frac{\partial}{\partial t}  {\mathcal B} \, \tilde{q} = {\mathcal A}_k \, \tilde{q}
$$

and the we have seen that the solution can be projected on the basis of eigenmodes (including discrete and continuous ones), namely :
$$
\tilde{q}(k,y,t) = \sum_{\omega_n \in Sp^d(k)} c_n(k) \, \hat{q}_n(k,y) \, e^{-i \omega_n(k) t} + \int_{Sp^c} c(k,\omega)  \hat{q}_\omega(k,y) e^{-i \omega t}
$$ 
where $c_n(k)$ and $c(k,\omega)$ can be obtained from the projection of the initial condition $q'(x,y,t=0)$.



To simplify, we keep only the leading eigenmode $(n=1)$ of the discrete spectrum   and discard the continuous spectrum, namely 
$\tilde{q}(k,y) \equiv c(k) \hat{q}(k,y) e^{-i \omega(k) t}$

Reintroducing this solution in the Fourrier transform leads to

$$
 q'(x,y,t) = \int_{-\infty}^{+\infty}  c(k) \hat{q}(k,y) e^{i k x} e^{-i \omega(k) t} \, dk
$$

where $c(k)$ can be obtained from the $x$-Fourier transform of the initial condition $q'(x,y,t=0)$.

## Asymptotic behavior

The asymptotic behaviour of the previous expression can be predicted using the *stationnary phase theorem* 

(Applicable if $Im ( \omega(k) ) \leq 0$ as $k \rightarrow \pm \infty$ ($k \in \mathbb R$), which physically means that the short-wave perturbations are not amplified).

**Theorem (stationnary phase)**

1. If there exists a (complex) $k_0$ such that $\partial \omega/\partial k = 0$ (corresponding to a *saddle point* of the complex function $\omega(k)$), and if there exists a *steepest descent path* in the complex $k$-plane linking $-\infty$ to $+\infty$ and passing through this saddle point; then
$$
 q'(x,y,t) = \sqrt{\frac{ 2 \pi}{ i \omega''(k_0) t} } 
   c(k_0)\hat{q}(k_0,y)  e^{ i k_0 x} e^{ - i \omega(k_0) t}
  \mathrm{ as } \quad t \rightarrow +\infty \quad (\mathrm{ with \, fixed \,  } x) 
$$
2. Otherwise
$$
 q'(x,y,t) = o(t^{-1/2}) \qquad \mathbf{ as } \quad t \rightarrow +\infty \quad (\mathrm{ with \, fixed \,  } x) 
$$

The demonstration of the theorem relies on the fact that the argument in the integral expession of $q'(x,y,t)$ is an analytical fonction of the complex variable $k$. Hence, a theorem of complex anaysis states that the contour of integration from $k=-\infty$ to $k=+\infty$ can be arbitrarily deformed (provided it does not cross singularities of the function). Taking the aformentionned steepest descent path contour leads to the result.

## Interpretation : convective/absolute instabilities

This results leads to the following conclusion :

1. If the complex frequency $\omega_0 = \omega(k_0)$ corresponding to the saddle point verifies $Im( \omega_0 <0)$ (or if there is no such saddle point), then the perturbation will decay with time at every fixed position $x$. The instability is said to be **convective**. In this case the transient response decays and it is allowed to use spatial stability analysis to study the permanent response to a harmonic forcing.

2. If the complex frequency $\omega_0 = \omega(k_0)$ corresponding to the saddle point verifies $Im( \omega_0 >0)$ , then the perturbation will grow with time at every fixed position $x$. The instability is said to be **absolute**. In this case the transient response always dominate and the spatial stability is not relevant.


## Asymptotic behaviour in a moving frame


The notion of a *convectively unstable* flow may sound paradoxal: this case implies that there exists a (real) value of $k$ such as perturbation proportional to $e^{ikx}$ is exponentially amplified;  BUT on the other hand at a given $x$ the response to ANY spatially localized perturbation decays.
 
To clarify this point let us reconduct the stationary phase argument in a frame of reference moving at a velocity $V$. Writing $x=Vt+x'$ allows to re-write the solution as 
$$
 q'(x,y,t) = \int_{-\infty}^{+\infty}  c(k) \hat{q}(k,y) e^{i k x'} e^{-i (\omega(k)-kV) t} \, dk
$$

Letting $t \rightarrow \infty$ BUT keeping $x'$ fixed and applying again the stationary phase therorem leads to:

$$
 q'(x,y,t) = \sqrt{\frac{ 2 \pi}{ i \omega''(k_V) t} } 
   c(k_V)\hat{q}(k_V,y) e^{i k_v x'} e^{  - i\omega_V  t}
  \mathrm{ as } \quad t \rightarrow +\infty \quad (\mathrm{ with \, fixed \,  } x'=x-Vt) 
$$

where $k_V$ is the solution of $\partial \omega / \partial k = V$ and $\omega_V = \omega(k_V) - k_V V$ the corresponding frequency (in the moving frame). If $Im(\omega_V) >0$, then the perturbation grows as $x,t \rightarrow \infty$ along the *ray* defined by $x/t = V$.



We can then define the interval $[V_1,V_2]$ corresponding to velocities of *rays* along which the perturbation grow. The rays with velocities $V_1$ and $V_2$ can be interpreted as the upstream an downstream *fronts* of the *unstable wavepacket* (region in the $x-t$ plane where the perturbation will grow).
 

1. For an *absolutely unstable* case, the interval $[V_1,V_2]$ comprises $V=0$, so an initially localized perturbation will tend to occupy a region in the $x-t$ plane  expanding both upstream ($V_1<0$) and dowstream $(V_2>0$). 

2. On the other hand, in a *convectively unstable* case, an initially localized perturbation will give rise to an *unstable wavepacket* which will be advected in the positive direction. At any fixed $x$, the perturbation will become large as the unstable wavepacket crosses the $x$ location, but will eventually decay.


##  Reexamining the spatial approach

We have seen that there exists a value $R=R_c$ at which there is a coincidence of the spatial branches ; namely there exists two solutions $k^+(\omega), k^-(\omega)$ such as $k^+(\omega_0) = k^-(\omega_0)$ for some *real* $\omega_0$.
This coincidence of two branches precisely corresponds to a *double root* of the function $\omega(k)$. The existence of a double root precisely corresponds to the condition $\partial \omega / \partial k =0$.

In complex analysis, a coincidence of two branches of solutions when a parameter $R$ crosses a value $R_c$ is always associated to an exchance of identities of the two branches, which reconnect differently for $R<R_c$ and $R>R_c$, respectively. We thus see that the transition between the situation with distinct $k^+$ and $k^-$ branches for $R<R_c$ and the situation where the branches is completely explained by the existence of a saddle point along the real $\omega$-axis for $R=R_c$.
  
# Exercices for lecture 9

### **Exercice 1** Absolute/convective transition for a model equation.

We consider a simple model characterised by the following dispersion relation :

$$ D(\omega,k) = (\omega - k U_m) -  i \Delta U k (1-k) = 0$$

1. Considering *temporal stability analysis*, compute $\omega$ (complex) as function of $k$ (real). Plot the amplification rate $\omega_i$ as function of $k/k_c$. Justify that the present model can be used as an aproximation for the stability properties of a finite-thickness shear layer.

2. Find the value $k_0$ corresponding to the saddle point of the relation and deduce the corresponding absolute frequency $\omega_0$. Show that 
$$
\omega_0 = \frac{\Delta U}{2}  + \frac{ i ( (\Delta U)^2-U_m^2) }{4 \Delta U}    $$

3. Deduce ranges of $R= \Delta U/U_m$ corresponding respectively to absolute and convective instabilities.
(it is admitted that the saddle point verifies the "steepest descent path" property) 

4. Compute the spatial stability branches $k(\omega)$ of this model and plot them in the convective and absolute cases (you may start from [this program]( 
http://basilisk.fr/sandbox/easystab/SpatioTemporal_ModelEquations.m) and adapt it).


### **Exercice 2** Numerical resolution of the spatial stability problem.

We consider modal perturbations of a parallel shear flow under the form :

$$ 
[u,v,p] = [ \bar{U}(y) ,0, 0 ] + [ \hat{u}, \hat{v}, \hat{p}] e^{i k x - i \omega t}
$$ 

1. Considering non-viscous theory, show that the linear stability problem can be set under the form of a generalized eigenvalue problem:

$$
k {\mathcal B} \hat{q} = {\mathcal A} \hat{q} \qquad \mathrm{with} \quad \hat q = [ \hat{u}; \hat{v}; \hat{p}]
$$
and give the expression of the matricial operators $\mathcal A$ and $\mathcal B$.

(Solution : see the program [KH_spatial_inviscid.m](/sandbox/easystab/KH_spatial_inviscid.m))

2. In the viscous case, show that the problem can be reduced to 

$$
k {\mathcal B}^v \hat{q}^v = {\mathcal A}^v \hat{q}^v \qquad \mathrm{with} \quad   \hat q = [ \hat{u}; \hat{v}; \hat{p} ; \hat{u}_1 ; \hat{v}_1] \qquad \mathrm{ ; } \quad \hat{u}_1 = k \hat{u}; \quad  \hat{v}_1 = k \hat{v}
$$

Give the expression of the matricial operators.

(Solution : see the program [KH_spatial_viscous.m](/sandbox/easystab/KH_spatial_viscous.m))




# Recommended lectures

For a more complete discussion of the spatio-temporal aspects of instabilities, the recommended reference is the lecture of Huerre & Rossi (1998), which is actually a chapter of a book edited by Godrèche \& Manneville (1998).
  
