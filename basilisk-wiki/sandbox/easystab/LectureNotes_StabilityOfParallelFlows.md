# Linear stability analysis of parallel flows : General formalism

The objective of this section is to expose the general formalism used in the stability analysis of parallel flows, hence forming the starting points of chapters 6,7,8,9.

## Hypotheses

* The mean flow is parallel : 
 $\vec{u} = \bar{U}(y) \vec{e}_x$  for $y \in [y_1,y_2]$ ; $y_1$ and $y_2$ can be finite or infinite. 
 
 (examples : shear layers, wakes, jets, boundary layers...)
 
* Incompressible, constant density $\rho \equiv 1$.

* Viscosity $\nu \equiv Re^{-1}$ (using convenient nodimensionalization).

## Flow expansion

We expand the flow as base flow plus pertubations 


$$ 
{\left[
\begin{array}{c} u \\ v \\ w \\ p 
\end{array} 
\right]} \,= \, {\left[
\begin{array}{c} \bar{U}(y) \\ 0 \\ 0 \\ \bar{P} 
\end{array} 
\right]} \quad + \quad  {\left[
\begin{array}{c} u'(x,y,z,t) \\ v'(x,y,z,t) \\ w'(x,y,z,t) \\ p'(x,y,z,t) 
\end{array} 
\right]}
$$

Which we write symbolically $q = q_0 + q'$, introducing the the *state vector* notation $q = [u;v;w;p]$ 

Here $\bar{P}$ is either a constant (if gravity is neglected) or a hydrostatic field.

## Perturbation equations
 
 Injecting this ansatz in the Navier-Stokes equations and linearizing lead to the following set of equations:
 
$$ 
\left\{
\begin{array}{rcl} 
\partial_t u' + \bar{U} \partial_x u' + v' \partial_y \bar{U}
&=& - \partial_x p' +
Re^{-1} ( \partial_x^2 + \partial_y^2 + \partial_z^2) u'
\\
\partial_t v' + \bar{U} \partial_x v'  
&=& - \partial_y p' +
Re^{-1} ( \partial_x^2 + \partial_y^2 + \partial_z^2) v'
\\
\partial_t w' + \bar{U} \partial_x w'  
&=& - \partial_z p' +
Re^{-1} ( \partial_x^2 + \partial_y^2 + \partial_z^2) w'
\\
\partial_x u' + \partial_y v' + \partial_z w' &=& 0 
\end{array}
\right.
$$ 
 
## Modal decomposition in space
 
### Notations and physical meaning

 Owing to the invariance with respect to $x$ and $z$ directions, the perturbation can be considered in *modal form*. We introduce $k$ and $\beta$ the *axial and transverse wavenumbers*. On the other hand at this stage we keep the temporal dependance. 
 
Hence :
$$
 q'(x,y,z,t) = \Re [\tilde{q}(y,t) e^{i k x + i \beta z} ]
$$
 (In the sequel the symbol $\Re$ for the real part will be ommitted ; it is understood that only the real part of any complex expression is relevant).

Note that when both $k$ and $\beta$ are nonzero (but real) this modal expression describes an *oblique wave* with wavefronts oriented at an angle $\Theta = atan( \beta / k)$ with respect to the axial directions.

The case $\beta=0$ corresponds to *transverse waves* (chapters 6 and 7), and in this case the velocity component $w'$ can be dropped.

The case $k=0$ corresponds to *axially aligned structures* and is relevant to describe axial streaks and vortices (see chapter 8).

The case where $k$ is a *complex number* is relevant to describe spatial of spatiotemporal growth of perturbations and will be considered in chapter 9.

### Equations 
When injecting into the starting equations, all axial derivatires $\partial_x$ are replaced by $ik$ and transverse derivatives $\partial_z$ are replaced by $i\beta$.
The linear equations for the perturbation can be then written as a linear algebraic problem as follows:
$$
\frac{\partial}{\partial t}  {\mathcal B} \, \tilde{q} = {\mathcal A} \, \tilde{q}
$$

with 
$$
{\mathcal B} = 
\left[
\begin{array}{cccc} 
1 & 0 & 0 & 0 \\ 
0 & 1 & 0 & 0 \\
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 
\end{array} 
\right] 
$$

$$
{\mathcal A} = 
\left[
\begin{array}{cccc} 
-i k \bar{U} + Re^{-1} ( \partial_y^2 - k^2 - \beta^2) & - \partial_y \bar{U} & 0 & - i k \\ 
0 & -i k \bar{U} + Re^{-1} ( \partial_y^2 - k^2 - \beta^2) & 0 & - \partial_y \\
0 & 0 & -i k \bar{U} + Re^{-1} ( \partial_y^2 - k^2 - \beta^2) & - i \beta \\
i k  & \partial_y & i \beta & 0 
\end{array} 
\right] 
$$

### Modal decomposition in space and time

Eigenmode analysis consists of searching for modal perturbations in space *and time*, i.e:
$$
 q'(x,y,z,t) = \Re [\hat{q}(y) e^{i k x + i \beta z- i\omega t} ]
$$

This leads to the matricial eigenvalue problem:

$$
- i \omega  {\mathcal B} \, \hat{q} = {\mathcal A} \, \hat{q}
$$

Note that the notation $c = \omega/k$ (where $c$ is complex) is often used in the litterature, hence the eigenvalue problem can also be written $\left[ {\mathcal A} + i k c  {\mathcal B} \right] \, \hat{q} = 0.$

## Inviscid temporal stability analysis of parallel flows - Kelvin-Helmholtz instability (chap. 6)

This chapter has been moved [here](/sandbox/easystab/LectureNotes_Inviscid.md)


## Viscous temporal stability analysis of parallel flows - Tollmien-Schlishting instability (chap. 7)

This chapter has moved [here](/sandbox/easystab/LectureNotes_Viscous.md)

## 3D and transient growth mechanism in shear flows (chap. 8)

This chapter has moved [here](/sandbox/easystab/LectureNotes_NonModal.md)

## Spatial and spatio-temporal stability analysis (chap. 9)
  
  This chapter has moved  [here](/sandbox/easystab/LectureNotes_SpatioTemporal.md)
 
 [Back to main page](http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md)

