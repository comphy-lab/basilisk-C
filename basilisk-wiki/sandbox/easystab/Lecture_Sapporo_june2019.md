

**Numerical methods for flow instability problems **

**D. Fabre, IMFT (Toulouse, France) **

**Lecture in Sapporo University, june 27, 2019 **

# Linear dynamical systems and linear PDE problems

## Linear dynamical systems

(extract from [david/LinearSystems.md]() )

- A linear dynamical system of dimension $N$ is as follows :

$$
\frac{d X}{dt} = A X
$$

Where $X(t) = [x_1(t);x_2(t);.. ;x_N(t)]$ is a vector (column-vector) of dimension $N$ and $A$ is a square matrix of dimension $N$.

- Eigenmode solutions are sought under the form 
$$
X(t) = \hat{X} e^{\lambda t}
$$

- The dynamical system thus results to a eigenvalue problem :

$$
\lambda \hat{X} = A \hat{X}
$$

The solutions $\lambda$ are the *eigenvalues* and the corresponding 
$\hat{X}$ the *eigenvectors*. The full set of eigenvalues ${\mathcal Sp} =
\{\lambda_n, n=1..N \}$ is called the *spectrum* of the matrix.

- The general solution of the dynamical system is (except in exceptional degenerated cases)


$$ 
X(t) = \sum_{\lambda_n \in {\mathcal Sp}} c_n \hat X_n e^{\lambda_n t} \qquad \mathrm{ (Eq. 1)}
$$

- The system is *linearly unstable* if there exists an eigenvalue verifying $Re (\lambda) > 0$



## Linear algebraic systems

- A *linear algebraic system* is defined as follows :

$$
B \frac{d X}{dt} = A X
$$

where $A$ and $B$ are square matrices of dimension $N$.

- If $B$ is inversible, this is equivalent to a dynamical system with 
$\frac{d X}{dt} = ( B^{-1} A ) X$

- If $B$ is non-inversible, one can still look for eigenvalue/eigenvector solutions by solving the *generalized eigenvalue problem*

$$
\lambda B \hat{X} = A \hat{X}
$$

There are again $N$ eigenpairs, but **some of the eigenvalues are infinite**

- The general solution can still be written as


$$ 
X(t) = \sum_{\lambda_n \in {\mathcal Sp}^*} c_n \hat X_n e^{\lambda_n t} 
$$

Where ${\mathcal Sp}^*$ is the set of *finite eigenvalues*

## Numerical resolution of eigenvalue (and generalized eigenvalue) problems

**Through the characteristic polynomial

Mathematically elegant but practically efficientless for $N>3$ !

**Direct method : matrix diagonalization.**

$$
B^{-1} A = U S U^{-1}
$$
where $S = diag(\lambda_n)$ is diagonal, and $U =  [\hat{X}_1, \hat{X}_2, ..., \hat{X}_N]$ contains the eigenvectors sorted in columns. 

-> Costly, and generally useless to compute all eigenvalues.

- Fortran/c++ : **Lapack** ;  Matlab/Octave/Python : **eig**

**Iterative methods **

- Power method 

The iteration scheme 
$$X^{n} = \frac{A X^{n-1}}{|| A X^{n-1} ||} $$ 
asymptotes to $X^{n} \approx \hat{X}_{\lambda_max}$

-> allows to compute a *single* eigenvalue/eigenvector pair $\lambda_{max}$ (the one with maximum absolute value)

- "Shifted inverse power method" 
The iteration scheme 
$$X^{n} = \frac{B ( A - \sigma B)^{-1} X^{n-1}}{||B ( A - \sigma B)^{-1} X^{n-1}||}$$ 
asymptotes to 
$X^{n} \approx \hat{X}_{\lambda_0}$

-> allows to compute the eigenvalue $\lambda_0$ closest to a specified *shift*
$\sigma$.

- Arnoldi "Shift-invert" method

Variant allowing to compute the $n$ eigenvalues closest to a given *shift* $\sigma$

-> This is the prefered method in Stability analyses

- Fortran/c++ : **Arpack** ; Matlab/Octave/Python : **eigs**

- For large systems handled with parallelization: **SLEPc** (no simple interface with Matlab/Octave/Python)

## Linear PDE problems

- Consider a linear partial difference equation (PDE) as follows :

$$
\frac{\partial F(x,t)}{ \partial t} = {\mathcal A} F(x,t)
$$

where $F(x,t)$ is defined on a continuous interval ($x\in [a,b] ; t \in [0, \infty]$  and $\mathcal A$ is a *linear operator* involving $x$-derivatives ( operators $\partial / \partial x$ , $\partial^2 / \partial x^2$, ... )

- Define a *mesh* $x = [x_1, x_2, ... x_N]$ spanning the interval $[a,b]$

- Define the discretized version of $F$ as 

$$
F = \left[ \begin{array}{c} 
F_1 \\  F_2 \\ ... \\ F_N 
\end{array}
\right] = \left[ \begin{array}{c} 
F(x_1) \\  F(x_2) \\ ... \\ F(x_N) 
\end{array}
\right]
$$ 

- Then the differential operators $\partial / \partial x$ and 
$\partial^2 / \partial x^2$ are approximated by **differentiation matrices**
noted *dx*, *dxx*

- These matrices are used to build the discretized version $A$ of the operator $\mathcal A$

## Discretization methods

- Finite difference method

-> mesh $x = [x_1, x_2, ... x_N]$ is *uniform* : $x_k = a+ (k-1)/(N-1)(b-a)$  

-> Differentiation matrices *dx*, *dxx* are *sparse* (bi- or tridiagonal)

-> Error is $E \approx N^{-d}$ ($d=2$ for 2nd order) 

- Chebyshev collocation method

-> *Pseudo-mesh* $x = [x_1, x_2, ... x_N]$  (actally collocation points) is non-uniform, with clusering around boundaries $a,b$.

-> Differenciation matrices *dx*, *dxx* are *full* 

-> Error is $E \approx e^{-N}$ : *spectral precision !*

- Implementation in **easystab** project : **dif1D.m** (including finite-difference, chebyshev and many more !)

## Example : heat equation

See example program here : [diffusion_eigenmodes.m]()

Remarks : 

- Chebyshev is more precise than finite difference

- the boundary conditions lead to two infinite eigenvalues (which are easily discarded) 

## the Easystab project

Easystab ([http://basilisk.fr/sandbox/easystab/README]())  is "a door open for you to study the stability and bifurcation of physical systems using octave/matlab, mainly in the domain of fluid mechanics".

- Initiated by Jérôme Hoepffer (UPMC, Paris)

- Hosted by the website of Basilisk code (Stéphane Popinet, UMPC, Paris)

- Maintained as a wiki : upload your codes, they will including the resulting figures and the comments in readable format.
(You can even run the codes through the website without launching octave/matlab !)

- The engine also allows to share lecture notes in .md format...

- This is a collaborative project. Please contribute !




# Local stability analysis for parallel flows

## Starting equations

Consider a flow field as follows :

$$
[{\bf u} , p] = [{\bf u}_0 , p_0] + \epsilon [{\bf u'} , p']
$$

- $[{\bf u}_0 , p_0]$ is the (time-independant) *base flow* 

- $[{\bf u}' , p']$ is the perturbation, which is governed by the Linearized Navier-Stokes equations (LNSE) :

$$
\partial_t 
\left[
\begin{array}{c} 
{\bf u}'
\\ 
0 
\end{array}
\right]
= 
\left[ 
\begin{array}{c} 
- \left( {\bf u}_0 \cdot \nabla {\bf u}' + {\bf u}' \cdot \nabla {\bf u}_0 \right)  - \nabla p' + \frac{2}{Re}  \nabla \cdot {\bf{D}}({\bf u}')
\\
\nabla \cdot {\bf u} 
\end{array}
\right]
\equiv 
\mathcal{ LNS}([{\bf u}' ; p'])
$$

(Here ${\bf{D}}({\bf u}) = 1/2( \nabla {\bf u}' + \nabla {\bf u}'^T )$ is 
the rate of deformation tensor)


## Local approach

- In the *local approach* the base-flow is assumed *parallel*, with the form

$$
[{\bf u}_0 , p_0] = [U(y) ,0,0]$$

- The eigenmode ansatz (for 2D perturbations) is as follows :

$$ 
[{\bf u}' , p'] = [\hat{u}(y) ,\hat{v}(y),\hat{p}(y) ] e^{i k x} e^{\lambda t}
$$

- The linearized Navier-Stokes equations (primitive form) are as follows :

$$
\lambda
\left[
\begin{array}{c} 
\hat{u}
\\ 
\hat{v}
\\
0
\end{array}
\right]
= 
\left[ 
\begin{array}{cccc} 
- i k U \hat{u} & - (\partial_y U) \hat{v} & - i k \hat{p} 
& + Re^{-1} (\partial_{yy} - k^2 ) \hat{u}
\\
& - i k U \hat{v}  & -\partial_y \hat{p} 
&+ Re^{-1} (\partial_{yy} - k^2 ) \hat{v}
\\
i k \hat{u} & + \partial_y \hat{v} &&
\end{array}
\right]
$$

- The equations can also be set in reduced form (Orr-Sommerfeld equation) by 
introducing a streamfunction $\hat \psi$ :


$$
\left( U - c \right) \left( \partial_{yy} - k^2 \right) \hat{\psi} 
= \frac{1}{i k Re} \left( \partial_{yy} - k^2 \right)^2 \hat{\psi} 
$$

(here $c = \omega/k \equiv i \lambda / k$)


## Illustration for the tanh shear layer

- Inviscid approach (from Raylegh equation) 
[KH_temporal_inviscid.m]() 

-> We find un unstable mode for $k<1$ with maximum amplification rate for $k \approx 0.5$.

- Viscous approach (from primitive equations)
[KH_temporal_viscous.m]()

-> We observe that the viscosity restabilizes the 

## Illustration for the Plane Poiseuille flow
[poiseuille_uvp.m]()

-> a viscous instability (Tollmien-Schlishting wave) exists for $Re>5772$ in accordance with litterature.


# Global stability approach

## Method



- In the *global approach*, the base flow depends on two spatial coordinates:

$$
[{\bf u}_0(x,y) , p_0(x,y)]
$$

This base flow is solution of Steady Navier-Stokes equations

$$
\mathcal{NS}\left(\left[
\begin{array}{c} 
{\bf u}_0
\\ 
p_0 
\end{array}
\right]\right)
\equiv 
\left[ 
\begin{array}{c} 
- \left( {\bf u}_0 \cdot \nabla {\bf u}_0 \right)  - \nabla p_0 + \frac{2}{Re}  \nabla \cdot {\bf{D}}({\bf u}_0)
\\
\nabla \cdot {\bf u}_0
\end{array}
\right]
= 
\left[
\begin{array}{c} 
0
\\ 
0 
\end{array}
\right]
$$



## Base flow computation

-> **Newton method**

- Suppose we know a *guess* $[{\bf u}_0^g,p_0^g]$ for the base flow

- Assume small perturbations : 
$[{\bf u}_0,p_0] \approx [{\bf u}_0^g,p_0^g] + [\delta {\bf u}_0,\delta p_0]$

- Injecting in NS equations leads to 

$$
\mathcal{ NS }(  {\bf u}_0^g,p_0^g) 
+ \mathcal{ LNS }  [\delta {\bf u}_0, \delta p_0] \approx [{\bf 0}, 0]
$$

- Inverting leads to 
$$
[\delta {\bf u}_0; \delta p_0] = - (\mathcal{ LNS } )^{-1} \mathcal{ NS }(  {\bf u}_0^g,p_0^g) 
$$

- Repeat until convergence

NB convergence is quadratic (typically less than 10 iterations)

## Stability analysis

Once base-flow is determined, we repeat the stability analysis


$$
[{\bf u} , p] = [{\bf u}_0 , p_0] + \epsilon [\hat{\bf u} , \hat{p}] e^{\lambda t} + c.c.
$$


- $[\hat{\bf u}(x,y) , \hat{p}(x,y)]$ is the eigenmode, which is governed by the Linearized Navier-Stokes equations (LNSE) :

$$
\lambda {\mathcal B}  
\left[
\begin{array}{c} 
\hat{\bf u}
\\ 
\hat{p} 
\end{array}
\right]
= 
\left[ 
\begin{array}{c} 
- \left( {\bf u}_0 \cdot \nabla \hat{\bf u} + \hat{\bf u} \cdot \nabla {\bf u}_0 \right)  - \nabla \hat{p} + \frac{2}{Re}  \nabla \cdot {\bf{D}}(\hat{\bf u})
\\
\nabla \cdot \hat{\bf u}
\end{array}
\right]
\equiv 
\mathcal{ LNS}\left[
\begin{array}{c} 
\hat{\bf u}
\\ \hat{p}
\end{array}
\right]
$$

Here
$$\mathcal{B} = \left[
\begin{array}{cc} 
1 & 0
\\ 
0 & 0
\end{array}
\right]
$$

Or in eigenvalue form :

$$( \mathcal{LNS} - \lambda \mathcal{B} )  [\hat{\bf u}, \hat{p}] = [{\bf 0}, 0]
$$

Note that matrix $\mathcal{LNS}$ is the same as for the Newton method !

## Numerical methods

- Prefered method : *Finite elements*  

weak form of NS equation :
$$\forall [\mathbf v, q], \qquad  \int_{\Omega} [\mathbf v, q] \cdot \mathcal{ NS }(  {\bf u},p) d S = 0
$$ 

after integration by parts a very convenient form is obtained :

$$\forall [\mathbf v, q], \qquad  \int_{\Omega} 
\left( - {\mathbf v} \cdot ( {\mathbf u} \cdot \nabla) \mathbf u) 
+ p \nabla \cdot \mathbf v + q \nabla \cdot \mathbf v - 2 Re^{-1}
{\bf{D}}(\hat{\bf u}) : {\bf{D}}(\hat{\bf v}) \right) d S = 0
$$ 

- Prefered software : [FreeFem](https://freefem.org/) (Université 
Paris-Sorbonne)

(
See how the weak form in implemented in [Newton solver for a 2D flow](https://gitlab.com/stabfem/StabFem/blob/master/SOURCES_FREEFEM/Newton_2D.edp) and [Eigenvalue solver for a 2D flow](https://gitlab.com/stabfem/StabFem/blob/master/SOURCES_FREEFEM/Stab2D.edp)
)


## StabFem

FreeFem is very convenient for stability studies (and used worldwide for this)  *but* it is not suited to functional programming (the language is interpreted, not compiled) and its graphical interface is limited.

-> We designed a Matlab/Octave to Freefem : **StabFem**

Developed as a collaborative project with same philosophy as easystab (please contribute !)


**Methods currently available :**

- Base-flow computation using Newton method (including arclength continuation)

- Eigenvalue solver (including adjoint modes, structural sensitivity, and interative spectrum exploration)  

- DNS (time-stepping integration)

- Nonlinear stability analysis (weakly nonlinear, Harmonic Balance/Self consistent)

- Powerful and versatile mesh adaptation

- (...)

**Models currently available :**

- Incompressible flows (bluff-body wakes, 

- Compressible flows (whistling jet instabilities, acoustic pipes, ...) 

- Fluid-structure interactions (VIV) and 

- Free surface flows (rising or attached bubbles, whirlpools,...)

- (Currently limited to 2D geometries)


StabFem Showroom Website : [https://stabfem.gitlab.io/StabFem/]()

StabFem Sources Website : [https://gitlab.com/stabfem/StabFem]()


**Documentation :**

- research paper by [Fabre et al. AMR, 2019](https://gitlab.com/stabfem/StabFem/blob/master/99_Documentation/ARTICLE_STABFEM/amr_070_06_060802.pdf) 

- A workshop presentation [FreeFem++ days 2018](
https://gitlab.com/stabfem/StabFem/blob/master/99_Documentation/PRESENTATIONS/Beamer_13dec2018_handout.pdf)

- A pdf documentation [Doc StabFem](https://gitlab.com/stabfem/StabFem/blob/master/99_Documentation/MANUAL/main.pdf) (in progress, currently not up to date...)

(- All codes are internally documented, type 'help' in matlab/octave)


## Illustration for the wake of a cylinder


- [Stability analysis for the wake of a cylinder](https://stabfem.gitlab.io/StabFem/stable/cylinder/CYLINDER_LINEAR.html)

- [DNS of the wake of a cylinder](https://stabfem.gitlab.io/StabFem/stable/cylinder_dns/SCRIPT_DNS_EXAMPLE.html)


## Other examples :


- [Polygonal patterns in a free-surface swirling flow](https://stabfem.gitlab.io/StabFem/stable/ROTATING_POLYGONS/SCRIPT_POLYGONS.html)

- [Acoustic radiation from an open pipe](https://stabfem.gitlab.io/StabFem/stable/acoustic_pipes/SCRIPT_DEMO_ACOUSTIQUE.html)

- (...)



# Conclusions

- Numerical resolution of stability problems is now becoming a game

- Collaborative projects are providing you nice set of programs to do this

- So come and play with us !


Thanks for your attention 

