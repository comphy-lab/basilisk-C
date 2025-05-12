
**Numerical methods for flow instability problems **

**D. Fabre, IMFT (Toulouse, France) **

This document belongs to the lecture notes for the [M2-DET](http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md) course, D. Fabre, nov. 2018-dec. 2023 : Lecture 3.

(This document is adapted from the lecture given at Sapporo university in june 2019 available [here](http://basilisk.fr/sandbox/easystab/Lecture_Sapporo_june2019.md) )

# Linear dynamical systems and linear PDE problems

## Linear dynamical systems

(extract from [LinearSystems.md]() )

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

Just as for usual linear dynamical systems one can look for eigenvalue/eigenvector solutions, by considering the *generalized eigenvalue problem*
$$
\lambda B \hat{X} = A \hat{X}
$$



- If $B$ is invertible, this is equivalent to a dynamical system with 
$\frac{d X}{dt} = ( B^{-1} A ) X$. Hence the eigenvalues of the generalized problem are simply the eigenvalues of the matrix $B^{-1} A$, and the general solution has the form given above (Eq. 1).

- If $B$ is non-invertible, one can still look for eigenvalue/eigenvector solutions. There are again $N$ eigenpairs, but **some of the eigenvalues are infinite**.

In this case the general solution can still be written as


$$ 
X(t) = \sum_{\lambda_n \in {\mathcal Sp}^*} c_n \hat X_n e^{\lambda_n t} 
$$

Where ${\mathcal Sp}^*$ is the set of *finite eigenvalues*.

The case where $B$ is non-invertible generally arises when considering a system of equation containing both dynamical equations and constraints on the variables. Then, the eigenvectors corresponding to infinite eigenvalues correspond to solutions which does not respect the constraint, and it is logical that they cannot appear in the general solution.

More explanations can be found in the companion document [LinearSystems.md](). 
See in particular 
[this example](http://basilisk.fr/sandbox/easystab/LinearSystems.md#an-example) to understand the significance of infinite eigenvalues.

## Numerical resolution of eigenvalue (and generalized eigenvalue) problems

We review here the main classes of numerical methods to solve eigenvalue problems. We restrict here to simple dynamical systems of the form $dX/dt = A X$ but most ideas can be generalized to generalized eigenvalue problems.


### Using the characteristic polynomial

It is well known that the eigenvalues are the root of the characteristic polynomial $P(\lambda) = det ( A - \lambda I)$ (or $det ( A - \lambda B)$ for generalized eigenvlaue problems).

This method is mathematically elegant but practically efficientless for $N>3$ !

### Direct method : matrix diagonalization.

The idea consists of factorizing the matrix in the following way.

$$
A = U S U^{-1}
$$
where $S = diag(\lambda_n)$ is diagonal, and $U =  [\hat{X}_1, \hat{X}_2, ..., \hat{X}_N]$ is a matrix containing the eigenvectors sorted in columns. 

There are direct algorithms to do this factorization. 

- Fortran/c++ : **Lapack** library ;  

- Matlab/Octave : **eig** (actually based on Lapack)

- Python : **scipy.linalg.eig**

However these algrithms are costly : the algorithmic complexity (number of operations) is of order $N^3$.
Hence a full diagonalization is only achievable for relatively small systems ($N<1000$).

Moreover, it is generally useless to compute all eigenvalues to decide for stability/instability, since only the one with largest real part is needed. 


### Iterative methods 

Iterative methods are alternative class of methods which are much faster and efficient when it only required to compute a single eigenvalue or a few ones, instead of the whole spectum.

The simplest is the power method defined by the folloqing algorithm 

#### Power method

Consider the ***Power method algorithm*** as follows:

**Require :** initial condition $X_0$

**Initialization :** $X \leftarrow X_0$

**Loop until convergence :**
 
$\quad X \leftarrow A X$ 
 
$\quad X \leftarrow \frac{X}{||X||}$

**End loop**

One can prove easily that the algorithm converges towards the eigenmode associated to the eigenvalue with largest absolute value.


#### Variants of the power method 

- A variant of the previous method is the ***inverse power method*** which consists of applying $\quad X \leftarrow A^{-1} X$ at each step. Then one can prove that the method converges towards the eigenmode associated to the eigenvalue with smallest absolute value.

- Another variant is to apply the inverse iteration to the *shifted problem* $(A- \sigma B) \hat{X} = (\lambda- \sigma) B \hat{X}$ where $\sigma$ is a given complex number called the *shift*. This leads to the following algorithm :

 ***"Shift-and-invert power method algorithm"*** 

**Require :** initial condition $X_0$ and shift $\sigma$

**Initialization :** $X \leftarrow X_0$
 
**Loop until convergence :**
 
$\quad X \leftarrow (A-\sigma B)^{-1} B X$ 
 
$\quad X \leftarrow \frac{X}{||X||}$

**End loop**

One can show again that the algorithm converges towards the eigenmode $\hat X$ corresponding to the eigenvalue such as $|\lambda - \sigma|$ is minimum. The convergence is fastest if the shift is carefully selected to be very close of an eigenvalue.
Hence this method is extremely efficient if one wants to compute a single eigenvalue and one already knows approximately this eignevalue.


#### Arnoldi methods ####

The drawback of the power method and its variants is that they provide only a single eigenvalue/eigenvector. Arnoldi methods are an extension of these methods which allow to compute a limited number of eivenvalues/eigenvectors (typically $n= 10$). Considering power iteration (***Arnoldi power method***), 
the idea is to consider not only the final result of the iteration ($A^n X_0$ for the power method) but the whole sequence $\{ X_0, A X_0, A^2 X_0, ... A^n X_0 \}$ which form a subspace of order $n$ called *Krylov space*. A diagonalization procedure is then applied to this space, allowing to extract an approximation of the $n$ eigenvalues with largest absolute value. 

The same idea can also be used replacing power iteration by inverse power iteration (***Arnoldi inverse method***) leading to the $n$ eigenvalues with smallest absolute value, or by the shift-and-inverse iteration (***Arnoldi shift-and-invert method***) leading to the $n$ eigenvalues closest to the shift. See [this page](https://en.wikipedia.org/wiki/Arnoldi_iteration) for more details.

These methods are currently the prefered methods in Stability analyses. Current implementations are: 

- Fortran/c++ : **Arpack** library ; 

- Matlab/Octave : **eigs** (actually based on Arpack). 

- Python : **scipy.sparse.linalg.eigs**

- For large systems handled with parallelization: **SLEPc** library (yet no simple interface with Matlab/Octave/Python)

## Linear PDE problems

### Definition
Consider a linear partial difference equation (PDE) as follows :

$$
{\mathcal B} \frac{\partial F(x,t)}{ \partial t}  = {\mathcal A} F(x,t)
$$

where $F(x,t)$ is defined on an interval $x\in \Omega = [a,b]$ (possibly infinite)  and for $t \in [0, \infty]$,  and $\mathcal A$ and $\mathcal B$ are *linear operators* involving $x$-derivatives ( operators $\partial / \partial x$ , $\partial^2 / \partial x^2$, ... ). 

If the interval is bounded then one also has to specify *boundary conditions*; typically either Diriclet, i.e. $F(x=\{a,b\},t) = 0$ or Neuman ($(\partial F/\partial x)_{x=\{a,b\}} = 0$).

<span style="color:green">
NB : to set the problem in a rigourous way, one will typically require the initial condition and the solution to belong to some functional space ${\mathcal H}$. A convenient choice is the *Sobolev space* of functions defined on $\Omega$ which are continuously derivable up to a certain order, are square-integrable, and verify the boundary conditions.
</span>

### Mathematical analysis

It is tempting to generalize the notion of eigenmodes for linear PDE problems. The situation is however more complex, and two kind of solutions can be encountered :

- First, on some situations, one can find *regular eigenfunctions* with the form $F(x,t) = \hat{F}(x) e^{\lambda t}$, which are the solutions of the *continuous eigenvalue problem*:

$$
\lambda {\mathcal B} \hat{F} = {\mathcal B} \hat{F}
$$

A simple example is the 1D heat equation in bounded domain, see [Exercice 2](http://basilisk.fr/sandbox/easystab/NumericalMethodsForEigenvalueProblems.md#d-heat-equation). 



<span style="color:green">
NB: There is a particular class of problems called **Sturm-Liouville problems**  for which it can be demonstrated that there is an countable infinity of such solutions noted $\{\lambda_n,\hat{F}_n(x)\}$, with $n= 1,2,...$.  Moreover in such problems the eigenfunctions form a *complete basis* of the space ${\mathcal H}$, meaning that the solution can always be expressed as 
$$ 
F(x,t) = \sum_{n=1}^{\infty} c_n \hat F_n(x) e^{\lambda_n t}  .
$$
The 1D heat equation example belongs to this class. 
See document [LinearSystems.md]() for more on this point.
</span>

- Secondly, in other situations, one can find *generalized eigenfunctions* with the form $F(x,t) = \tilde{F}_\lambda(x) e^{\lambda t}$, when $\lambda$ belongs to some continuous set called the *continuous spectrum*.
These are weak solutions of the problem, in the sense that they do not belong to the functional space ${\mathcal H}$  (they are either discontinous, non-square-integrable, or don't verify the boundary conditions). 
These solutions are not acceptable solutions of the problem when considered alone but a continuous sum of such solutions is acceptable.

A simple example is studied in [Exercice 4](http://basilisk.fr/sandbox/easystab/NumericalMethodsForEigenvalueProblems.md#an-example-involving-a-continuous-spectrum-transport-equation-in-an-infinite-domain). 

- In general, the solution of a linear PDE problem can be expressed as

$$ 
F(x,t) = \sum_{\lambda_n \in {\mathcal Sp_d}} c_n \hat F_n(x) e^{\lambda_n t} 
+
\int_{\lambda \in Sp_c} c(\lambda) \tilde{F}_\lambda(x)  e^{\lambda t} d \lambda
$$

Where $Sp_d$ and $Sp_c$ are the **discrete spectrum** and **continuous spectrum**, respectively. Depending on the cases, the discrete spectrum can be infinite-denombrable, finite or empty; 
the continuous spectrum is infinite and non-denombrable (but may also be empty). 

*Note that the continuous spectrum, if present, is generally located in the half-plane $Re(\lambda) \le 0$, so it cannot lead to an instability. The contribution of the continuous spectrum is generally decaying in time (the decay being algebraic instead of exponential).*

<span style="color:green">
Mathematically, the establishment of the general solution given above is based on the theory of *Laplace transform* in time.
</span>

## Discretization of a continuous problem

### Procedure

The **discretization** is the reduction of a continuous problem to a discrete problem of large dimension $N$.

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

- These matrices are used to build the discretized version $B$ of the operators $\mathcal A$ and $\mathcal B$.

*Remark* : 

The discretized problem necessarily admits a finite number $N$ of eigenmodes which are used to build an approximation of the solution of the continuous problem. In practice the computed eigenvalues will be of two kinds :

- A small number of *physical eigenvalues* which do not depend upon the discretization and correspond to the discrete modes with simplest structure.

- A large number of *spurious eigenvalues* which depend upon the disctetization, and correspond either to badly converged discrete modes or to a discretized version of the continuous modes *(if the discretization method is well designed the spurious modes are generally in the half-plane $Re(\lambda) \le 0)$ so they do not lead to unstable behaviour).*


### Discretization methods

- Finite difference method

-> mesh $x = [x_1, x_2, ... x_N]$ is *uniform* : $x_k = a+ (k-1)/(N-1)(b-a)$  

-> Differentiation matrices *dx*, *dxx* are *sparse* (bi- or tridiagonal)

-> Error is $E \approx N^{-d}$ ($d=2$ for 2nd order) 

- Chebyshev collocation method

-> *Pseudo-mesh* $x = [x_1, x_2, ... x_N]$  (actally collocation points) is non-uniform, with clusering around boundaries $a,b$.

-> Differenciation matrices *dx*, *dxx* are *full* 

-> Error is $E \approx e^{-N}$ : *spectral precision !*

- Implementation in **easystab** project : **dif1D.m** (including finite-difference, chebyshev and many more !)

### Example : heat equation

See example programs here : 

[diffusion_eigenmodes.m]() (using finite-difference discretization)

[diffusion_eigenmodes_Chebyshev.m]() (using Chebyshev pseudo-spectral discretization)

Remarks : 

- Chebyshev is MUCH more precise than finite difference

- the boundary conditions lead to two infinite eigenvalues (which are easily discarded) 

## the Easystab project

Easystab ([http://basilisk.fr/sandbox/easystab/README]())  is "a door open for you to study the stability and bifurcation of physical systems using octave/matlab, mainly in the domain of fluid mechanics".

- Initiated by Jérôme Hoepffer (UPMC, Paris)

- Hosted by the website of Basilisk code (Stéphane Popinet, UMPC, Paris)

- Maintained as a wiki : upload your codes, they will including the resulting figures and the comments in readable format.
(You can even run the codes through the website without launching octave/matlab !)

- The engine also allows to share lecture notes in .md format...

- This is a collaborative project. Please contribute !

# Exercices for lecture 3:

## A. Linear algebraic systems

### Exercice 1.

Consider a differential algebraic system defined as 
$$\frac{d}{d t} B X
= A X.
$$

  a. Compute the eigenvalues and eigenvectors, and give the general solution $X(t)$, for the following case :
$$
A = \left[ \begin{array}{cccc} 
1 & 0 & 0 & 0  \\ 1 & -2 & 1 & 0 \\ 0 & 1 & -2 & 1  \\ 0 & 0 & 0 & 1
\end{array}
\right] 
\quad 
B = \left[ \begin{array}{cccc}
0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 &1 & 0 \\ 0 & 0 & 0 & 0
\end{array} \right]
$$
  b. Same questions for the case:
$$
A = \left[ \begin{array}{ccc} 
-1 & 1 & 1 \\ 1 & -2 & 1 \\ 1 & 1 & 0
\end{array}
\right] 
\quad 
B = \left[ \begin{array}{ccc} 
1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0
\end{array}
\right]
$$

  c. (preparation for lecture 4) Same questions for the case:

$$
A = \left[ \begin{array}{ccc} 
k & 0 & 0 \\ 0 & -k & 0 \\ 0 & 0 & (\rho_2-\rho_1) g - \gamma k^2 
\end{array}
\right] 
\quad 
B = \left[ \begin{array}{ccc} 
0 & 0 & 1 \\ 0 & 0 & 1 \\ \rho_2 & \rho_2 & 0
\end{array} \right]
$$

It is advised to work these exercices in three ways : $(i)$ by writing down a characteristic polynomial, $(ii)$ by reducing the order of the system by incorporating the constraints (as in [this example](http://basilisk.fr/sandbox/easystab/LinearSystems.md#an-example) ) and $(iii)$ With a computer, using **eig** and/or **eigs** under Matlab/Octave/Python as in [this program](http://basilisk.fr/sandbox/easystab/demo_eig.py)

For case $(c)$ use $\rho_1 = 1/2,\rho_2 = 1, g = 1, \gamma = 1, k = 1/2$.

## B. Continuous problems

### 2. 1D heat equation 

Consider the heat equation in a finite domain for $T(x,t)$:

$$
\frac{\partial T}{\partial t} = \mu \frac{\partial^2 T}{\partial x^2} 
\qquad \mathrm{ with } \quad x \in [0,L]
$$

with boundary conditions 
$$
T(0,t) = T(L,t) = 0
$$

and initial condition $T(x,0)= T_0(x)$.

a. Considering eigenmode solutions of the form $T(x,t) = \hat{T}_n(x) e^{\lambda_n t}$, give the expression of the eigenvalues $\lambda_n$ and corresponding "eigenfunctions"  $\hat{T}_n(x)$.

b. Deduce that the general solution is:

$$
T(x,0) = \sum_{n=0}^\infty c_n sin \left( \frac{n \pi x}{L} \right) 
\exp \left( - \frac{\mu n^2 \pi^2 t }{L^2} \right).
$$

c. Explain how the coefficients $c_n$ can be deduced from the initial condition $T_0(x)$.



### 3. Heat conduction in cylindrical coordinates

Consider the heat equation, for a material with density $\rho$, thermal capacity $c_p$ and conductivity $\kappa$:
$$
\rho c_p \frac{d T}{d t} = \nabla \cdot \left( \kappa \nabla T \right)
$$
We want to solve this problem for the case of a cylindrical rod of radius $R$. We suppose that the properties of
the rod are inhomogeneous but respect the cylindrical geometry ($\rho = \rho(r) ; c_p =  c_p(r), \kappa = \kappa(r)$. 
We thus expect that the solution $T$ will respect the same symmetries and this look for a solution under the form
$T = \hat T(r) e^{\lambda t}$. Assume homogeneous Dirichlet conditions at the boundary of the rod ($T(R,t)=0$).

a. Write a differential eigenvalue problem verified by the eigenmodes/eigenvalues ($\hat T ; \lambda$).

b. Is this problem of Sturm-Liouville type ? Deduce the sign of the eigenvalues. 

c. In the case of a uniform material ($\rho,c_p,\kappa$ are constant), show that the eigenmodes $\hat T(r)$ can be expressed
in terms of the Bessel function $J_0(\xi r)$ where $\xi$ is a "radial wavenumber" directly related to the eigenvalue $\lambda$.
Using the boundary condition at $r=R$, give an expression of the eigenvalues $\lambda_n$ in terms of the zeroes of the bessel function (noted $j_{0,n}$). 

d. Write a program to solve numerically this problem
(you may start from the example program [diffusion_eigenmodes.m]() and adapt it for polar coordinates)

### 4. An example involving a continuous spectrum : Transport equation in an infinite domain

Consider the transport equation in a 1D infinite domain :

$$
\frac{\partial Y}{\partial t} + c \frac{\partial Y}{\partial x} = 0
\qquad \mathrm{ with } \quad x \in [-\infty,\infty]
$$

With a initial conditions $Y(x,0) = Y_0(x)$ ; $Y_0 \in {\mathcal H}^2$ where ${\mathcal H}^2$ is the Sobolev space of functions which are continuous, derivable up to order 2, finite, and square-integrable. 

a. Show that the general solution is $Y(x,t) = Y_0(x-c/t)$. 

b. Consider now eigenmode-like solutions of the form 
$Y(x,t) = \hat{Y}_\lambda(x) e^{\lambda t}$. Give the expression of  the eigenfunctions $\tilde{Y}_\lambda(x)$. Are these functions physically acceptable solutions (do they belong to the Sobolev space ${\mathcal H}^2$) ? 
Justify that if $\lambda$ is pure imaginary ($\lambda = -i \omega$ with $\omega \in R$) these functions are "nearly acceptable".

c. Justify that the general solution can be expanded as 
$Y(x,t) = \int_{\mathcal C} c(\lambda) \tilde{Y}_\lambda(x) e^{\lambda t} d \lambda 
= \int_{-\infty}^{+\infty}   c(\omega)  e^{i (k x - \omega t) } d \omega.$ 

where $\mathcal C$ is a contour going from $-i\infty$ to $+i\infty$.

<span style="color:green">
NB : This solution can also be obtained by starting with a *Fourier transform in space*. The approach proposed in this exercice is closer to a *Laplace transform in time*.
</span>
