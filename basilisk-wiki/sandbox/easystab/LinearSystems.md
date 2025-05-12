**(This document belongs to the lecture notes for the [M2-DET](http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md) course, D. Fabre, nov. 2018-dec. 2023)**

This documents gives mathematical support for the study of linear dynamical systems. The notions will
be used throughout the course, especially in lecture 1 and lecture 9.

<span style="color:green">
Paragraphs in green are mostly useful for lectures 9 and 10, and can be skipped at first lecture.
</span>



#Linear dynamical systems

In this section we investigate the behaviour of an exactly linear dynamical system  possessing a fixed point $X_0 =0$, written with the form:

$$
\frac{d X}{dt} = A X
$$

Where $X(t) = [x_1(t);x_2(t);.. ;x_N(t)]$ is a vector (column-vector) of dimension $N$ and $A$ is a square matrix of dimension $N$.

The matrix $A$ (which coincides with the gradient of the flow of the underlying dynamical system) defines a linear application defined on  $\mathbb{R}^N \rightarrow \mathbb{R}^N$ (if the matrix is real) or on  $\mathbb{C}^N \rightarrow \mathbb{C}^N$ (if the matrix is complex). However, even if the matrix is real, the solution of the dynamical system may have to be expressed in terms of complex solutions, so it is better to consider it as a complex application in all cases.

## A few definitions

In the sequel we will require to define a scalar product $<X,Y>$ and a norm $||X||$ on $\mathbb{C}^N$. The simplest choice (canonical) is as follows (where the overbar denotes the complex conjugate):
$$
<X,Y> = \overline{X}^T \cdot Y = \sum_{i=1}^N \overline{x}_i y_i,
$$
$$
||X|| = {\left(<X,X>\right)}^{1/2} = {\left( \sum_{i=1}^N |{x}_i|^2 \right)}^{1/2}.
$$
<span style="color:green">
We call $A^\dag$ the *adjoint* of the matrix $A$, defined by the property :
$$
\forall (X,Y), \quad   < A^\dag x, y> = <x, A y>.
$$
It is clear that $A^\dag = \overline{A}^T$, namely the adjoint matrix is the conjugate transpose of the matrix 
(also called the hermitian transpose).    
A matrix is said to be *self-adjoint* (or Hermitian) if $A^\dag = A$.   
For real matrices this notion is equivalent to saying that the matrix is *symmetric*.
</span>

## Eigenvalues and eigenmodes

We look for solutions under the "eigenmode" ansatz :
$$
X(t) = \hat{X} e^{\lambda t}
$$

The dynamical system thus results to a eigenvalue problem :

$$
\lambda \hat{X} = A \hat{X}
$$

The solutions $\lambda$ are the *eigenvalues* and the corresponding 
$\hat{X}$ the *eigenvectors*. The full set of eigenvalue is called the *spectrum* of the matrix.

We suppose in this section and in the next one that **the eigenvalues are all distinct**.
The exceptional case where two eigenvalues (or more) are identical requires specific treatment,
but can be set aside for the moment.  

There are $n$ couples of eigenvalue/eigenvector which can be numbered  $(\lambda_n;\hat X_n)$ for $n=1..N$. 

The set of eigenvectors ${\{ \hat X_n \}}_{n=1..N}$  forms a *basis* of $\mathbb{R}^N$ (or $\mathbb{C}^N$)
They can be normalised by  $||\hat{X}_n||=1$, leading to a normed basis (not necessarily orthogonal).


- If the matrix A is symmetric (or hermitian in the complex case), the eigenvalues 
are all real and the eigenmodes form an orthogonal basis. Namely, they satisfy the orthogonality condition :
$$
<\hat X_n, \hat X_m> = \delta_{mn} 
$$

- If $A$ is real, nonsymmetric, the eigenvalues may be either real values (noted $\lambda_n$) 
or complex conjugate pairs (noted $\lambda_n = \sigma_n \pm i \omega_n $).
The eigenmode associated to complex eigenvalues can be written under the form $\hat{X}_n = U_n \pm i V_n$
where $U_n$ and $V_n$ are real vectors which can be chosen as orthogonal (but not necessarily of same norm).

- <span style = "color:orange"> If the starting problem is conservative (in the Hamiltonian sense), the full spectrum is symmetric with respect to the real and imaginary axes. This means that the eigenvalues are either couples of opposite real numbers, ($\pm \lambda_r$), couples of  opposite imaginary numbers §\pm i \lambda_i$ or quadruplets of complex numbers 
$\pm \lambda_r \pm \lambda_i$.
</span>

- <span style="color:green"> In the general case, we can define the *adjoint eigenmodes* defined as the solutions 
of the adjoint eigenvalue problem : 

<span style="color:green">
$$ 
\overline \lambda \hat X^\dag = A^\dag \hat X^\dag 
$$
It can be proved that the eigenvalues of this adjoint problem are idential to those of the matrix $A$ (and hence also noted $\lambda_n$), while that the associated adjoint eigenmodes satisfy the generalized orthogonality condition:
$$ <\hat X^\dag_n, \hat X_m> = \delta_{mn}$$
The set of adjoint eigenmodes forms a basis of $\mathbb{R}^N$ (or $\mathbb{C}^N$)called the *adjoint basis*, also called the *dual basis* of the eigenmode basis.
<\span>


## Solution of the initial-value problem

Given an initial condition $X(0) = X_0$, the solution can be projected along the basis of eigenmodes.   

$$ 
X(t) = \sum_{n=1}^N c_n \hat X_n e^{\lambda_n t} \qquad \mathrm{ (Eq. 1)}
$$

- If the matrix A is real, nonsymmetric, the solution (Eq. 1) can also be written under the equivalent form :

$$ 
X(t) = \sum_{\lambda_n \in \mathbb{R} } c_n \hat X_n e^{\lambda_n t} + \sum_{\lambda_n \in \mathbb{C} } \left(c_{n,c} U_n \cos \omega_n t 
+ c_{n,s} V_n \sin \omega_n t \right) e^{\sigma_n t}.
$$

-  If the matrix A is real, symmetric, the orthogonality condition allows to compute the coefficients
$c_n$ by simply projecting the initial condition along the corresponding eigenmodes :
$$
c_n = < \hat{X}_n , X_0 >
$$


- <span style="color:green">  In  the general case, we have to use the *Adjoint basis* to project the initial condition over the basis of eigenmodes.
This leads to :
$$
c_n = \frac{< \hat{X}_n^\dag , X_0 >}{< \hat{X}_n^\dag , \hat{X}_n >}
$$
</span>

### Stability criteria


We supose that the eigenvalues are sorted in the order of decaying real part, namely
$$
Re(\lambda_1) \geq  Re(\lambda_2) \geq  ... \geq  Re(\lambda_N)
$$



- The system is unstable if there is at least one eigenvalue with positive real part, namely : if $Re(\lambda_1)>0$

- Reversely, the system is stable (more precisely exponentially stable) if all eigenvalues have negative real part, namely : if $Re(\lambda_1)<0$.

- If $Re(\lambda_1)=0$ the linear system is said to be marginally stable. 

<span style="color:green">
Note that if $A$ is symmetric, the condition for stability implies monotonous stability. If $A$ is nonsymmetric, on the other hand, the system may be exponentially stable but not monotonously stable. This situation corresponds to "transient growth" and will be reviewed in lecture 9.
<\span>


### Particular case of two zero eigenvalues.

As said above, the case of degenerated eigenvalues requires specific treatment. The only case which deserve attention in the perspective of stability studies is the case of two zero eigenvalues. We suppose here that $\lambda_1=0$ is an eigenvalue of order 2
and that all other eigenvalues are stable ($0>Re(\lambda_3)\ge Re(\lambda_4)\ge ...$).

In this case the eigenspace associated to the eigenvalue $0$ is of dimension 2. In general there exists one 
*true eigenmode* $\hat X_1$ associated to the eigenvalue (hence $A \hat X_1 =0$), and in addition we can find a 
*generalised eigenmode* $\tilde X_2$ so that $A \tilde X_2 = \alpha \hat X_1$.

In that case, the initial-value problem solution can be generalized to:

$$ 
X(t) = (c_1+a c_2 t) \hat X_1 + c_2 \tilde X_2 
+
\sum_{n=3}^N c_n \hat X_n e^{\lambda t} 
$$

If $\alpha \ne 0$, the solution is said to be algebraically unstable.


### Generalization : expression of the solution using matrix exponentials.
<span style="color:green">
In the general case, a systematic way to express the solution is to make use of the 
*diagonalisation* of the matrix, which is is the following expression :
$$
A = U D U^{-1}
$$
Here :   

- <span style="color:green"> $D$ is diagonal (if all eigenvalues are distinct) or includes *Jordan blocks* (if some eigenvalues are degenerated).    
- <span style="color:green"> $U$ is a matrix constructed by arranging the eigenvalues in columns.   
- <span style="color:green"> $U^{-1}$ is the matrix obtained by arranging the adjoint eigenvalues in lines.
(it is clear that if $A$ is symmetric (or hemitian), then $U^{-1} = U^T$).


<span style="color:green">
The solution can thus be expressed using the matrix exponential :
$$
X(t) = e^{At} X_0 = U e^{Dt} U^{-1} X_0.
$$
which is equivalent to (Eq. 1) if the eigenvalues are distinct. If not, it contains terms with the form $t e^\lambda_n t$ due to the Jordan blocks. 




## Numerical resolution of eigenvalue problems

Let us review the main numerical methods available to solve eigenvalue problems:

- Direct computation of characteristic polynomial : impossible, unless for very small dimensions $N$ !
- Diagonalisation of the matrix (function  **eig** of Matlab/Octave) to compute all eigenvalues.
- Power method to compute the eigenvalue with largest absolute value.
- Arnoldi methods : generalisation of the power method to compute a limited number of eigenvalues
(function  **eigs** of Matlab/Octave)

# Linear algebraic differential problems 

We call "linear algebraic differential problems" systems of equations which can be set under the form 

$$ B \frac{d X}{dt} = A X$$

Where $E$ is a square matrix ("weight matrix") of dimension N.

The eigenmode analysis leads to a generalised eigenvalue problem with the form :

$$
\lambda B \hat{X} = A \hat{X}
$$

- If $B$ is inversible, then the problem is equivalent to a dynamical system
$$
\frac{d X}{dt} = B^{-1} A X
$$

All results of the previous section, based on eigenvalues/eigenmodes of matrix $B^{-1} A$, apply.
We can use the matlab functions *eig(A,B)* and *eigs(A,B)* to solve such generalised problems.

- If $B$ is not inversible, we can still look for solutions in eigenmode form,  but it turns out
that there is at least one infinite eigenvalue.

In this case, the expression of the initial value problem given by (Eq. 1) can still be used,
but only eigenmodes corresponding to finite eigenvalues have to be retained in the solution


### An example

Consider the differential algebraic sytem $B dX/dt = A X$, with

$$
A = \left[ \begin{array}{ccc} 
1 & 0 & 0 \\ 0 & 2 & 0 \\ 1 & 1 & 1  
\end{array}
\right] 
\quad 
B = \left[ \begin{array}{ccc} 
1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0
\end{array}
\right] 
$$

This can be split into the $2\times 2$ dynamical system


$$
\frac{d}{dt}  \left[ \begin{array}{c} x_1 \\ x_2 \end{array} \right] = A 
 \left[ \begin{array}{c} x_1 \\ x_2 \end{array} \right]
 \qquad \mathrm{ with } \quad  
 A = \left[ \begin{array}{cc} 
1 & 0  \\ 0 & 2  
\end{array}\right]
$$
and the constraint $x_1+x_2+x_3 = 0$.

The general solution is 
$$
 \left[ \begin{array}{c} x_1 \\ x_2 \\x_3 \end{array} \right] 
 =
 \left[ \begin{array}{c} c_1 e^t \\ c_2 e^{2t} \\ -c_1 e^t-c_2 e^{2t} \end{array} \right] 
 =
 c_1
\left[ \begin{array}{c} 1 \\ 0 \\ -1 \end{array} \right]  e^{t} 
+ c_2 \left[ \begin{array}{c} 0 \\ 1 \\ -1 \end{array} \right]  e^{2t} 
$$

A Python program solving this problem using both **scipy.linalg.eig** (direct method) and **scipy.sparse.linalg.eigs** (iterative method) can be found 
[here](http://basilisk.fr/sandbox/easystab/demo_eig.py).


The program returns the following set of eigenvalues : 
$$\left[ \lambda_1, \lambda_2 , \lambda_3 \right] = \left[ \infty ,1, 2\right]$$ 
with the corresponding eigenvectors :
$$
\hat X_1 = \left[ \begin{array}{c} 0 \\ 0 \\  1 \end{array} \right],
\quad
\hat X_2 = \left[ \begin{array}{c} \sqrt{1/2} \\ 0 \\ -\sqrt{1/2} \end{array} \right], 
\quad 
\hat X_3 = \left[ \begin{array}{c} 0 \\ \sqrt{1/2} \\  -\sqrt{1/2} \end{array} \right]. 
$$ 
which is equivalent to the general solution given above up to a renormalization. Of course the eigenvector $\hat X_1$ does not verify the constraint so it cannot be part of the general solution. 





# Linear partial-difference problems.

In hydrodynamic stability  we will encounter problems for a function $U(x,t)$, with the form

$$
\partial U/\partial t = F(U, \partial U/\partial x, \partial^2 U/\partial x^2,...)
$$

where $F$ is generally a nonlinear operator containing the partial derivatives with respect to $x$.

The function $U(x,t)$ is generally defined on an interval $x\in [a,b]$ which may be infinite. 
The resolution of the problem consists of finding the solution $U(x,t)$ for $t>0$, for specified *initial condition* $U(x,0)= U_0(x)$ and *boundary conditions* at $x=a$ and $b$.

In the case where $F$ depends linearly on $U$ and its spatial derivatives, we can consider the problem
as a *generalised linear dynamical system* , namely :

$$
\partial U/ \partial t = A \cdot  U 
$$
where $A$ is a linear operator (corresponding to the linearisation of $F$) acting on the function $U(x)$.

In this paragraph we adress the two following questions :

- It is possible to generalises the notion of linear eigenmode under the form 
$$
U(x,t) = \hat{U}(x) e^{\lambda t}.
$$
where $\hat U(x)$ are the "eigenfunctions" and $\lambda$ the eigenvalues ?

- Is it possible to use such solutions to construct a general solution of the problem ?


### A few definitions


- We first define a vectorial space $S$ defined as the space of functions $U(x)$ 
defined over the range $x\in [a,b]$, with values on $\mathbb{R}$ (or $\mathbb{C}$).
To enunciate mathematical results, we define more precisely $S$ as the space of function 
which are differenciable up to second order (class $C^2$) and have a finite norm ($\int_a^b |U(x)|^2 < \infty$)
(**Sobolev space**). 

-  We can define a scalar product and a norm on $S$ as follows :
$$
< U,V >  = \int_a^b \bar{U}(x) V(x) w(x) dx,
$$ 
$$
|| U || = \left( \int_a^b |U(x)|^2  w(x) dx\right)^{1/2}.
$$
where $w(x)$ is a "weight" function such as $w(x) >0$ (the simplest choice is $w=1$).




## Sturm-Liouville problems.

Definition :

An eigenvalue problem can be set in Sturm-Liouville form if there exist three functions $p(x) > 0$, $w(x) > 0$ and $q(x)$ 
such that :

$$
\lambda w(x) U(x) 
= \frac{\partial}{\partial x} \left( p(x) \frac{\partial U(x)}{\partial x} \right)
+ q(x) U(x)
$$

With boundary conditions at $x=a$ and $x=b$ which can be either Dirichlet, Neumann or Robin.


In this case, it can be proven that the operator $A$ admits an eigenvalue decomposition which generalises 
the case of a symmetric, real matrix. Namely :


- There exist a infinite number of *real* eigenvalues $\lambda_n$  verifying 
$$\lambda_1>\lambda_2 > ... \lambda_n, \quad \, \mathrm{ and } \, 
\lambda_n  \rightarrow - \infty \, \mathrm{ as } \, n \rightarrow \infty
$$

- The corresponding eigenvalues (eigenfunctions) verify an orthogonality condition 

$$
< \hat{U}_n, \hat{U}_m > = \int_a^b \hat{U}_n(x) \hat{U}_m(x) = \delta_{mn} 
$$

- The general solution of the intial value problem can be written as follows :
$$
U(x,t) = \sum_{n=0}^{\infty} c_n \hat{U}_n(x) e^{\lambda_n t} ,\quad \mathrm{ with } \quad c_n = < \hat{U}_n , U_0 > 
$$


Moreover, if $q(x)\le 0$, all the eigenvalues are negative, and the corresponding 
dynamical problem is stable.

### An example

The simplest example is the heat equation for $T(x,t)$  with $x \in [0,L]$ :
$$
\frac{\partial T}{\partial t} =  \frac{\partial^2 T}{ \partial x^2}
$$
with Dirichlet conditions $T(0,t)=0$ and $T(L,t)=0\ \forall t$.

This problem is of Sturm-Liouville form with $w(x)=1,p(x)=1,q(x)=0$. The "eigenmode" solutions 
are classically found as 
$$
\lambda_n = \frac{n^2\pi^2}{L^2} ; \qquad 
\hat T_n (x) = sin (n \pi x/L)
$$



## Non Sturm-Liouville problems

Many problems arising in stability theory are non Sturm-Liouville. There is no theorem predicting the existence of eigenmode solutions. In practise many situations may be encountered :

- The problem may possess "regular eigenfunctions" $\hat U_n(x)$ for values of $\lambda_n$ belonging to the *Discrete spectrum* $S_d$. The number of such regular eigenfunctions may be zero, a finite number or and infinite collection.

- The problem may possess *generalized eigenfunctions* $\tilde U_\lambda(x)$ which are either non-differientiable or non-summable
(hence do not belong to the Sobolev Space), for  values of $\lambda$ belonging to a continuous interval $S_c$ called the *continuous spectrum*.

The solution to the initial-value problem can be expressed as a a sum of the contributions of the regular and generalized eigenfunctions :

$$
U(x,t) = \sum_{\lambda_n \in S_d} c_n \hat{U}_n(x) e^{\lambda_n t} 
+ \int_{\lambda \in S_d} c(\lambda)  \tilde{U}_\lambda(x) e^{\lambda t} d \lambda
$$



# Going further : 
- Charru, chapter 11
- Charru, Iooss & Léger
- Glendinning
- Bender & Orszag (for differential eigenvalue problems)




# Exercices 

## 1. Harmonic oscillator

Consider the damped harmonic oscillator, with the form :
$$m \ddot x = - k x - \mu \dot x
$$

a. Write this problem under the form of a linear dynamical system of order 2, with the form $X = [x,\dot x]$.

b. By mathematical analysis of the problem, identify the five different regimes of solutions (amplified pseudo-pediodic, periodic, damped pseudo-periodic, critical, aperiodic) and characterize their range of existence as function of the damping rate $\mu$.

c. Using the program [PhasePortrait_Linear.m](), draw phase portraits for this problem for the five regimes of solution.

## 2. **Linear differential algebraic systems**

Consider a differential algebraic system defined as 
$$\frac{d}{d t} B X
= A X
$$

a. Compute the eigenvalues and eigenvectors, and give the general solution $X(t)$, for the following matrices :

$$
A = \left[ \begin{array}{cc} 
1 & 2 \\ 3 & 4 
\end{array}
\right] 
\quad 
B = \left[ \begin{array}{cc} 
1 & 0 \\ 0 & 0
\end{array}
\right] 
$$

b. Same questions for the matrices:
$$
A = \left[ \begin{array}{ccc} 
0 & 1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 0
\end{array}
\right] 
\quad 
B = \left[ \begin{array}{ccc} 
1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0
\end{array}
\right] 
$$

c. Same questions for the matrices:
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

## 3. Heat conduction in a cylindrical rod

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
















