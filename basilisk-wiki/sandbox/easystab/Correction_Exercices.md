**Solutions for some of the exercices of the M2R-DET stability course**

### Exercice 1.2


a. Setting $x=x_0$ and $y=y_0$ into the system we get

$$ A + x_0^2 y_0 - (B+1) x_0 = 0 $$

$$ B x_0 - x_0^2 y_0 = 0$$

whose solution is $[x_0;y_0] = [A ; B/A]$.

b. Injecting into the system leads to

$$ dx'/dt = A + (x_0+ x')^2 (y_0 + y') - (B+1) (x_0 + x')$$

$$ dy'/dt = B (x_0 + x') - (x_0+ x')^2 (y_0 + y') $$

Replacing $x_0,y_0$ by their values and retaining only linear terms with respect to $(x',y')$ leads to

$$ dx'/dt = (B-1) x' + A^2 y'$$

$$ dy'/dt = - B x' - A^2 y' $$

Looking for eigenmodes, i.e. $[x';y'] = [\hat{x};\hat{y}] e^{\lambda t}$ leads to the eigenvalue problem

$$
\lambda 
\left[ \begin{array}{c} \hat{x} \\ \hat{y} \end{array} \right] 
= 
\left[ \begin{array}{cc} B-1 & A^2 \\ -B & - A^2 \end{array} \right] 
\left[ \begin{array}{c} \hat{x} \\ \hat{y} \end{array} \right] 
$$
whose characteristic polynomial is

$$
\lambda^2 - (B-1-A^2) \lambda + A^2 = 0
$$

The solution is 
$$
\lambda = \frac{ B - (1+A^2) + \sqrt{ (B-1-A^2)^2 - 4 A^2}}{2 A^2}
$$

The fixed point is stable for $B<B_c$ with $B_c = (1+A^2)$, and is either
a focus or a node depending on the sign of the discriminent.

It becomes unstable for $B>B_c$ and again either a focus or a node.

Close to $B_c$ one has

$$\lambda = \frac{B - (1+A^2)}{2A^2} \pm \frac{1}{A}$$

This is a Hopf bifurcation with characteristic frequency $\omega_0 = 1/A$.

### Exercice 2.2

b. Starting with the equation of motion:

$$
 m L^2 \ddot \theta = - \mu_f \dot \theta - m g L \sin \theta + m L^2 \Omega^2 \sin \theta \cos \theta 
$$

We write $\theta(t) = \theta^* x$ where $x$ is the rescaled variable and $\theta^*$ an amplitude. Similarly we redefine time as $t = t^* \tilde{t}$ where $\tilde{t}$ is the nondimensional time variable and $t^*$ the time scale. . Injecting this variable change in the equation, and developing in series the sin and cosine, leads to:

$$
\frac{m L^2}{{t^*}^2} \theta^* \ddot x = - \frac{\mu_f}{t^*} \theta^* \dot x - m g L \left( \theta^* x - \frac{{\theta^*}^3 x^3}{3} \right)  + m L^2 \Omega^2 \left( \theta^* x - \frac{{\theta^*}^3 x^3}{3} \right) \left( 1 - \frac{ {\theta^*}^2 x^2}{2} \right) + {\mathcal O}(x^5)
$$

Where $\dot x$ and $\ddot x$ are now derivatives with respect to the nondimensional time $\tilde{t}$

This can be regrouped as
$$
\frac{m L^2}{{t^*}^2} \ddot x + \frac{\mu_f}{t^*} \dot x 
= \left( m L^2 \Omega^2 - m g L \right) x - \left(  \frac{5 m L^2 \Omega^2}{6} -\frac{mgL}{2} \right) {\theta^*}^2 x^3 
+ {\mathcal O}(x^5)$$

To get the desired form, we have to set 
$$
t^* = \frac{ m L^2}{\mu_f}
$$
$$
r = \frac{ m L^2}{\mu_f^2} \left( m L^2 \Omega^2 - mgL \right)
$$
$$
{\theta^*} = \left[\frac{ m L^2}{\mu_f^2} \left(  \frac{5 m L^2 \Omega^2}{6} -\frac{mgL}{2} \right) \right]^{-1/2}
$$




### Exercice 6.0

We use the same convention as in Charru : $\bar{U}(y) = U_1$ for $y<0$ and $\bar{U}(y) = U_2$ for $y>0$.

#### A. Demonstration using velocity potential

*This first demonstation uses the velocity potential, and the approach is very similar to the treatment of the Rayleigh-Taylor instability done in chapter 4.*

The assumption of *potential flow* can be justified as follows.
For both $y>0$ and $y<0$, the initial flow (i.e. the base flow) is irrotational. Neglecting diffusion effets, the flow remains potential for all times (Kelvin theorem). Hence we can introduce a velocity potential satisfying $\Delta \phi = 0$.

Assuming modal dependence with form $e^{i k x} e^{-i \omega t}$, the solution decaying at $|y| \rightarrow \infty$ is

$$
\phi(x,y,t) = \left\{ \begin{array}{ll} 
A e^{-ky} e^{i k x} e^{-i \omega t} & \qquad (y>0), \\ 
B e^{+ky} e^{i k x} e^{-i \omega t} & \qquad (y<0). \\ 
\end{array}
\right.
$$


Moreover we assume a displacement of the shear layer with the form:
$$
y = \eta(x,t) = C e^{i k x} e^{-i \omega t}.
$$

*This starting point is the same as for the Rayleigh-Taylor instability in chapter 4, but the pressure field and the boundary conditions are different.*

The pressure field is deduced from the potential by the unsteady Bernoulli equation :
$$ p' + \rho \frac{\partial \phi}{\partial t} + \rho | {\bf u} ^2|/2 = Cte.$$
Considering $|{\bf u} ^2| = (U_{1,2} {\bf e_x} + \nabla \phi)^2 \approx U_{1,2}^2 + 2 U_{1,2} \partial \phi/ \partial x$, we are lead to:

$$
p'(x,y,t) = \left\{ \begin{array}{ll} 
- \rho ( i k U_1 - i \omega) B e^{+ky} e^{i k x} e^{-i \omega t} & \qquad (y<0). \\ 
- \rho ( i k U_2 - i \omega) A e^{-ky} e^{i k x} e^{-i \omega t} & \qquad (y>0), \\ 
\end{array}
\right.
$$

The dynamical boundary condition  $p(x,y=0^-,t)= p(x,y=0^+,t)$ leads to:
$$
(k U_1 - \omega ) B - (k U_2 - \omega ) A = 0 \qquad \mathrm{(Eq. 1)}
$$

The kinematical boundary condition is 
$d\eta/dt = v'(x,y=0^-,t) =  v'(x,y=0^+,t)$ and leads to:
$$
i (k U_1 - \omega ) C = + k  B; \qquad i (k U_2 - \omega ) C = - k A.
$$

Hence 

$$
\frac{k}{k U_1 - \omega } B + \frac{k}{k U_2 - \omega } A = 0 \qquad \mathrm{(Eq. 2)}
$$
Taking the determinent of the system (Eq.1,Eq.2) for $(A,B)$ leads to:

$$
(k U_1-\omega)^2 + (k U_2 - \omega)^2 =0
$$

Whose solution is $$
\omega = k \left( \frac{U_1+U_2}{2} \pm i \frac{U_1-U_2}{2} \right)
$$

#### B. Demonstration using streamfunction

*This demonstration is more in the line as done in Charru, and starts from the Rayleigh equation in terms of the streamfunction. The use of a streamfunction relies on the assumption of incompressible flow which is easier to justify than the hypothesis of irrotationnal flow. *

For both $y>0$ and $y<0$, the Rayleigh equation reduces to 
$\Delta \hat{\psi} = (\partial_{yy} - k^2) \hat{\psi} = 0$.

Assuming modal dependence with form $e^{i k (x-ct) }$, 
the solution decaying at $|y| \rightarrow \infty$ is

$$
\hat{\psi}(y) = \left\{ \begin{array}{ll} 
A_2 e^{-ky}& \qquad (y>0), \\ 
B_1 e^{+ky}  & \qquad (y<0). \\ 
\end{array}
\right.
$$

Moreover we assume a displacement of the shear layer with the form:
$$
y = \eta(x,t) = C e^{i k (x-ct)}.
$$

To obtain an expression for the pressure field we go back to the linearised equation for $u'$ : 
$$i k \rho (\bar{U}-c) \hat{u} + \bar{U}'(y) \hat{v} = - i k \hat{p}$$ 

Considering that ${\bar{U}'(y)}=0$ for $y \ne 0$, and 
$\hat{u} = \partial \hat{\psi}/\partial y$,
this leads to:

$$
\hat{p}(y) = \left\{ \begin{array}{ll} 
- \rho ( U_2 - c) A_2 \left[ -k  e^{-ky} \right]  & \qquad (y>0), \\ 
- \rho ( U_1 - c) B_1 \left[ +k e^{+ky} \right]  & \qquad (y<0). \\ 
\end{array}
\right.
$$

The dynamical boundary condition  $p'(x,y=0^+,t)= p'(x,y=0^-,t)$ leads to:
$$
(U_2 - c ) A_2 + ( U_1 - c ) B_1 = 0 \qquad \mathrm{(Eq. 1)}
$$

The kinematical boundary condition is 
$d\eta/dt = v(x,y=0^+,t) =  v(x,y=0^-,t)$, and leads to:
$$
i k (U_2 - c ) C = -i k  A_2; \qquad i k (U_1 - c ) C = - i k B_1.
$$

Hence 

$$
\frac{1}{U_2 -c  } A_2 - \frac{1}{U_1 - c } B_1 = 0 \qquad \mathrm{(Eq. 2)}
$$
Taking the determinent of the system (Eq.1,Eq.2) for $(A,B)$ leads to:

$$
(U_1-c)^2 + (U_2 - c)^2 =0
$$

Whose solution is 
$$
c = \left( \frac{U_1+U_2}{2} \pm i \frac{U_1-U_2}{2} \right)
$$




### Exercice 6.1
 
 
We work in cylindrical coordinates and introduce the state-vector $q = [u,v,p]$ where $u$ and $v$ are the radial and azimuthal components of the velocity.

We expand the flow as follows:
 $$ q = q_0 + \hat{q} e^{i m \theta} e^{-i \omega t}$$ 
 where 
 $q_0 = [0,\bar{V}(r),\bar{P}(r)]$ is the "base flow" and $\hat{q} = [\hat{u},\hat{v},\hat{p}]$ is the perturbation in eigenmode form.

Injecting in the Euler equations leads to:
$$
- i \omega \hat{u} = - i \frac{m}{r} \hat{u} +  \frac{2 \bar{V}(r)}{r} \hat{v} -\partial_r \hat{p} +\underbrace{ \frac{\bar{V}(r)^2}{r} - \partial_r \bar{P}}_{=0}
 $$
 
$$
- i \omega \hat{v} =  - i \frac{m}{r} \hat{v} - \left( \frac{\bar{V}(r)}{r} + \partial_r \bar{V}(r) \right) \hat{u} -i \frac{m}{r} \hat{p}
$$

$$
0 = \left(\frac{\partial}{\partial r}  + \frac{1}{r}\right) \hat{u} +i \frac{m}{r} \hat{v}
$$
 
 (the underbraced term in the first expression cancels because the base-flow is an equilibrium solution of the Euler equations)
 
 Introducing the streamfunction and combining the two first expression leads to the desired result:
 
 $$
(m \Omega(r) - \omega) 
\left( \partial_r^2 + r^{-1} \partial_r -m^2/r^2 \right) 
\hat{\psi} - m/r \partial_r \Xi(r) \hat{\psi} = 0.
$$

To demonstrate the criterion, we multiply by $r \frac{\hat{\psi}^*}{(m \Omega(r) - \omega)}$ and integrate over $[r_1,r_2]$. 

Tip : $(\partial_r^2 + r^{-1} \partial_r)  = \frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{\partial}{\partial r} \right)$.

This leads to:
$$
- \int_{r_1}^{r^2}  \left[ \partial_r \hat \psi \right|^2 \, r dr + \underbrace{{\left[ r \hat \psi^* \partial_r \psi \right]}_{r_1}^{r^2}}_{=0} - m^2 \int_{r_1}^{r^2} \left| \hat \psi \right|^2 r^{-1} dr = m \int_{r_1}^{r^2} \frac{ \left| \hat \psi \right|^2 [\partial_r \Xi(r)] }{(m \Omega(r)-\omega)} \, dr
 $$
 
 Note that the boundary term arising from the integration by parts cancels because if the boundaries at $r_1$ and $r_2$ are walls the boundary condition is $\hat{u} = i m /r \hat{\psi} = 0$, hence $\hat{\psi}(r_1) = \hat{\psi}(r_2) = 0$.
  
 The imaginary part leads to:
 

$$
0 = m \omega_i \int_{r_1}^{r^2} \frac{ \left| \hat \psi \right|^2 [\partial_r \Xi(r)] }{(m \Omega(r)-\omega_r)^2 + \omega_i^2} \, dr
$$

Which allows to conclude that if there exists and unstable mode ($\omega_i>0$), then the term $[\partial_r \Xi(r)]$ necessarily changes sign along the interval. 
  
## Exercice 6.2


 We look for solution of the Rayleigh equation (which reduces to $\Delta \hat{\phi} = 0$ in each layer) as follows:
 
$$
\begin{array}{rcll}
 \hat{\phi}_1 &=& B_1 e^{ky}, \qquad \qquad & y<-\delta
 \\
 \hat{\phi}_0 &=& B_0 e^{ky} + A_0 e^{-ky}, & -\delta<y<\delta
\\
\hat{\phi}_2 &=& A_2 e^{-ky}, \qquad \qquad & y>\delta
\end{array}
$$

Writing $K = k \delta$, the kinematic and dynamic matching conditions at $y= \pm \delta$ lead to a matricial system
 
$$[A] \left[ \begin{array}{c} B_1 \\ A_0 \\ B_0 \\ A_2
\end{array}\right]$$
 where
 $$
 [A] = \left[\begin{array}{cccc} 
 exp(-K) & -exp(-K) & -exp(K) & 0 \\
     K (U_2-c) exp(-K)& (-K (U_2-c)+\Delta U) exp(-K)
      & (K(U_2-c)+\Delta U ) exp(K)& 0 \\
       0 &  (-K(U_1-c)+\Delta U) exp(K) & (K(U_1-c)+\Delta U) exp(-K) &  -K(U_1-c) exp(-K) \\
       0 &  -exp(K) & -exp(-K) & exp(-K)
\end{array}
\right]
$$

The determinant of this matrix leads to the expected result.

We can then write the dispersion relation as
$$
(U_m - c)^2 - (\Delta U)^2 F(K) = 0
$$
where 
$$
F(K) = \frac{(2 K -1)^2- e^{-4 K}}{4 K^2}
$$ 

The solution is then $c = U_m \pm \Delta U \sqrt{F(K)}$.

We can check that :

- $F(K) \approx -1$ as $K \approx 0$, so in this limit the dispersion relation is equivalent to that of the zero-thickness shear layer,

 
-  $F(K) < 0$ for $0<K<0.6932$, so in this range of long-wavelength the problem is unstable (two complex solutions for $c$ and one of them veries $c_i >0$).

- $F(K) >0$ for $K<0.6932$, so in this range of short-wavelength the problem is stable (two real solutions for $c$ corresponding to two neutral waves).


## Exercice 8.1

**1. Case $\epsilon = 0$.**

In this case the system is $dx/dt = y$, $dy/dt=0$. The solution, considering initial conditions $x_0,y_0$, is
$$
x(t) = x_0+ y_0 t , \quad y(t) = y_0.
$$

The enegy growth $G(t) = \frac{||X(t)||^2}{||X_0||}$ behaves, as $t\rightarrow \infty$, as
$$
G(t) \approx \frac{y_0^2}{x_0^2+y_0^2} t^2
$$

One observes that in this case displays {\em algebraic growth}. 
According to the general definitions of chapter 1, the system is unstable nut  {\em not exponentially unstable}. ( the two eigenvalues of the problem are zero).

**2. Eigenvalues/eigenvectors **

The eigenvalue-eigenvectors of the matrix 
$A =  \left[\begin{array}{cc} -\epsilon &1  \\ 0 & - 2 \epsilon \end{array}\right]$
are as follows :

$$
\lambda_1 = - \epsilon ; \quad 
\hat{X}_1 = \left[\begin{array}{c} 1  \\ 0   \end{array}\right]; 
\quad 
\lambda_2 = - 2 \epsilon ; \quad 
\hat{X}_2 = \frac{1}{\sqrt{1+\epsilon^2}}\left[\begin{array}{c} -1  \\ \epsilon  \end{array}\right].
$$

According to the general definitions of chapter 1, the system is exponentially stable (two negative eigenvalue). Accordingly, for $t\rightarrow \infty$ the energy growth $G(t)$ ultimately tentds to zero. However we will show that it may allow very large transient growth for intermediate values of $t$.

**3. Initial value problem **

The general solution can be expanded upon the basis of eigenvectors, under the form 
$$
X(t) = c_1 \hat{X}_1 e^{\lambda_1 t} + c_2 \hat{X}_2 e^{\lambda_2 t}.
$$
There are 3 possible methods to compute the coefficients $c1,c2$:

**3.a First method : direct determination **

The condition $X(t) = X_0$ leads, when separating $x$ and $y$ components, to two equations :
$$
c_1 - \frac{c_2}{\sqrt{1+\epsilon^2}} = x_0, \qquad \frac{\epsilon c_2}{\sqrt{1+\epsilon^2}} = y_0.
$$
The resolutions yieds $c_1$ and $c_2$ as functions of $x_0$ and $y_0$.

NB this method is easy for a dimension-2 problem, but not applicable to larger problems. We demonstrate two other methods which are directly generalisable to larger problems.

**3.b Second method : using adjoint basis **

The dual basis is defined by the property $\hat X_i^\dag \cdot  \hat X_j = \delta_{ij}$ (see the reference document [LinearSystems.md](http://basilisk.fr/sandbox/easystab/LinearSystems.md).

This leads to :
$$
\hat{X}_1^\dag = \frac{1}{\epsilon}\left[\begin{array}{c} \epsilon   \\ 1   \end{array}\right]; 
\quad 
\hat{X}_2^\dag = \frac{\sqrt{1+\epsilon^2}}{\epsilon} \left[\begin{array}{c} 0  \\ 1  \end{array}\right].
$$


The initial condition $X_0 = [x_0 ; y_0]$ can be expanded on the base of eigenvectors as $X_0 = c_1 \hat{X}_1 + c_2 \hat{X}_2$.
The coefficients $c1,c_2$ can be determined thanks to the dual basis:

$$c_1 = (\hat{X}_1^\dag \cdot X_0) = x_0 - \frac{y_0}{\epsilon} ; \qquad 
c_2 = (\hat{X}_2^\dag \cdot X_0) = \frac{\sqrt{1+\epsilon^2}}{\epsilon}y_0$$

This leads to the desired expression.


**3.c Thirds method : matrix exponentials**

The general solution for initial condition $X_0$ can be written as 
$$
X(t) = e^{A t} X_0
$$

To evaluate the matrix exponential we first diagonalize it by writing 
$$ A = U D U^{-1}$$

Here $D$ is a diagonal matrices with the eigenvalues on the diagonal, and $U$ a square matrix contaioning the eigenvectors in column, namely $U = [\hat{X}_1, \hat{X}_2]$.


$$
D = 
 \left[\begin{array}{cc} - \epsilon & 0  \\ 0 & - 2 \epsilon \end{array}\right]
$$


$$
U = 
 \left[\begin{array}{cc} 1 & \frac{-1}{\sqrt{1+\epsilon^2}}  \\ 0 & \frac{\epsilon}{\sqrt{1+\epsilon^2}}  \end{array}\right]
$$


Inverting this matrix leads to 


$$
U^{-1} = 
 \left[\begin{array}{cc} 1 & \frac{1}{\epsilon}  \\ 0 & \frac{\sqrt{1+\epsilon^2}}{\epsilon}  \end{array}\right]
$$

Then we can use the identity $e^{A t} = U e^{ D t } U^{-1}$ to obtain the final result.


Note that the matrix $U^{-1}$ is actually composed by the adjoint eigenvectors $\hat{X}_1^\dag$ and $\hat{X}_2^\dag$ written in lines. This property allows to verify that the two methods are equivalent.

**4. Initial growth**

For $t \ll \epsilon^{-1}$, a taylor development leads to $e^{-\epsilon t} \approx 1-\epsilon t$, $e^{-2 \epsilon t} \approx 1-2 \epsilon t$. One can verify that the solution becomes equivalent to the case considered in question 1. Thus, the transient growth is initially algebraic.

**5. optimal perturbation** 

An initial condition with norm 1 can be written $[x_0,y_0] =[\cos \theta, \sin \theta]$ without loss of generality.


The corresponding energy growth is 
$$
G(t) = \left[ \cos \theta e^{-\epsilon t} + \frac{\sin \theta}{\epsilon} \left( e^{-\epsilon t} - e^{-2 \epsilon t} \right) \right]^2 
+ ( \sin \theta e^{-2 \epsilon t} )^2
$$
Keeping only the leading order term considering $\epsilon \ll 1$ :
$$ G(t) \approx 
\left(\frac{\sin \theta}{\epsilon}\right)^2 \left( e^{-\epsilon t} - e^{-2 \epsilon t} \right)^2
$$

One can  deduce that the initial condition leading to the largest growth corresponds to the direction $\sin \theta = 1$ (i.e. $[x_0,y_0]  = [0,1]$, which physicvally corresponds to an initial condition consisting of vortices only).

**6. Optimal growth** 

Starting from the expression of $G(t)$ written above (in the case $\epsilon \ll 1$), the time corresponding to maximum can be obtained from 
$$\partial G/\partial t =
\left(\frac{\sin \theta}{\epsilon}\right)^2  2 \left( e^{-\epsilon t} - e^{-2 \epsilon t} \right) \left( -\epsilon e^{-\epsilon t} + 2 \epsilon e^{-2 \epsilon t} \right) 
$$
Factorizing a number of terms, the condition $\partial G/\partial t=0$ is equivalent to $\left( -1 + 2 e^{- \epsilon t} \right) =0$.
The solution is $t=t^{opt} = \frac{\log(2)}{\epsilon}$. 
The maximum growth is then 

$$
G^{opt} = max_{(t,\theta)} G(t) = G(t^{opt},\theta = \pi/2) = \frac{1}{16 \epsilon^2}.
$$



## Exercice 9.1

1. *Temporal stability analysis* (solving for $\omega$ as function of $k$
leads to 
$$ \omega \equiv \omega(k) = k U_m + i \Delta U k (1-k)$$

We can see that : a/ the phase speed is $c_r = \omega_r/k = U_m$ is constant ;
b/ the amplification rate $\omega_i$ is positive for $0<k<1$ and maximum for $k= 1/2$ (with value $\Delta U/4$).

This looks much like the properties of the Kelvin-Helmholtz instability of a shear layer with tanh or piecewise velocity profiles (see chapter 6) ; hence this simple model can be used to study qualitatively this case.

2. Spatio-temporal stability analysis consists of finding the $k_0$ and $\omega_0$ which are solution of the dispersion relation and the statiuonary phase condition $\partial \omega / \partial k = 0$.

Using the model one has
$$\frac{\partial \omega}{\partial k} = U_m + \Delta U i (1- 2k) = 0$$
Solving for $k$ leads to 
$$k \equiv k_0 = \frac{U_m - i \Delta U}{2 \Delta U}$$
injecting back in the dispersion relation leads to the expected result:
$$
\omega_0 = \omega(k_0) \frac{\Delta U}{2}  + \frac{ i ( (\Delta U)^2-U_m^2) }{4 \Delta U}    $$

3/ From the previous expression we can see that if $R =\Delta U/U_m < 1$ then the growth rate $\omega_{0,i}$ of the stationnary-phase perturbation is negative: *convective instability*.
On the other hand, if $R =\Delta U/U_m > 1$ then the growth rate $\omega_{0,i}$ of the stationnary-phase perturbation is positive : *absolute instability*.

4/ Direct resolution of the second-order polynomial for $k$ leads to:

$k(\omega) = \frac{(U_m + i \Delta U) \pm \sqrt{  (U_m + i \Delta U)^2 - 4 i \Delta U \omega}}{ 2 i \Delta U}$

However this expression is difficult to interpret because of the determination of the complex square root. Plotting using a computer shows that if $R<1$ one gets two well-separated branches $k^+(\omega)$ and $k^-(\omega)$ (relevant to study spatial instability in the convective case) while if $R>1$ the two branches exchange their identities and can no longer be associated to spatial stability properties.



[//]: # (It’s a little bizarre, but it works with MacDown and Pandoc.)


