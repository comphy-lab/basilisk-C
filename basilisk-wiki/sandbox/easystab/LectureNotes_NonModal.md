
**3D and transient growth mechanism in shear flows (chap. 8)**

We have seen that for important classes of flows belonging to class C, 
the classical theory 
 (*linear* *modal*  *2D*) cannot explain the onset of turbulence which is observed for $Re \gtrsim 2000$.
 
Even for boundary layers flows belonging to class B, experiments show that transition to turbulence is not always associated to TS waves. 

In some cases instead of TS waves with a 2D structure (transverse vortices) we observe *streamwise* structures, of two kinds:

* "Streaks" (modulation of the axial velocity $u$) 

* Axial Vortices (associated to transverse velocities $v$ and $w$)


[Illustration : streamwise vortices in boundary layers (see figure 10  of this paper)](https://engineering.jhu.edu/zaki/wp-content/uploads/2018/12/Lee_CF_2018.pdf)

# Modal stability for 3D perturbations.


## General equations 

Let us reconsider the case of *Oblique waves* with modal behavior 
$[u',v',w',p'] = \hat{q} e^{ i k x + i \beta z - i \omega t}$ with $\beta \ne 0$. 

The governing equations are derived in [this introductory document  ](LectureNotes_StabilityOfParallelFlows.md) and can be written in matricial form as follows:

$- i \omega {\mathcal B}_{3D} \hat{q} = {\mathcal A}_{3D}(k,\beta,Re) \hat{q}$ 


## Squire theorem

[Exercice associated to lecture 7](LectureNotes_Viscous.md#exercices)


Starting from the general equations, one can introduce the following change of variables (squire transformation):
$$
\tilde{k} = \sqrt{k^2 + \beta^2}, \quad \tilde{\omega} = \frac{\tilde k}{k} \omega, \quad \tilde{Re} =  \frac{k}{\tilde k} Re, \quad  
$$
$$
\tilde{\hat{u}} = (k/\tilde{k})\,  \hat{u} + (\beta/\tilde{k})\,  \hat{w} , \quad \tilde{\hat{u}} = \hat{v}, \quad \tilde{\hat{p}} = (k/\tilde{k}) \, \hat{p}
$$

Then the problem can be reduded to an equivalent 2D problem :

 $- i \tilde{\omega} {\mathcal B}_{2D} \tilde{\hat{q}} = {\mathcal A}_{2D}(\tilde{k},\tilde{Re}) \tilde{\hat{q}}$

with $\hat{q} = [ \tilde{\hat{u}},\tilde{\hat{v}},\tilde{\hat{p}}]$


An important consequence is that if for some $(k,\beta,Re)$ there exists an 3D unstable eigenmode with growth rate $\omega_i$, there necessarily exists a  2D eigenmode for a *smaller value of Re* (given by $\tilde{Re}$) which *is more unstable* (with growth rate given by $\tilde{\omega_i}$).

A corollary is that the critical value $Re_c$ above which instability exists is smaller for 2D than for 3D.

Another corollary is that perturbations with $k=0$, $\beta \ne 0$ are always stable.


Hence 3D modal instability is a wrong direction to resolve the puzzle !


# Non-modal theory 

## Equations
Let us try a second idea : abandon the modal expansion and go back to temporal equations. Let us consider for simplicity the case $k=0$ corresponding to *streamwise structures*. 
Consider the Couette flow for which $\bar{U}(u) = U_0 y/d$ (for $y\in [0,d]$)

The temporal equations written in [the introductory document  ](LectureNotes_StabilityOfParallelFlows.md) reduce to:

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
Re^{-1} ( \partial_y^2 - \beta^2) & - \partial_y \bar{U} & 0 & 0 \\ 
0 &  Re^{-1} ( \partial_y^2 - \beta^2) & 0 & - \partial_y \\
0 & 0 & Re^{-1} ( \partial_y^2 - \beta^2) & - i \beta \\
0  & \partial_y & i \beta & 0 
\end{array} 
\right] 
$$


 
## A simple model

Let us introduce a simple model mimmicking the full set of equations :

$$ 
\dot x = - Re^{-1} x  + y 
$$
$$ 
\dot y = - 2 Re^{-1} y
$$

In this equation $x$ models the amplitude of an array of "streaks"  (associated to velocity component $u$) and $y$ models the amplitude of an array of vortices.

The study of this equation is left as an [exercice](
http://basilisk.fr/sandbox/easystab/LectureNotes_NonModal.md#exercice-1-transient-growth-in-a-nonnormal-linear-system.). One can show that :

* Modal analysis shows that the system is linearly stable, with two eigenvalues given by $\lambda_1 = -1/R$, $\lambda_2 = -2/R$.
* However, strong transient growth are possible ! 

To characterize the transient growth, one can define the optimal growth as 
$$
G^{opt}(t) = max_{||X(0)||=1 } \frac{||X(t)||}{||X(0)||} \quad \mathrm{ with } \,  X(t) = [x(t);y(t)].
$$
and subsequently define the maximum growth as 
$$
G^{max} = max_{t\ge 0} G^{opt}(t) 
$$

One can demonstrate that  $G^{max} = Re^2/16$ and that this maximum energy growth corresponds to an initial condition $[x(t);y(0)] \approx [0,1]$ (close to pure vortices) and to an optimal time $t_{max} = Re \ln(2)$. The state observed at $t_{max}$ has $||x(t_{max}|| \gg ||y(t_{max}||$ (close to pure streaks).


## Extension to nonlinear model

The very large transient growth suggests that nonlinearities could become active and lead to a "bypass transition".
To investigate this, Trefethen (1993) introduced the following model:

$$ 
\dot x = -(1/R) x  + y -\sqrt{x^2+y^2}  y;
$$
$$ 
\dot y = -(2/R) y + \sqrt{x^2+y^2} x;
$$

Numerical resolution of this problem shows that although infinitely small perturbations always end up by decaying, in accordance with the linear theory, there exist some small but finite initial condition (with norm $\sqrt{x_0^2+y_0^2} = \mathcal{O} (Re^{-1})$) which do not go back to zero but end up to another attractor (corresponding to a fixed point $[x,y] \approx [0,1]$).

## The "bypass transition" scenario
 
This model is at the origin of a "bypass transition" scenario explaining the transition to turbulence. This scenario postulates that when starting from the laminar state (represented by the fixed point $[x,y] = [0,0]$ in the model), small-amplitude perturbations corresponding to vortices (represented by $y$ in the model) can lead to the transient growth of streaks (represented by $x$ in the model), and that once these streaks become of sufficiently large amplitudes, nonlinearities take over, eventually leading to the fully turbulent state (represented by the second fixed point $[x,y] \approx [0,1]$ in the model).


## Optimal amplification of a slightly unstable mode
  
  For oblique waves ($\beta \ne 0$, $k\ne 0$) the Tollmien-schlishting (modal) and transient growth (non-modal) mechanism can coexist.
 To investigate such a sitution consider the model system:
 
 $$ 
\dot x = + Re^{-1} x  + y 
$$
$$ 
\dot y = - 2 Re^{-1} y
$$

This system displays a (slightly unstable) eigenmode with $\lambda_1 = + Re^{-1}$ and a stable mode with $\lambda_2 = - 2 Re^{-1}$. In this case one can show (exercice) that the optimal growth is $G^{opt}(t) = \frac{Re^2}{9} e^{2 Re^{-1} t}$. Compared to an initial condition proportional to the most unstable mode (which simply gives an energy growth $e^{2 Re^{-1} t}$), the optimal perturbation thus substantially boosts the growth, possibly leading to a faster transition to turbulence.


# Exercices for chapter 8
 
### **Exercice 1** Transient growth in a nonnormal linear system.

[Solution](/sandbox/easystab/Correction_Exercices.md#exercice-8.1) 

*Warning : this exercice is close to the case considered by Charru (p. 31-32) but there is an error in Eq. 1.43. Will you find it ?*

*This exercice is also treated (considering a slightly different system) in Rieutort, p. 228 ; Schmid & Henningson, p. 522.*



We consider the following two-dimensional system :

$$\frac{dx}{dt} = -\epsilon x + y $$

$$\frac{dy}{dt} = - 2 \epsilon y$$

Or in matricial form :

$$
\frac{dX}{dt} = \left[\begin{array}{cc} -\epsilon &1  \\ 0 & - 2 \epsilon \end{array}\right]  X.
$$

<span style="color:red">
The purpose of the exercice is to study the  energy growth defined by $G(t) = \frac{||X(t)||}{||X_0||}$.

1. Give the solution of this system for $\epsilon = 0$. What can be observed ?


2. Determine the eigenvalues ($\lambda_1, \lambda_2)$ and normalised eigenvectors $(\hat X_1, \hat X_2$) of the system.

3. Show that the general solution of this system witg initial condition $[x_0,y_0]$ is as follows:

$$
x(t) = x_0 e^{-\epsilon t} + \frac{y_0}{\epsilon} \left( e^{-\epsilon t} - e^{-2 \epsilon t} \right)
$$
$$
y(t) = y_0 e^{-2 \epsilon t}
$$

Tip : you may introduce the *adjoint* (or dual) basis $(\hat X_1^\dag, \hat X_2^\dag$)
defined by the property $\hat X_i^\dag \cdot  \hat X_j = \delta_{ij}$.

4. Study the behaviour of this expression when $t \ll \epsilon^{-1}$.  Compare with the solution for $\epsilon = 0$ (question 1).

5.  Among the initial conditions with $||X_0|| = \sqrt{x_0^2+y_0^2}=1$, 
which one leads to the largest value of the energy growth factor
$G(t) = \frac{ x(t)^2+y(t)^2}{x_0^2+ y_0^2}$ ?
(this initial condition is called the *optimal perturbation*).

6. Show that the energy growth factor is maximum for $t = \frac{\log(2)}{\epsilon}$, and leads to an optimal growth $G^{max} = \frac{1}{16 \epsilon^2}$.



<!--
5. Among the initial conditions with $||X_0|| = \sqrt{x_0^2+y_0^2}=1$, 
which one leads to the maximum 
(this initial condition is called the *optimal perturbation*).


6. Show that the energy growth associated to the optimal perturbation is given by

$$G = \frac{x(t)^2+y(t)^2}{x_0^2+y_0^2} = e^{-2 \epsilon t} 
+ \left(\frac{e^{-\epsilon t} - e^{-2 \epsilon t}}{\epsilon} \right)^2
$$


7. Show that when $\epsilon \ll 1$ the energy growth reaches a maximum for $t_{opt} \approx ln/\epsilon$, with $G_{max} = G(t_{opt}) \approx 1/(16 \epsilon^2)$.


8. Show that for $t \ll t_{opt}$ the energy growth is algebraic. Interpret this by comparing with the solution of the linear system with $\epsilon = 0$.
-->

### **Exercice 2** bypass transition in a stable, nonnormal system.

*Sources : Schmid & Henningson, p. 525; Trefethen et al. , Science, 1993.*


We consider the following two-dimensional system :


$$
\frac{dX}{dt} = \left[\begin{array}{cc} -\epsilon & 1 \\ 0 & - 2 \epsilon \end{array}\right]  X 
+ ||X|| \left[\begin{array}{cc} 0 & -1 \\ 1 & 0 \end{array}\right] X
$$

Using the program [PhasePortrait_NonLinear.m](), study the phase portrait of this problem.

Show that the combined effect of transient growth and nonlinearities can lead to a "bypass" transition to a non-zero solution,
even for very small initial conditions.

### **Exercice 3:** Optimal excitation of a slightly unstable mode.

We consider the following two-dimensional system :

$$\frac{dx}{dt} = +\epsilon x + y $$

$$\frac{dy}{dt} = - 2 \epsilon y$$

Or in matricial form :

$$
\frac{dX}{dt} = \left[\begin{array}{cc} \epsilon &1  \\ 0 & - 2 \epsilon \end{array}\right]  X.
$$


1. Determine the eigenvalues ($\lambda_1, \lambda_2)$ and normalised eigenvectors $(\hat X_1, \hat X_2$) of the system.

2. Show that the general solution of this system witg initial condition $[x_0,y_0]$ is as follows:

$$
x(t) = x_0 e^{\epsilon t} + \frac{y_0}{3 \epsilon} \left( e^{\epsilon t} - e^{-2 \epsilon t} \right)
$$
$$
y(t) = y_0 e^{-2 \epsilon t}
$$

3. Show that the optimal energy growth factor is given by 
$$G^{opt}(t) = max_{||X_0||\ne 0} \frac{||X(t)||^2}{||X_0||^2} 
\approx \frac{1}{9 \epsilon^2} e^{2 \epsilon t} $$.
Characterize the corresponding initial condition $X_0$ (called the optimal perturbation).

4. Application : consider $\epsilon = 10^{-3}$ and an initial condition $||X_0|| = 10^{-6}$. Estimate the time scale $t_s$ leading to $||X_0|| \approx 1$ (i.e. the time scale at which nonlinearity is expected to dominate) consindering that the initial condition is $(a)$ the unstable eigenmode and ($b$) the optimal perturbation. 

 [Back to main page](http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md)
