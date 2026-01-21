/**
# One-dimensional Stefan solidification (benchmark)


## Problem description
This case provides a minimal reproducible benchmark for validating the improved phase-field model against the classical one-dimensional Stefan solidification problem. The ice front moves purely due to thermal diffusion, and the effect of convection is neglected.

The solid--liquid interface is located at $x=x_s(t)$ with the initial position $x_s(0)=x_0$.

## Governing equation (dimensionless temperature)
In each phase, the dimensionless temperature satisfies the diffusion equation
$$
\frac{\partial \theta}{\partial t} = \frac{\partial^2 \theta}{\partial x^2}.
$$

## Domain and boundary conditions
The semi-infinite domain is truncated to a finite interval $x\in[0,2]$.
Boundary conditions are prescribed as
left boundary: $\theta(0,t)=0$,
right boundary: $\theta(2,t)=\theta_m$,
with $\theta_m=1$.

At the sharp interface,
 $\theta(x_s(t),t)=\theta_m$,
the Stefan condition holds
$$
\mathrm{St}\,\frac{dx_s}{dt} =
\left.\frac{\partial \theta}{\partial x}\right|_{x=x_s^{-}}
-
\left.\frac{\partial \theta}{\partial x}\right|_{x=x_s^{+}}.
$$

## Analytical solution for validation
The analytical solution used for validation is
$$
x_s(t)=2\Lambda\sqrt{t+t_0},\qquad
\theta(x,t)=
\begin{cases}
\displaystyle \frac{\mathrm{erf}\!\left(\frac{\Lambda x}{x_s(t)}\right)}{\mathrm{erf}(\Lambda)}, & x<x_s(t),\\[8pt]
\theta_m, & x\ge x_s(t),
\end{cases}
$$
where $\Lambda$ satisfies
$$
\sqrt{\pi}\Lambda e^{\Lambda^2}\mathrm{erf}(\Lambda)=\frac{1}{\mathrm{St}}.
$$

In this setup, the initial interface location is $x_s(0)=x_0=0.2$, corresponding to a time shift
$$
t_0=\left(\frac{x_0}{2\Lambda}\right)^2.
$$
The initial temperature field is set consistently with the analytical solution
*/

//#include "axi.h" 
#include "CFDing_pp.h" 
#include "CFDing_centered.h"         
#include "CFDing_two-phase.h" 

T[left] = dirichlet(0.);
T[right] = dirichlet(1.);


int main() {	
  size(L_domain);
  init_grid(num_initgrid);
  Cn = 0.5*L0/numofcell ;
  D_AC = 1.0;	
  DT = DT_Max;
  TOLERANCE = 1e-4; 
  run();
}

// the initialization of u, phi,T;
#include "CFDing_init.h"

//tecplot output
#include "CFDing_output.h"

// AC equation is advanced in "CFDing_AC2phase.h"; 
#include "CFDing_AC2phase.h"


// temperature is advanced in "CFDing_temperature.h";
#include "CFDing_temperature.h"