/**

# Liquid film flowing along the wall

#### Theory :

We are trying to obtain the flowrate of a liquid film flowing down a wall under the action of the gravity.

With the dimensional analysis, we can start by saying that the flowrate  $Q$ should depend on

$$ Q = f(\rho, \mu, h, \gamma, \rho_{air}, \mu_{air}, g) $$

Considering that the surronding air has no influence on the dynamic of the fluid, and taking into account that the interface is flat, we can then say that the surface tension $\gamma$ will have no influence, we can then simplify the expression as follows :

$$ Q = f(\rho, \mu, h, g) $$

By the $\pi$-theorem, we can deduce the following expression :

$$ Q = h^{\frac{3}{2}}g^{\frac{1}{2}}F(\frac{1}{Ar}) $$

avec $Ar = \frac{\rho\sqrt{gh}h}{\mu}$

## Hydrodynamic analysis : 

Using the natural scale $h$, $g$ and $\rho$ of the problem, we write the nondimensionnal variables :

$$ 
x= h\,\bar x, \\[.5em]
y = h\,\bar y, \\[.5em]
u = \sqrt{gh}\, \bar u\\[.5em]
v = \sqrt{gh}\, \bar v\\[.5em]
p = \rho g h\, \bar p 
$$


## Governing equation


Upon using the steadiness of the flow and the invariance in the $y$-direction, the momentum equation reads :

$$
\bar u\frac{\partial \bar u}{\partial \bar x} = - \frac{\partial \bar p}{\partial \bar x} + \frac{1}{Ar} (\frac{\partial^2 \bar u}{\partial \bar x^2} + \frac{\partial^2 \bar u}{\partial \bar y^2})
$$


Subsequantly, the mass conservation reads :

$$ \frac{\partial \bar u}{\partial \bar x} = 0  $$

## Boundary conditions 

And the boundary conditions are :

* Adherence at wall $\bar x = 0$

$$ \bar u = \bar v = 0 \quad \text{ at the wall } \bar x = 0$$

* Stress continuity at the free surface, at $\bar x  = 1$

$$
\bar p = 0 \quad \text{ (normal stress continuity)} \\[.5em]
\frac{\partial \bar v}{\partial \bar x} = 0 \quad \text{ (tangential stress continuity)} 
$$

### Parabolic flow profile 

Upon using mass conservation, we see $\bar u$ does not depend on $\bar x$, and therefore $\bar u = \bar u (\bar x = 0) = 0$. As a consequence the $\bar x$-momentum equation reduces to:
$$
\frac{\partial \bar p}{\partial \bar x} = 0
$$
Similarly we find $\bar p$ vanishes throughout the film.

We readily obtain the following solution for $\bar v$:
$$
\bar v = \frac{1}{2} \text{Ar} \bar x^2 + B \bar x + C,
$$
where the constants can be identified with the boundary conditions:
$$
\bar v = \mathrm{Ar}\, \bar x\left(\frac{1}{2} \bar x -1 \right).
$$

### The flowrate

From the preceding analysis we get:
$$
\bar Q = \int_0^1 \bar v (\bar x) \, \mathrm d \bar x = - \frac{1}{3} \text{Ar}
$$
or, under dimensional form:
$$
Q = -\frac{1}{3} \sqrt{g h} h \text{Ar},
$$
which is indeed of the form:
$$
Q = h^{3/2} g^{1/2} F(\text{Ar}),
$$
with $F(x) = -\frac{1}{3} x$.

Using the initial parameters, the flowrate is:
$$
Q = -\frac{1}{3} \frac{\rho g h^3}{\mu}
$$
*/
