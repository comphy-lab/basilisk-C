/**
# Skalak law for biological membranes

In this file we define the Skalak strain energy function, used to compute
the elastic force on each Lagrangian node of the membrane.

The strain energy function is given by, e.g.,
[Barthès-Biesel](#barthes2016motion):

$$
W^{\text{Skalak}} = \frac{E_s}{4} \left( I_1^2 + 2I_1 - 2I_2 + CI_2^2 \right),
$$

where $E_s$ is the elastic modulus of the membrane, $C > -\frac{1}{2}$ is the
area dilatation modulus preventing strong area changes, and $I_{1,2}$ are two
strain invariants: $I_1 = \lambda_1^2 + \lambda_2^2 - 2$ and $I_2 = \lambda_1^2
\lambda_2^2 - 1 = J^2 - 1$, with $\lambda_{1,2}$ the principal stretches and
$J$ the surface Jacobian representing the ratio of deformed over undeformed
area.

Differentiating $W^{\text{Skalak}}$ with respect to the principal stretches
yields:
$$
\frac{\partial W^{\text{Skalak}}}{\partial \lambda_1} = E_s \left( \lambda_1 (
\lambda_1^2 - 1) + C \lambda_1 \lambda_2^2 (\lambda_1^2 \lambda_2^2 - 1) \right)
$$
and $\frac{\partial W^{\text{Skalak}}}{\partial \lambda_2}$ is obtained by
switching $\lambda_1$ and $\lambda_2$.

## Implementation

In this file we only have to specify the derivatives of the strain function
with respect to the stretches. The computation of the membrane tensions is
implemented in [elasticity-ft.h](elasticity-ft.h).
*/

#ifndef E_S
  #define E_S 1.
#endif
#ifndef AREA_DILATATION_MODULUS
  #define AREA_DILATATION_MODULUS 1.
#endif

#define DWDL1(L1, L2) (E_S*(L1*(sq(L1) - 1.) \
  + AREA_DILATATION_MODULUS*L1*sq(L2)*(sq(L1*L2) - 1.)))
#define DWDL2(L1, L2) (E_S*(L2*(sq(L2) - 1.) \
  + AREA_DILATATION_MODULUS*L2*sq(L1)*(sq(L1*L2) - 1.)))

#include "elasticity-ft.h"

/**
## References
~~~bib
@Article{barthes2016motion,
  author    = {Barthes-Biesel, Dominique},
  journal   = {Annual Review of fluid mechanics},
  title     = {Motion and deformation of elastic capsules and vesicles in flow},
  year      = {2016},
  pages     = {25--52},
  volume    = {48},
  file      = {:barthes2016motion - Motion and Deformation of Elastic Capsules and Vesicles in Flow.pdf:PDF},
  groups    = {Biological flows},
  publisher = {Annual Reviews},
}
~~~
*/
