/**
# Acceleration term in a hyperelastic incompressible material

In a hyperlastic material, the stresses appears when a deformation
occurs in the body. Usually this deformation is commonly described in
a Lagrangian framework although a Eulerian description, that is the
one commonly used for fluids, is possible. In effect, the
[Mooney--Rivlin constitutive law](https://en.wikipedia.org/wiki/Mooney%E2%80%93Rivlin_solid) 
for the hyperelastic stress tensor $\mathbf{S}_h$ writes,
$$ 
\mathbf{S}_h = 2c_1 \mathbf{B} +
2c_2 \, (tr(\mathbf{B}) \, \mathbf{B}-\mathbf{B} \cdot \mathbf{B}) + 
4c_3(tr(\mathbf{B})-3)\mathbf{B}
$$
where $tr()$ denotes the trace of a tensor and the dot denotes tensor
multiplication. $\mathbf{B}$ is the left Cauchy-Green deformation
tensor whose temporal evolution is governed by the upper convective
derivative,
$$ 
UCD(\mathbf{B}) = \partial_t \mathbf{B} 
+ \mathbf{u} \cdot \nabla \mathbf{B} 
- (\nabla \mathbf{u})^T  \cdot \mathbf{B} 
- \mathbf{B} \cdot (\nabla \mathbf{u}) = 0  
$$
$c_1$, $c_2$ and $c_3$ are characterizing coefficients of the
Mooney--Rivlin hyperelastic material. Particular cases of an
(incompressible) Mooney--Rivlin material are: 

* The neo-Hookian material: $c_1 = G/2$; $c_2=c_3 = 0$. 
$G$ is the shear modulus.  
* The Saint-Venant--Kirchoff material: $c_1 = \mu_L$, 
$c_2= -\mu_L/2$ and $c_3 = (\lambda_L +2 \mu_L)/8$.
$\mu_L$ and $\lambda_L$ are the Lame constants.


We will reorder the constitutive equation for the stresses in the
following form,
$$ \mathbf{S}_h = \beta_1 \mathbf{B} + \beta_2
TR(\mathbf{B}) \, \mathbf{B} + \beta_3 \mathbf{B} \cdot \mathbf{B} 
$$
being in 2D $\beta_1 = 2c_1 +2c_2-8c_3$, $\beta_2 = 2c_2+4c_3$ and
$\beta_3 = -2c_3$ and $TR(\mathbf{B}) = B_{xx} + B_{yy}$. Note that in
2D $tr(\mathbf{B}) = B_{xx}+B_{yy} + 1$ since $B_{zz} = 1$. In 3D
$\beta_1$ would be $\beta_1 = 2c_1 -12c_3$.*/

#include "upper.h"

tensor B[];
(const) scalar beta1 = zeroc, beta2 = zeroc, beta3 = zeroc;

event defaults (i = 0) {
  foreach ()
    foreach_dimension ()
      B.x.x[] = 1.;
}

event init (i = 0) {
  if (is_constant(a.x))
    a = new face vector;
}

event tracer_advection (i++) {
  upper_convected_derivative (u, B, dt = dt, theta = zeroc);
}

event acceleration (i++)
{
  
  /**
  The stress tensor $\mathbf{S}_h$ for the Mooney--Rivlin hyperelastic
  material is calculated. */  
  
  tensor S[];
  foreach ()
    foreach_dimension() {
      S.x.x[] = (beta1[] + beta2[]*B.y.y[])*B.x.x[] 
        + beta3[]*B.x.y[]*B.y.x[] + (beta2[] + beta3[])*sq(B.x.x[]);
      S.x.y[] = (beta1[] + (B.x.x[] + B.y.y[])*(beta2[]+beta3[]))*B.x.y[];
    }
  boundary ((scalar *) {S.x, S.y});

  face vector av = a;  
  foreach_face() {
    double shear = (S.x.y[0,1]*cm[0,1] + S.x.y[-1,1]*cm[-1,1] -
                    S.x.y[0,-1]*cm[0,-1] - S.x.y[-1,-1]*cm[-1,-1])/4.;
    av.x[] +=  (fm.x[] == 0. ? 0 :(shear +
                                   cm[]*S.x.x[]-cm[-1,0]*S.x.x[-1,0])
                *alpha.x[]/(fm.x[]*Delta));
  }
}
