/**
# Darcy and Forchheimer terms
This file implements the Darcy and Forchheimer terms for flow in porous media.
They allow to account for the resistance to flow due to the presence of a porous matrix.

In this step, we want to solve the following equation:
$$
\frac{\partial \mathbf{u}}{\partial t} = - \frac{1}{\rho_g}
\left[ \frac{\mu_g \epsilon_g \mathbf{u}_g}{\mathbf{K}} + 
\rho_g\frac{1.75}{\sqrt{150\epsilon_g^3}} \frac{\epsilon_g
|\mathbf{u}_g|\mathbf{u}_g}{\sqrt{\mathbf{K}}} \right]
$$
Note that the permeability tensor **K** (Da in the code) can reach very small values.
Therefore, an implicit treatment of the Darcy and Forchheimer
terms must and is here implemented.

Extern variables defined elsewhere:

+ *porosity*: scalar field representing the porosity of the medium
+ *f*: scalar field representing the volume fraction of the presudo-phase
+ *rhoGv_S*, *muGv_S*: scalar fields for variable density and viscosity
+ *rhoG*, *muG*: double values for constant density and viscosity
*/

#ifndef F_ERR
  #define F_ERR 1e-10
#endif

extern scalar porosity;
extern scalar f;

#ifdef VARPROP
extern scalar rhoGv_S, muGv_S;
#else
extern double rhoG, muG;
#endif

/**
The  permeability tensor **Da** is defiened as a constant vector.
This is to allow for the possibility to have different values in each direction
to account for anisotropic porous media. Units: m^2
*/

coord Da = {1e-10, 1e-10};

/**
## Viscous term event
After the viscous term is computed, this event modifies the velocity field
to account for the Darcy and Forchheimer resistance.
*/

event viscous_term (i++) {
/**
We leave the option to apply a correction before and after the Darcy
as this may be useful in some cases (e.g., high external flows).
*/
#ifndef NO_DARCY_CORRECTION
  correction(dt);
#endif
  foreach() {
    if (f[] > F_ERR) {
      double e = porosity[]/f[];
      double F = 1.75/sqrt (150.*pow (e, 3));
      double Umag = norm(u);

      double muGh, rhoGh;
      #ifdef VARPROP
      muGh = muGv_S[];
      rhoGh = rhoGv_S[];
      #else
      muGh = muG/e; // effective viscosity
      rhoGh = rhoG;
      #endif

      foreach_dimension() {
        double A = muGh*e/(Da.x*rhoGh);     // Darcy term
        double B = F*Umag*e/sqrt(Da.x);     // Forchheimer term
        u.x[] *= exp(-(A + B)*dt*f[]);
      }
    }
  }
#ifndef NO_DARCY_CORRECTION
  correction(-dt);
#endif
}

/**
## Previous implementation (commented out)
Here we used the default acceleration field *a* to add the Darcy and Forchheimer
contributions explicitly. This appoach is more integrated with the existing framework but
less efficient due to the explicit nature of the terms.
*/

// Old approach, works but explicit
// event defaults (i = 0) {
//   if (is_constant(a.x)) {
//     a = new face vector;
//     foreach_face() {
//       a.x[] = 0.;
//       dimensional (a.x[] == Delta/sq(DT));
//     }
//   }
// }

// event acceleration (i++){
//   face vector av = a;
//   foreach_face() {
//     double ff = face_value(f, 0);
//     if (ff > F_ERR) {
//       double ef = face_value(porosity, 0);
//       double F  = 1.75/pow (150*pow (ef, 3), 0.5);
//       double Umag = norm (uf);

//       // Darcy contribution, weighted by the face fraction of the interface
//       av.x[] -= alpha.x[]/(fm.x[] + SEPS)* (muG*ef/Da.x) *uf.x[] *ff; 

//       // Forcheimer contribution
//       av.x[] -= alpha.x[]/(fm.x[] + SEPS)* (F*ef*rhoG/pow(Da.x,0.5)) *Umag  *uf.x[] *ff;
//     }
//   }
// }