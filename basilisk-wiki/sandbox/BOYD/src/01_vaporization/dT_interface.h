
scalar sT_L[];
scalar sT_V[];
#include "intgrad.h"
//trace
void compute_interface_gradient(vector n_vec,
                                scalar interfacial_point) {
  foreach () {
    dT_interface_L[] = 0.0;
    dT_interface_V[] = 0.0;
  }
  
  scalar fL[];
  scalar fG[];

  face vector fsL[];
  face vector fsG[];

  foreach() {
    f[] = clamp (f[], 0., 1.);
    //
    fL[] = f[];
    fG[] = 1. - f[];
  }

  // We compute the value of volume fraction f on the cell-faces using 
  // a geometric approach (necessary for interface gradients and diffusion equations).

  face_fraction (fL, fsL);
  face_fraction (fG, fsG);

  // The calculation of the interface gradients is used also for the calculation of 
  // the source terms for the diffusion equation of the temperature fields.
  foreach() {
    sT_V[] = 0.;
    sT_L[] = 0.;

    if INTERFACIAL {
      coord n = facet_normal (point, fL, fsL), p;
      double alpha = plane_alpha (fL[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

      double gtrgrad = ebmgrad (point, T_V, fL, fG, fsL, fsG, true, T_sat, false);

      double ltrgrad = ebmgrad (point, T_L, fL, fG, fsL, fsG, false, T_sat, false);

      dT_interface_L[] = -ltrgrad;
      dT_interface_V[] = gtrgrad;

      double lheatflux = lambda_L*ltrgrad;
      double gheatflux = lambda_V*gtrgrad;

#ifdef AXI
      sT_L[] = lheatflux/rho_L/cp_L*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sT_V[] = gheatflux/rho_V/cp_V*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
      sT_L[] = lheatflux/rho_L/cp_L*area/Delta*cm[];
      sT_V[] = gheatflux/rho_V/cp_V*area/Delta*cm[];
#endif
    }
  }

  return;
}
