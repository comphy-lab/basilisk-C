#include "./../../utils/mathematica.h"

double D = 1e-2;  // diffusion coef
double s=-1.;  // reaction coef \lambda

double exact(double x, double y, double z) {  // analytical solution
  double kappa = 1./(50.*M_PI);
  return exp(-sq((x*y+z*z)/kappa));
}

double src(double x, double y, double z) {  // source term for the Poisson Helmholtz equation
   double u;
   double kappa = 1./(50.*M_PI);
   
   u =   (Power(kappa,4)*s + 4*D*Power(x*y + Power(z,2),2)*(Power(x,2) + Power(y,2) + 4*Power(z,2)) - 2*D*Power(kappa,2)*(Power(x + y,2) + 6*Power(z,2)))/(Power(E,Power(x*y + Power(z,2),2)/Power(kappa,2))*Power(kappa,4));

     return u;
}




  





