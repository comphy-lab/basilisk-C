double D = 1e-2;  // diffusion coef
double s=-1.;  // reaction coef \lambda

double kappa=1e-2;  // sharpness of the boundary layer

double exact(double x, double y, double z) {  // analytical solution
   return x*sq(y) - sq(y)*exp(2*(x-1)/kappa) - x*exp(3*(y-1)/kappa) + exp( (2*(x-1) + 3*(y-1)) /kappa);
}



double src(double x, double y, double z) {  // source term for the Poisson Helmholtz equation
   
  /* double kappa = D; */
  return  s*x*y*y + 2*D*x + exp(2*(x-1)/kappa)*(-s*sq(y) - 4*D/sq(kappa)*sq(y) - 2*D) \
   + exp(3*(y-1)/kappa)*(-s*x - 9*D/sq(kappa)*x) \
   + exp(2*(x-1)/kappa + 3*(y-1)/kappa)*(s + 13*D/sq(kappa));

	
}



double exact_x_neumann(double x, double y, double z) {  // analytical solution
   return sq(y) - sq(y)*2./kappa*exp(2*(x-1)/kappa) - exp(3*(y-1)/kappa) + 2./kappa*exp( (2*(x-1) + 3*(y-1)) /kappa);
}

double exact_y_neumann(double x, double y, double z) {  // analytical solution
   return 2.*x*y - 2.*y*exp(2*(x-1)/kappa) - 3./kappa*x*exp(3*(y-1)/kappa) + 3./kappa*exp( (2*(x-1) + 3*(y-1)) /kappa);
}






