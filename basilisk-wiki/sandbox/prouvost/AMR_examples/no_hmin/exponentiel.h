double D = 1e-2;  // diffusion coef
double s=-1.;  // reaction coef \lambda

double exact(double x, double y, double z) {  // analytical solution
   double c = 3.*M_PI/50.;
   double d = 50.*M_PI;
   return exp(-sq((d*(x*y-c))));
}


double src(double x, double y, double z) {  // source term for the Poisson Helmholtz equation
   double u;
   double c = 3.*M_PI/50.;
   double d = 50.*M_PI;
   
   double f = exp(-sq((d*(x*y-c))));
   u = -2.*D*(sq(d*x) + sq(d*y))*(1. - 2.*d*d*sq(x*y-c)) + s;
   return u*f;
}
