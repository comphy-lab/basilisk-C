/*
# Common waves

This file gather expressions of common wave field.

You can look for a surface elevation function, u.x, u.y, u.z

NOTE: g_ has to be either #define or declared !

*/



/*
## Monochromatic linear wave
*/

double wave_monolin(double t, double x, double a, double k){
  double omega = sqrt(g_*k);
  return a*sin(k*x-omega*t);
}

double u_x_monolin(double t, double x, double z, double a, double k){
  double omega = sqrt(g_*k);
  return a*omega*exp(k*z)*sin(k*x-omega*t);
}

double u_y_monolin(double t, double x, double z, double a, double k){
  double omega = sqrt(g_*k);
  return a*omega*exp(k*z)*cos(k*x-omega*t);
}
/*
## 3rd order Stokes wave
Same as src/examples/test/stokes.h but with all necessary argument in the function declaration
*/

double wave_stokes (double x, double y, double ak, double k_, double h_)
{
  //NOTE: 'y' is usually depth
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + ak*eta2 + sq(ak)*eta3 - y;
}

double u_x_stokes (double x, double y, double ak, double k_, double h_)
{
  //NOTE: 'y' is usually depth
  double alpa = 1./tanh(k_*h_);
  double a_ = ak/k_;
  double sgma = sqrt(g_*k_*tanh(k_*h_)*
		     (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					(sq(alpa) - 1.) + sq(alpa))));
  double A_ = a_*g_/sgma;
  return A_*cosh(k_*(y + h_))/cosh(k_*h_)*k_*cos(k_*x) +
    ak*3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
    cosh(2.0*k_*(y + h_))*2.*k_*cos(2.0*k_*x)/cosh(2.0*k_*h_) +
    ak*ak*1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
    (9.*sq(alpa) - 13.)*
    cosh(3.*k_*(y + h_))/cosh(3.*k_*h_)*a_*a_*k_*k_*A_*3.*k_*cos(3.*k_*x);
}

double u_y_stokes (double x, double y, double ak, double k_, double h_)
{
  //NOTE: 'y' is usually depth
  double alpa = 1./tanh(k_*h_);
  double a_ = ak/k_;
  double sgma = sqrt(g_*k_*tanh(k_*h_)*
		     (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					(sq(alpa) - 1.) + sq(alpa))));
  double A_ = a_*g_/sgma;
  return A_*k_*sinh(k_*(y + h_))/cosh(k_*h_)*sin(k_*x) +
    ak*3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
    2.*k_*sinh(2.0*k_*(y + h_))*sin(2.0*k_*x)/cosh(2.0*k_*h_) +
    ak*ak*1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
    (9.*sq(alpa) - 13.)*
    3.*k_*sinh(3.*k_*(y + h_))/cosh(3.*k_*h_)*a_*a_*k_*k_*A_*sin(3.*k_*x);
}
