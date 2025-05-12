#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static const double T_wall = T_sat;
static const double alpha_g = lambda_V / (rho_V * cp_V);
static const double alpha_l = lambda_L / (rho_L * cp_L);
#define const_term (cp_V * (T_wall - T_sat) / hlg)
double beta_guess = 0.7677;

double interface_location_func(double t, double beta) {
  double location = 2.0 * beta * sqrt(alpha_g * t);
  return location;
}

double time_given_interface_func(double x, double beta) {
  double t = sq(x / (2.0 * beta)) / alpha_g;
  return t;
}

double T_vap_func_no_lim(double time, double x, double beta) {
  double T_v = T_wall + (T_sat-T_wall)/(erf(beta))*erf(x/(2.0*sqrt(alpha_g*time)));
  return T_v;
}
double T_vap_func(double time, double x, double beta) {
  double T_v = T_vap_func_no_lim(time, x, beta);
  double interface = interface_location_func(time, beta);
  if (x > 1.01*interface || T_v<T_sat) {
    T_v = T_sat;
  }
  return T_v;
}


double T_liq_func_no_lim(double time, double x, double beta) {
  double erfc1 = erfc(beta * (rho_V * sqrt(alpha_g)) / (rho_L * sqrt(alpha_l)));
  double erfc2 = erfc(x / (2.0 * sqrt(alpha_l * time)) +
           beta * (rho_V - rho_L) / rho_L * sqrt(alpha_g / alpha_l));
  double T_l = T_inf - (T_inf - T_sat) / (erfc1)*erfc2;
  return T_l;
}
//// SHOA et al 2018
double T_liq_func(double time, double x, double beta) {
  double T_l = T_liq_func_no_lim(time, x, beta);
  double interface = interface_location_func(time, beta);
  if (x < 0.99*interface || T_l<T_sat) {
    T_l = T_sat;
  }
  return T_l;
}

double vel_interface_func(double time, double beta) {
  double dt_ = 1.0e-6*time;
  double x1 = interface_location_func(time, beta);
  double x2 = interface_location_func(time + dt_, beta);
  double vel = (x2 - x1) / (dt_);
  return vel;
}

double vel_liq_func(double time, double x, double beta) {
  double interface = interface_location_func(time, beta);
  double vel = vel_interface_func(time, beta);
  if (x < 0.99*interface) {
    vel = 0.0;
  }
  return vel;
}

double dT_interface_func(double time, double beta) {
  double x1 = interface_location_func(time, beta);
  double x2 = x1 + 1.0e-6*x1;
  double T_1 = T_liq_func_no_lim(time, x1, beta);
  double T_2 = T_liq_func_no_lim(time, x2, beta);
  double dT = (T_2 - T_1) / (x2 - x1);
  return dT;
}

double T_liq_linear_test(double time, double x, double beta) {
  double dTdx = 1.0*dT_interface_func(time, beta);
  double interface = interface_location_func(time, beta);

  double x1 = interface;
  double x2 = x1 + (T_inf-T_sat)/dTdx;

  double T_l = dTdx*(x-x1) + T_sat;

  if (x < interface || T_l<T_sat) {
    T_l = T_sat;
  }else if (x > x2 || T_l>T_inf) {
    T_l = T_inf;
  }
  return T_l;
}

//// SHOA et al 2018
double suck_term(double x) {
  double erfc1 = erfc(x * (rho_V * sqrt(alpha_g)) / (rho_L * sqrt(alpha_l)));
  double top = (T_inf - T_sat) * cp_V * lambda_L * sqrt(alpha_g) *
      exp(-x * x * (rho_V * rho_V * alpha_g) / (rho_L * rho_L * alpha_l));
  double bot = hlg * lambda_V * sqrt(M_PI * alpha_l) * erfc1;
  double y = x - top / bot;
  return y;
}

// SHOA et al 2018;
double func_beta(double beta) {
  double y = exp(beta * beta) * erf(beta) * suck_term(beta) - const_term / sqrt(M_PI);
  return y;
}

#include "../common_functions/bisection.h"
#include "../common_functions/secant.h"

double find_beta() {
  double beta = bisection(func_beta, 0.01 * beta_guess, 10.0 * beta_guess, 1.0e-8);
  fprintf(stdout, "beta found to be %g\n", beta);
  return beta;
}

double find_start_time(double interface_init, double beta) {
  double time_init = time_given_interface_func(interface_init, beta);
  fprintf(stdout, "time_init found to be %g\n", time_init);
  return time_init;
}

