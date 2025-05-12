#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef T_inf
#define T_inf T_wall
#endif
static const double dT_const = T_inf - T_sat;
static const double alpha_l = lambda_L / (rho_L * cp_L);
double beta_guess = 0.78;

double bubble_radius(double t, double beta) {
  double radius = 2.0 * beta * sqrt(alpha_l * t);
  return radius;
}
double time_given_radius(double radius_, double beta) {
  double time = pow(0.5 * radius_ / beta, 2.0) / alpha_l;
  return time;
}
double integrand(double x, double beta) {
  double pow_ = pow(max((1. - x), 1.0e-12), -2.);
  double const_ = (pow_ - 2. * (1. - rho_V / rho_L) * x - 1.);
  double val = exp(-beta * beta * const_);

  return val;
}

#include "../common_functions/trapezoid2.h"
double T_liq_func_no_lim(double time, double x, double beta) {
  double Radius = bubble_radius(time, beta);
  double T_integral = qtrap(integrand, 1.0 - Radius / x, 1.0, beta);
  double const1 = (rho_V * (hlg + (cp_L - cp_V) * dT_const)) / (rho_L * cp_L);
  double T_L = T_inf - 2. * beta * beta * const1 * T_integral;
  return T_L;
}
//// Malan 2021
double T_liq_func(double time, double x, double beta) {
  double Radius = bubble_radius(time, beta);

  if (x < 0.99*Radius) {
    return T_sat;
  }
  double T_L = T_liq_func_no_lim(time, x, beta);
  if (T_L < T_sat) {
    return T_sat;
  }
  return T_L;
}
double dT_interface_func(double time, double beta) {
  double x1 = bubble_radius(time, beta);
  double x2 = x1 + 1.0e-6*x1;
  double T_1 = T_liq_func_no_lim(time, x1, beta);
  double T_2 = T_liq_func_no_lim(time, x2, beta);
  double dT = (T_2 - T_1) / (x2 - x1);
  return dT;
}

double bubble_func(double beta) {
  double top = rho_L * cp_L * dT_const;
  double bot = rho_V * (hlg + (cp_L - cp_V) * dT_const);
  double iintegral = qtrap(integrand, 0.0, 1.0, beta);
  double y = top / bot - 2.0 * beta * beta * iintegral;
  return y;
}

double func_beta(double x) {
  double y = bubble_func(x);
  return y;
}

#include "../common_functions/secant.h"
double find_beta() {
  double beta = secant(func_beta, 0.5 * beta_guess, 1.5 * beta_guess, 1.0e-6);
  fprintf(stdout, "beta found to be %g\n", beta);
  return beta;
}

