#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// double T_inf = 378.15;
static const double T_inf = T_wall;
static const double alpha_g = lambda_V / (rho_V * cp_V);
static const double alpha_l = lambda_L / (rho_L * cp_L);
#define const_term (cp_V * (T_inf - T_sat) / (hlg * sqrt(M_PI)))
double beta_guess = 0.0667;

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
//// SHOA et al 2018
double T_vap_func(double time, double x, double beta) {
    double T_v = T_vap_func_no_lim(time, x, beta);
    double interface = interface_location_func(time, beta);
    if (x > 1.01*interface || T_v<T_sat) {
        T_v = T_sat;
    }
    return T_v;
}

double dT_interface_func(double time, double beta) {
    double x2 = interface_location_func(time, beta);
    //printf("%d\n", x1);
    double x1 = x2 - 1.0e-6*x2;
    double T_1 = T_vap_func_no_lim(time, x1, beta);
    double T_2 = T_vap_func_no_lim(time, x2, beta);
    double dT = (T_2 - T_1) / (x2 - x1);
    return dT;
}

//// SHOA et al 2018
double suck_term(double x) {
    double erfc1 = erfc(x * (rho_V * sqrt(alpha_g)) / (rho_L * sqrt(alpha_l)));
    double top =
        (T_inf - T_sat) * cp_V * lambda_L * sqrt(alpha_g) *
        exp(-x * x * (rho_V * rho_V * alpha_g) / (rho_L * rho_L * alpha_l));
    double bot = hlg * lambda_V * sqrt(M_PI * alpha_l) * erfc1;
    double y = x - top / bot;
    return y;
}

// SHOA et al 2018;
double func_beta(double x) {
    double y = x * exp(x*x) * erf(x) - const_term;
    return y;
}

#include "../common_functions/bisection.h"
#include "../common_functions/secant.h"

double find_beta() {
    // double beta = secant(func_beta, 0.5 * beta_guess, 1.5 *
    // beta_guess, 1.0e-6);
    double beta =
        bisection(func_beta, 0.01 * beta_guess, 10.0 * beta_guess, 1.0e-8);
    fprintf(stderr, "beta found to be %g\n", beta);
    // print("found beta ", beta);
    return beta;
}

double find_start_time(double interface_init, double beta) {
    double time_init = time_given_interface_func(interface_init, beta);
    fprintf(stderr, "time_init found to be %g\n", time_init);
    return time_init;
}

