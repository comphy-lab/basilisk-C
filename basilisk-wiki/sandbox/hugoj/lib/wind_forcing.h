/*
# Wind forcing

This file implements a wind forcing for a multilayer simulation.

See the test case [here](https://basilisk.fr/sandbox/hugoj/test_windinput/linear_wave_wind_input.c)


Note:

- This file needs to be imported after 'hydro.h' and 'nh.h'
- Ideally this file should be merged with the a_baro overload into a hydro-wind.h
 **/

double pp0 = 0. [-1, -2, 1];

// target energy
double Etarget = 0. [];
// Controler parameters
double K_1 = 1. [0,0,0];
double relax_dt = 0. [0,1];

double K_2 = 0. [0,0,0];
double T_f = 1.;
double sum_I = 0.;

double K_3 = 0. [0,0,0];
double last_dE=0.;

// Low pass filter
double dE_f = 0.;
double sum_f = 0.;
double sum_abs_f = 0.;
double tau_f = 1.;

// p0 clip
double p0_max = 100. [-1,-2,1];


// PI
// response time = 2π·sqrt(Tu/Ku)
// damping = sqrt(Ku·Tu)/2

/**
  ## Forcing using a controler

  IN 2D (x,z):
  input_power = -sum_xy(dp/dx Q)dx
             with p = p0.deta/dx
                  Q = sum_layer(u*h)
      
        so input_power = - p0* sum_xy(sum_layer(d2eta/dx2*u*h*dx))


  Energy deficit = E_ini - E(t)
  */

event update_p0 (i++)
{
  /** Computing the error on wave energy */
  double etam = 0., etavar_current = 0., pe = 0., ke = 0.;
  energy_stats (&etam, &etavar_current, &pe, &ke);
  double dE = Etarget - (pe + ke);
  //double dE = 0., etam=0., etavar_current=0.;
  /**
  Computing input power, P = p0*sum
  */
  double sum = 0.; 
    foreach(reduction(+:sum)) {
    double etaxx = (eta[1] + eta[-1]- 2*eta[])/sq(Delta);
    foreach_layer()
      sum += - h[] * u.x[] * etaxx * dv(); 
  }
  /**
  Low pass filter on error and sum
  */
  dE_f += (dE  - dE_f ) * dt / tau_f;
  sum_f += (sum - sum_f) * dt / tau_f;

  /**
  safeguard: never divide by near zero sum  
  */
  sum_abs_f += (fabs(sum) - sum_abs_f) * dt / tau_f;
  double sum_den = (sum_f >= 0. ? 1. : -1.) * fmax (fabs(sum_f), 0.1*sum_abs_f);
  /**
  PID Control

  p0 = 1/sum * K1( dE + 1/Ti*int_t dE dt + Td*dE/dt ) */
  double rhs = 0.;
  
  /** P: proportional constrain */
  if (i>1){
    rhs += K_1 * dE_f;
  }

  /** I: integral constrain */
  if (fabs((rhs + K_2*sum_I)/sum_den) <= p0_max)
  { 
    // avoid to take into account spikes of energy
    sum_I += dE_f * dt;
  }
  //sum_I += dE_f * dt;
  rhs += K_2 * sum_I;

  /**  D: derivative constrain */  
  double dedt = 0.;  
  dedt = (dE_f-last_dE)/dt;
  rhs += K_3*dedt;
  last_dE = dE_f;
  /** Finally we compute p0 */
  double dp0 = rhs/sum_den;

  p0 = pp0 + dp0;

  /** safeguard 2: p0 limiter */
  if (p0 >  p0_max) p0 =  p0_max;
  if (p0 < -p0_max) p0 = -p0_max;

  // if (i%10==0){
  //   fprintf(stderr, "i=%d etam=%.3e var=%.6e dE=%g sum=%g sum_den=%g p0=%g\n",
  //       i, etam, etavar_current, dE, sum, sum_den, p0);
  // }
  fprintf(stderr, "i=%d etam=%.3e var=%.6e dE=%g sum=%g sum_den=%g p0=%g\n",
        i, etam, etavar_current, dE, sum, sum_den, p0);
}


/** ## TODO 

* document the filters

*/
