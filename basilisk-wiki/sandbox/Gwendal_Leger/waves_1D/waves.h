#include "dalloc.h"
#include "waves_input.h"


typedef struct {
  int N; // Number of waves
  double * A; // Amplitudes
  double * lambda; // Wavelengths
  double * theta; // Angle (radians) of spatial diffusion
} wave_spectrum;




wave_spectrum init_waves () {
  wave_spectrum waves;
  int ii;

  waves.N = waves_N;
  waves.A = dalloc (waves.N);
  waves.lambda = dalloc (waves.N);
  waves.theta = dalloc (waves.N);

  double cumsum = 0.;
  for (ii = 0; ii < waves.N; ii++) {
    waves.A[ii] = waves_A[ii];
    waves.lambda[ii] = waves_lambda[ii];
    waves.theta[ii] = waves_theta[ii];
    cumsum += waves.A[ii];
  }

  for (ii = 0; ii < waves.N; ii++) {
    if (waves_normalize)
      waves.A[ii] /= cumsum;
    waves.A[ii] *= waves_amp;
  }

  FILE * fp = fopen ("wavedata", "w");
  fprintf (fp,
	   "\n+-------------WAVES-------------+\n"
	   "Number of waves = %i\n"
	   "Amplification factor = %g\n",
	   waves.N, waves_amp);
  if (waves_normalize)
    fprintf (fp, "Normalized wave amplitudes.\n");
  else
    fprintf (fp, "Un-normalized wave amplitudes.\n");
  fprintf (fp, "A  lambda  theta\n");
  for (ii = 0; ii < waves.N; ii++)
    fprintf (fp, "%g  %g  %g\n", waves.A[ii], waves.lambda[ii], waves.theta[ii]);
  fprintf (fp, "+-------------------------------+\n");
  fclose (fp);
  
  return waves;
}


double VeloX (double x, double y, double z, double t, wave_spectrum waves, double H) {
  double v = 0.;
  int ii;
  for (ii = 0; ii < waves.N; ii++) {
    double k = 2.*pi/waves.lambda[ii];
    double om = sqrt (G*k*tanh(k*H));
    v += waves.A[ii]*om * cosh(k*z)/sinh(k*H) * cos(waves.theta[ii]) * sin(om*t-k*(cos(waves.theta[ii])*x + sin(waves.theta[ii])*y));
    //v += waves.A[ii] * cos(waves.theta[ii]) * sin(om*t - k*(cos(waves.theta[ii])*x + sin(waves.theta[ii])*y));
  }
  //fprintf (stderr, "x=%g  y=%g  z=%g  t=%g  H=%g  vx=%g\n", x, y, z, t, H, v);
  return v;
}


double VeloY (double x, double y, double z, double t, wave_spectrum waves, double H) {
  double v = 0.;
  int ii;
  for (ii = 0; ii < waves.N; ii++) {
    double k = 2.*pi/waves.lambda[ii];
    double om = sqrt (G*k*tanh(k*H));
    v += waves.A[ii] * sin(waves.theta[ii]) * sin(om*t - k*(cos(waves.theta[ii])*x + sin(waves.theta[ii])*y));
  }
  return v;
}


double VeloZ (double x, double y, double z, double t, wave_spectrum waves, double H) {
  double v = 0.;
  int ii;
  for (ii = 0; ii < waves.N; ii++) {
    double k = 2.*pi/waves.lambda[ii];
    double om = sqrt (G*k*tanh(k*H));
    v += waves.A[ii]*om * sinh(k*z)/sinh(k*H) * cos (om*t - k*(cos(waves.theta[ii])*x + sin(waves.theta[ii])*y));
    //v += waves.A[ii] * cos(om*t - k*(cos(waves.theta[ii])*x + sin(waves.theta[ii])*y));
  }
  return v;
}


double SurfElev (double x, double y, double t, wave_spectrum waves, double H) {
  double h = 0.;
  int ii;
  for (ii = 0; ii < waves.N; ii++) {
    double k = 2.*pi/waves.lambda[ii];
    double om = sqrt (G*k*tanh(k*H));
    h += waves.A[ii] * sin(om*t - k*(cos(waves.theta[ii])*x + sin(waves.theta[ii])*y));
    //h += waves.A[ii] * cos(om*t - k*(cos(waves.theta[ii])*x + sin(waves.theta[ii])*y));
  }
  h += H;
  return h;
}

