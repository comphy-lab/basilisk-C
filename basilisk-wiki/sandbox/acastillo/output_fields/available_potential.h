/**
# Functions for compute the available potential energy

## *reference_height()*: Compute the reference height field using a cumulative distribution function

This function computes the reference height `yref` for a scalar field `c`. 
It uses the cumulative distribution function (CDF) of `c` to determine the 
reference height as proposed by [Tseng & Ferzinger (2001)](#tseng2001).

The arguments and their default values are:

*yref*
: Scalar field to store the computed reference height.

*c*
: Scalar field for which the CDF is computed.

*cmin*
: Lower boundary of *c*. Default is `0`.

*cmax*
: Upper boundary of *c*. Default is `1`. 

*store*
: Boolean flag to store the output in a file. Default is `true`. 

*H0*
: Reference height. Default is `L0`. 

The output `ymix` is the reference height at the midpoint of the CDF.

*/

#include "acastillo/output_fields/histograms1D.h"

double reference_height(scalar yref, scalar c, double cmin=0, double cmax=1, bool store=true, double H0=L0, bool reverse=true){

  /* Get the Cumulative distribution function of c. p_pdf isn't used here,
     but probability_distribution_1D() always computes it, so it needs a
     backing array -- heap-allocated since it's discarded right after. */
  double p_range[NBIN], p_sum[NBIN];
  double * p_pdf = malloc(NBIN*sizeof(double));
  probability_distribution_1D(c, cmin, cmax, NBIN, p_range, p_pdf, p_sum,
                               "reference_state_gsl.asc", store);
  free(p_pdf);

  /* Use the CDF to compute the reference height Yref(x,y,z,t) */
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, NBIN);
  gsl_spline_init(spline, p_range, p_sum, NBIN);
  
  double ymix;
  double offset;
#if dimension == 2
  offset = Y0;
#else
  offset = Z0;
#endif

  if (reverse) {
    foreach()
      yref[] = (1.0 - gsl_spline_eval (spline, c[], acc))*(H0/2) + offset;
    ymix = (1.0 - gsl_spline_eval (spline, (cmax+cmin)/2, acc))*(H0/2) + offset;
  } else {
    foreach()
      yref[] = gsl_spline_eval (spline, c[], acc)*(H0/2) + offset;
    ymix = gsl_spline_eval (spline, (cmax+cmin)/2, acc)*(H0/2) + offset;
  }

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return ymix;
}

/**
# References

~~~bib
@article{tseng2001,
  title={Mixing and available potential energy in stratified flows},
  author={Tseng, Yu-heng and Ferziger, Joel H},
  journal={Physics of Fluids},
  volume={13},
  number={5},
  pages={1281--1293},
  year={2001},
  publisher={AIP Publishing},
  doi={10.1063/1.1358307}
}
~~~
*/
