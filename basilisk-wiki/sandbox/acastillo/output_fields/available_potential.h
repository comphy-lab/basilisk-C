/**
# Functions for compute the available potential energy

## *get_distribution_function()*: Cumulative distribution function of an scalar

This function calculates the cumulative distribution function of a scalar field
*c* within a specified region defined by *f* using histograms. Fields are
weighted by the cell volume. The function uses parallel reduction to accumulate
values and weights, then normalizes the accumulated values. 

The arguments and their default values are:

*f*
: Scalar field defining the region-of-interest.

*c*
: Scalar field for which the CDF is computed.

*cmin*
: Lower boundary of *c*. Default is `0`.

*cmax*
: Upper boundary of *c*. Default is `1`. 

*nbin*
: Number of points (bins) between `cmin` and `cmax`. Default is `220`.

*p_range*
: Array to store the discretized bin ranges.

*p_sum*
: Array to store the cumulative function values.

*store*
: Boolean flag to store the output in a file. Default is `true`. 
*/


#include <gsl/gsl_spline.h>
#include <gsl/gsl_histogram.h>
#pragma autolink -lgsl -lgslcblas
#define NBIN 220


void get_distribution_function(scalar f, scalar c, double cmin=0, double cmax=1, int nbin=NBIN, double p_range[nbin], double p_sum[nbin], bool store=true){

  double crange = (cmax - cmin);
  /* Set the number of bins */
  gsl_histogram * h = gsl_histogram_alloc (nbin);
  gsl_histogram_set_ranges_uniform (h, cmin-0.05*crange, cmax+0.05*crange);

  double p_sum_loc[nbin];
  for (int ii = 0; ii < nbin; ii++){
    p_range[ii] = h->range[ii];
    p_sum_loc[ii] = 0;
  }

  /* Obtain the total volume per pid */
  double vol_cells=0.0, vol_cells_pid=0.0;
  foreach(serial, noauto)
     vol_cells_pid += f[]*dv();

  /* Populate a 1D histogram using fields weighted by volume */
  if (vol_cells_pid > 0){
    foreach()
      gsl_histogram_accumulate(h, c[], f[]*dv());

    /* Obtain the cumulative distribution from the histogram */
    gsl_histogram_pdf * p = gsl_histogram_pdf_alloc (nbin);
    gsl_histogram_pdf_init (p, h);

    /* We multiply the function by the volume per pid */
    for (int ii = 0; ii < nbin; ii++){
      p_sum_loc[ii] = p->sum[ii]*vol_cells_pid;
    }
    gsl_histogram_pdf_free (p);
  }
  gsl_histogram_free (h);

  /* Sum the values over the ensemble of sub-domains */
@if _MPI
  MPI_Allreduce (&p_sum_loc[0],  &p_sum[0],  nbin, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&vol_cells_pid, &vol_cells, 1   , MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
@else
  vol_cells = vol_cells_pid;
  for (int ii = 0; ii < nbin; ii++)
    p_sum[ii] = p_sum_loc[ii];
@endif
  /* We divide by the total volume to ensure the cdf has values \in [0,1] */
  for (int ii = 0; ii < nbin; ii++)
    p_sum[ii] /= vol_cells;

  for (int ii = 0; ii < nbin; ii++){
    if ( fabs(p_range[ii] - cmin) < 1e-10 ){
      p_sum[ii] = 0.5*(p_sum[ii]+p_sum[ii+1]);
    }
    if ( fabs(p_range[ii] - cmax) < 1e-10 ){
      p_sum[ii] = 0.5*(p_sum[ii]+p_sum[ii-1]);
    }
  }

  if ((pid() == 0) & (store)) {
    char name[80];
    sprintf(name, "reference_state_gsl_%d.asc", pid());
    FILE * fp = fopen(name, "a") ;
    for (int i = 0; i < nbin; i++)
      fprintf (fp, "%g %g %g \n", t, p_range[i], p_sum[i]);
    fputs ("\n", fp);
    fclose (fp);
  }
}

/**
## *reference_height()*: Compute the reference height field using a cumulative distribution function

This function computes the reference height `yref` for a scalar field `c` within
a specified region defined by the scalar field `f`. It uses the cumulative
distribution function (CDF) of `c` to determine the reference height as prposed 
by [Tseng & Ferzinger (2001)](#tseng2001).

The arguments and their default values are:

*yref*
: Scalar field to store the computed reference height.

*f*
: Scalar field defining the region-of-interest.

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

double reference_height(scalar yref, scalar f, scalar c, double cmin=0, double cmax=1, bool store=true, double H0=L0){

  /* Get the Cumulative distribution function of c*/
  double p_range[NBIN], p_sum[NBIN];
  get_distribution_function(f, c, cmin, cmax, NBIN, p_range, p_sum, store);

  /* Use the CDF to compute the reference height Yref(x,y,z,t) */
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, NBIN);
  gsl_spline_init(spline, p_range, p_sum, NBIN);
  foreach()
    yref[] = (1.0 - gsl_spline_eval (spline, c[], acc))*(H0/2) + Y0;

  double ymix = (1.0 - gsl_spline_eval (spline, (cmax+cmin)/2, acc))*(H0/2) + Y0;
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