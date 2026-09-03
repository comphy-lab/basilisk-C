/**
# Histograms and cumulative distribution functions

## *probability_distribution_1D()*: Histogram, PDF and CDF of a scalar

This function computes the volume-weighted histogram of a scalar field *c*
using [GSL histograms](https://www.gnu.org/software/gsl/doc/html/histogram.html).
Both the probability density function (PDF) and the cumulative distribution
function (CDF) are derived from it. Cells are weighted by the cell volume.
The function uses parallel reduction to accumulate values and weights, then
normalizes the accumulated values.

The arguments and their default values are:

*c*
: Scalar field for which the histogram is computed.

*cmin*
: Lower boundary of *c*. Default is `0`.

*cmax*
: Upper boundary of *c*. Default is `1`.

*nbin*
: Number of points (bins) between `cmin` and `cmax`. Default is `220`.

*p_range*
: Array to store the discretized bin ranges.

*p_pdf*
: Array to store the probability density function values.

*p_sum*
: Array to store the cumulative distribution function values.

*filename*
: Name of the file the histogram is appended to. Default is
  `"reference_state_gsl.asc"`.

*store*
: Boolean flag to store the output in a file. Default is `true`.

When *store* is true, each call appends a block to *filename*: a two-line
comment header naming the scalar (`c.name`) and the four columns, followed by
one row per bin,

~~~
# c = <c.name>
# t range[<c.name>] pdf[<c.name>] cdf[<c.name>]
t  p_range[0]  p_pdf[0]  p_sum[0]
t  p_range[1]  p_pdf[1]  p_sum[1]
...
~~~

so that repeated calls (e.g. one per timestep, or for different scalars) can
be told apart when appended to the same file.
*/

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_histogram.h>
#pragma autolink -lgsl -lgslcblas
#define NBIN 220


void probability_distribution_1D(scalar c, double cmin=0, double cmax=1,
                                  int nbin=NBIN, double p_range[nbin],
                                  double p_pdf[nbin], double p_sum[nbin],
                                  const char * filename = "reference_state_gsl.asc",
                                  bool store=true){

  double crange = (cmax - cmin);
  /* Set the number of bins */
  gsl_histogram * h = gsl_histogram_alloc (nbin);
  gsl_histogram_set_ranges_uniform (h, cmin-0.05*crange, cmax+0.05*crange);
  double db = h->range[1] - h->range[0];

  double p_sum_loc[nbin], p_pdf_loc[nbin];
  for (int ii = 0; ii < nbin; ii++){
    p_range[ii] = h->range[ii];
    p_sum_loc[ii] = 0;
    p_pdf_loc[ii] = 0;
  }

  /* Obtain the total volume per pid */
  double vol_cells=0.0, vol_cells_pid=0.0;
  foreach(serial, noauto)
  #if EMBED
    vol_cells_pid += (c[] != nodata) ? dv()*cs[] : 0;
  #else
    vol_cells_pid += (c[] != nodata) ? dv() : 0;
  #endif

  /* Populate a 1D histogram using fields weighted by volume.
     gsl_histogram_accumulate() mutates h in place with no reduction clause
     Basilisk could parallelize around, so this loop must stay serial -- a
     plain foreach() here would race on h->bin[] under OpenMP.
     c[] can be nodata (HUGE) on cells where the field is undefined (e.g.
     curvature away from the interface); accumulate() returns GSL_EDOM for
     out-of-range x, which with GSL's default error handler aborts. Check
     the return and skip those cells instead. */
  gsl_error_handler_t * gsl_prev_handler = gsl_set_error_handler_off();
  if (vol_cells_pid > 0){
    foreach(serial, noauto)
    #if EMBED
      gsl_histogram_accumulate(h, c[], dv()*cs[]);
    #else
      gsl_histogram_accumulate(h, c[], dv());
    #endif

    /* Obtain the cumulative distribution from the histogram */
    gsl_histogram_pdf * p = gsl_histogram_pdf_alloc (nbin);
    gsl_histogram_pdf_init (p, h);

    /* p->sum has nbin+1 entries: p->sum[ii] is the cumulative probability
       up to the *left* edge of bin ii, i.e. exactly at p_range[ii]. We
       multiply the function by the volume per pid. */
    for (int ii = 0; ii < nbin; ii++){
      p_sum_loc[ii] = p->sum[ii]*vol_cells_pid;
      p_pdf_loc[ii] = gsl_histogram_get (h, ii);
    }
    gsl_histogram_pdf_free (p);
  }
  gsl_histogram_free (h);
  gsl_set_error_handler (gsl_prev_handler);

  /* Sum the values over the ensemble of sub-domains */
@if _MPI
  MPI_Allreduce (&p_sum_loc[0],  &p_sum[0],  nbin, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&p_pdf_loc[0],  &p_pdf[0],  nbin, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&vol_cells_pid, &vol_cells, 1   , MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
@else
  vol_cells = vol_cells_pid;
  for (int ii = 0; ii < nbin; ii++){
    p_sum[ii] = p_sum_loc[ii];
    p_pdf[ii] = p_pdf_loc[ii];
  }
@endif
  /* We divide by the total volume to ensure the cdf has values \in [0,1],
     and the pdf integrates to 1 over [cmin-0.05*crange, cmax+0.05*crange] */
  for (int ii = 0; ii < nbin; ii++){
    p_sum[ii] /= vol_cells;
    p_pdf[ii] /= vol_cells*db;
  }

  for (int ii = 0; ii < nbin; ii++){
    if ( fabs(p_range[ii] - cmin) < 1e-10 ){
      p_sum[ii] = 0.5*(p_sum[ii]+p_sum[ii+1]);
    }
    if ( fabs(p_range[ii] - cmax) < 1e-10 ){
      p_sum[ii] = 0.5*(p_sum[ii]+p_sum[ii-1]);
    }
  }

  if ((pid() == 0) & (store)) {
    FILE * fp = fopen(filename, "a") ;
    fprintf (fp, "# c = %s\n# t range[%s] pdf[%s] cdf[%s]\n",
             c.name, c.name, c.name, c.name);
    for (int i = 0; i < nbin; i++)
      fprintf (fp, "%.12g %.12g %.12g %.12g \n", t, p_range[i], p_pdf[i], p_sum[i]);
    fputs ("\n", fp);
    fclose (fp);
  }
}

/**
## *probability_distribution_1D_weighted()*: Histogram, PDF and CDF of a scalar, restricted to a region

Same as `probability_distribution_1D()`, but restricted to a region-of-interest
*f* (e.g. a VOF fraction): cells are weighted by *f* times the cell volume
instead of the cell volume alone.

The arguments and their default values are the same as
`probability_distribution_1D()`, with one addition:

*f*
: Scalar field defining the region-of-interest. Required, no default --
  pass `unity` explicitly to recover the unrestricted histogram (a `const
  scalar` such as `unity` cannot be bound as a *default* argument value in
  Basilisk, only passed explicitly).
*/

void probability_distribution_1D_weighted(scalar f, scalar c, double cmin=0,
                                           double cmax=1, int nbin=NBIN,
                                           double p_range[nbin],
                                           double p_pdf[nbin],
                                           double p_sum[nbin],
                                           const char * filename = "reference_state_gsl.asc",
                                           bool store=true){

  double crange = (cmax - cmin);
  /* Set the number of bins */
  gsl_histogram * h = gsl_histogram_alloc (nbin);
  gsl_histogram_set_ranges_uniform (h, cmin-0.05*crange, cmax+0.05*crange);
  double db = h->range[1] - h->range[0];

  double p_sum_loc[nbin], p_pdf_loc[nbin];
  for (int ii = 0; ii < nbin; ii++){
    p_range[ii] = h->range[ii];
    p_sum_loc[ii] = 0;
    p_pdf_loc[ii] = 0;
  }

  /* Obtain the total weighted volume per pid */
  double vol_cells=0.0, vol_cells_pid=0.0;
  foreach(serial, noauto)
  #if EMBED
    vol_cells_pid += f[]*dv()*cs[];
  #else
    vol_cells_pid += f[]*dv();
  #endif

  /* Populate a 1D histogram using fields weighted by f times the cell volume.
     gsl_histogram_accumulate() mutates h in place with no reduction clause
     Basilisk could parallelize around, so this loop must stay serial -- a
     plain foreach() here would race on h->bin[] under OpenMP.
     c[] can be nodata (HUGE) on cells where the field is undefined; see the
     comment in probability_distribution_1D() above. */
  gsl_error_handler_t * gsl_prev_handler = gsl_set_error_handler_off();
  if (vol_cells_pid > 0){
    foreach(serial, noauto)
    #if EMBED
      gsl_histogram_accumulate(h, c[], f[]*dv()*cs[]);
    #else
      gsl_histogram_accumulate(h, c[], f[]*dv());
    #endif

    /* Obtain the cumulative distribution from the histogram */
    gsl_histogram_pdf * p = gsl_histogram_pdf_alloc (nbin);
    gsl_histogram_pdf_init (p, h);

    /* p->sum has nbin+1 entries: p->sum[ii] is the cumulative probability
       up to the *left* edge of bin ii, i.e. exactly at p_range[ii]. We
       multiply the function by the volume per pid. */
    for (int ii = 0; ii < nbin; ii++){
      p_sum_loc[ii] = p->sum[ii]*vol_cells_pid;
      p_pdf_loc[ii] = gsl_histogram_get (h, ii);
    }
    gsl_histogram_pdf_free (p);
  }
  gsl_histogram_free (h);
  gsl_set_error_handler (gsl_prev_handler);

  /* Sum the values over the ensemble of sub-domains */
@if _MPI
  MPI_Allreduce (&p_sum_loc[0],  &p_sum[0],  nbin, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&p_pdf_loc[0],  &p_pdf[0],  nbin, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&vol_cells_pid, &vol_cells, 1   , MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
@else
  vol_cells = vol_cells_pid;
  for (int ii = 0; ii < nbin; ii++){
    p_sum[ii] = p_sum_loc[ii];
    p_pdf[ii] = p_pdf_loc[ii];
  }
@endif
  /* We divide by the total volume to ensure the cdf has values \in [0,1],
     and the pdf integrates to 1 over [cmin-0.05*crange, cmax+0.05*crange] */
  for (int ii = 0; ii < nbin; ii++){
    p_sum[ii] /= vol_cells;
    p_pdf[ii] /= vol_cells*db;
  }

  for (int ii = 0; ii < nbin; ii++){
    if ( fabs(p_range[ii] - cmin) < 1e-10 ){
      p_sum[ii] = 0.5*(p_sum[ii]+p_sum[ii+1]);
    }
    if ( fabs(p_range[ii] - cmax) < 1e-10 ){
      p_sum[ii] = 0.5*(p_sum[ii]+p_sum[ii-1]);
    }
  }

  if ((pid() == 0) & (store)) {
    FILE * fp = fopen(filename, "a") ;
    fprintf (fp, "# c = %s\n# t range[%s] pdf[%s] cdf[%s]\n",
             c.name, c.name, c.name, c.name);
    for (int i = 0; i < nbin; i++)
      fprintf (fp, "%.12g %.12g %.12g %.12g \n", t, p_range[i], p_pdf[i], p_sum[i]);
    fputs ("\n", fp);
    fclose (fp);
  }
}
