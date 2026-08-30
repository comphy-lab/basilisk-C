/**
# 2D histograms and joint probability density functions

## *probability_distribution_2D()*: Joint histogram and PDF of two scalars

This function computes the volume-weighted joint histogram of two scalar
fields *u* and *v* using [GSL 2D histograms](https://www.gnu.org/software/gsl/doc/html/histogram.html#d-histograms).
The joint probability density function (PDF) is derived from it. Cells are
weighted by the cell volume. The function uses parallel reduction to
accumulate values and weights, then normalizes the accumulated values.

Unlike the 1D case, a joint histogram has no single natural cumulative
distribution function (the GSL `gsl_histogram2d_pdf.sum` array is a
row-major flattened cumulative sum over bins, not a marginal CDF), so only
the joint PDF is computed and stored here.

The arguments and their default values are:

*u*
: First scalar field.

*v*
: Second scalar field.

*umin*, *umax*
: Boundaries of *u*. Default is `0`, `1`.

*vmin*, *vmax*
: Boundaries of *v*. Default is `0`, `1`.

*nx*, *ny*
: Number of bins along *u* and *v*. Default is `110`.

*p_xrange*
: Array of size *nx* to store the discretized bin ranges along *u*.

*p_yrange*
: Array of size *ny* to store the discretized bin ranges along *v*.

*p_pdf*
: Array of size *nx* x *ny* (row-major, index `i*ny+j`) to store the joint
  probability density function values.

*filename*
: Name of the file the histogram is appended to. Default is
  `"reference_state_gsl_2d.asc"`.

*store*
: Boolean flag to store the output in a file. Default is `true`.

When *store* is true, each call appends a block to *filename*: a two-line
comment header naming the scalars (`u.name`, `v.name`) and the columns,
followed by one row per bin (in row-major, i.e. blank-separated by *i*,
suitable for gnuplot's `pm3d`/`splot` with `every :::i::i` or simply piped
through `pm3d map`),

~~~
# u = <u.name>  v = <v.name>
# t xrange[<u.name>] yrange[<v.name>] pdf[<u.name>,<v.name>]
t  p_xrange[0]  p_yrange[0]  p_pdf[0,0]
t  p_xrange[0]  p_yrange[1]  p_pdf[0,1]
...
t  p_xrange[0]  p_yrange[ny-1]  p_pdf[0,ny-1]

t  p_xrange[1]  p_yrange[0]  p_pdf[1,0]
...
~~~

so that repeated calls (e.g. one per timestep, or for different scalar
pairs) can be told apart when appended to the same file.
*/

#include <gsl/gsl_histogram2d.h>
#pragma autolink -lgsl -lgslcblas
#define NBIN2D 110


void probability_distribution_2D(scalar u, scalar v, double umin=0, double umax=1,
                                  double vmin=0, double vmax=1,
                                  int nx=NBIN2D, int ny=NBIN2D,
                                  double p_xrange[nx], double p_yrange[ny],
                                  double p_pdf[nx*ny],
                                  const char * filename = "reference_state_gsl_2d.asc",
                                  bool store=true){

  double urange = (umax - umin);
  double vrange = (vmax - vmin);
  gsl_histogram2d * h = gsl_histogram2d_alloc (nx, ny);
  gsl_histogram2d_set_ranges_uniform (h, umin-0.05*urange, umax+0.05*urange,
                                          vmin-0.05*vrange, vmax+0.05*vrange);
  double du = h->xrange[1] - h->xrange[0];
  double dv_ = h->yrange[1] - h->yrange[0];

  double p_pdf_loc[nx*ny];
  for (int ii = 0; ii < nx; ii++)
    p_xrange[ii] = h->xrange[ii];
  for (int jj = 0; jj < ny; jj++)
    p_yrange[jj] = h->yrange[jj];
  for (int kk = 0; kk < nx*ny; kk++)
    p_pdf_loc[kk] = 0;

  /* Obtain the total volume per pid */
  double vol_cells=0.0, vol_cells_pid=0.0;
  foreach(serial, noauto)
  #if EMBED
    vol_cells_pid += dv()*cs[];
  #else
    vol_cells_pid += dv();
  #endif

  /* Populate a 2D histogram using fields weighted by volume.
     gsl_histogram2d_accumulate() mutates h in place with no reduction clause
     Basilisk could parallelize around, so this loop must stay serial -- a
     plain foreach() here would race on h->bin[] under OpenMP. */
  if (vol_cells_pid > 0){
    foreach(serial, noauto)
    #if EMBED
      gsl_histogram2d_accumulate(h, u[], v[], dv()*cs[]);
    #else
      gsl_histogram2d_accumulate(h, u[], v[], dv());
    #endif

    for (int ii = 0; ii < nx; ii++)
      for (int jj = 0; jj < ny; jj++)
        p_pdf_loc[ii*ny+jj] = gsl_histogram2d_get (h, ii, jj);
  }
  gsl_histogram2d_free (h);

  /* Sum the values over the ensemble of sub-domains */
@if _MPI
  MPI_Allreduce (&p_pdf_loc[0],  &p_pdf[0],  nx*ny, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&vol_cells_pid, &vol_cells, 1    , MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
@else
  vol_cells = vol_cells_pid;
  for (int kk = 0; kk < nx*ny; kk++)
    p_pdf[kk] = p_pdf_loc[kk];
@endif
  /* We divide by the total volume to ensure the pdf integrates to 1 over
     [umin-0.05*urange, umax+0.05*urange] x [vmin-0.05*vrange, vmax+0.05*vrange] */
  for (int kk = 0; kk < nx*ny; kk++)
    p_pdf[kk] /= vol_cells*du*dv_;

  if ((pid() == 0) & (store)) {
    FILE * fp = fopen(filename, "a") ;
    fprintf (fp, "# u = %s  v = %s\n# t xrange[%s] yrange[%s] pdf[%s,%s]\n",
             u.name, v.name, u.name, v.name, u.name, v.name);
    for (int i = 0; i < nx; i++){
      for (int j = 0; j < ny; j++)
        fprintf (fp, "%.12g %.12g %.12g %.12g \n", t, p_xrange[i], p_yrange[j], p_pdf[i*ny+j]);
      fputs ("\n", fp);
    }
    fputs ("\n", fp);
    fclose (fp);
  }
}

/**
## *probability_distribution_2D_weighted()*: Joint histogram and PDF of two scalars, restricted to a region

Same as `probability_distribution_2D()`, but restricted to a region-of-interest
*f* (e.g. a VOF fraction): cells are weighted by *f* times the cell volume
instead of the cell volume alone.

The arguments and their default values are the same as
`probability_distribution_2D()`, with one addition:

*f*
: Scalar field defining the region-of-interest. Required, no default.
*/

void probability_distribution_2D_weighted(scalar f, scalar u, scalar v,
                                           double umin=0, double umax=1,
                                           double vmin=0, double vmax=1,
                                           int nx=NBIN2D, int ny=NBIN2D,
                                           double p_xrange[nx], double p_yrange[ny],
                                           double p_pdf[nx*ny],
                                           const char * filename = "reference_state_gsl_2d.asc",
                                           bool store=true){

  double urange = (umax - umin);
  double vrange = (vmax - vmin);
  gsl_histogram2d * h = gsl_histogram2d_alloc (nx, ny);
  gsl_histogram2d_set_ranges_uniform (h, umin-0.05*urange, umax+0.05*urange,
                                          vmin-0.05*vrange, vmax+0.05*vrange);
  double du = h->xrange[1] - h->xrange[0];
  double dv_ = h->yrange[1] - h->yrange[0];

  double p_pdf_loc[nx*ny];
  for (int ii = 0; ii < nx; ii++)
    p_xrange[ii] = h->xrange[ii];
  for (int jj = 0; jj < ny; jj++)
    p_yrange[jj] = h->yrange[jj];
  for (int kk = 0; kk < nx*ny; kk++)
    p_pdf_loc[kk] = 0;

  /* Obtain the total weighted volume per pid */
  double vol_cells=0.0, vol_cells_pid=0.0;
  foreach(serial, noauto)
  #if EMBED
    vol_cells_pid += f[]*dv()*cs[];
  #else
    vol_cells_pid += f[]*dv();
  #endif

  /* Populate a 2D histogram using fields weighted by f times the cell volume.
     gsl_histogram2d_accumulate() mutates h in place with no reduction clause
     Basilisk could parallelize around, so this loop must stay serial -- a
     plain foreach() here would race on h->bin[] under OpenMP. */
  if (vol_cells_pid > 0){
    foreach(serial, noauto)
    #if EMBED
      gsl_histogram2d_accumulate(h, u[], v[], f[]*dv()*cs[]);
    #else
      gsl_histogram2d_accumulate(h, u[], v[], f[]*dv());
    #endif

    for (int ii = 0; ii < nx; ii++)
      for (int jj = 0; jj < ny; jj++)
        p_pdf_loc[ii*ny+jj] = gsl_histogram2d_get (h, ii, jj);
  }
  gsl_histogram2d_free (h);

  /* Sum the values over the ensemble of sub-domains */
@if _MPI
  MPI_Allreduce (&p_pdf_loc[0],  &p_pdf[0],  nx*ny, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&vol_cells_pid, &vol_cells, 1    , MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
@else
  vol_cells = vol_cells_pid;
  for (int kk = 0; kk < nx*ny; kk++)
    p_pdf[kk] = p_pdf_loc[kk];
@endif
  /* We divide by the total volume to ensure the pdf integrates to 1 over
     [umin-0.05*urange, umax+0.05*urange] x [vmin-0.05*vrange, vmax+0.05*vrange] */
  for (int kk = 0; kk < nx*ny; kk++)
    p_pdf[kk] /= vol_cells*du*dv_;

  if ((pid() == 0) & (store)) {
    FILE * fp = fopen(filename, "a") ;
    fprintf (fp, "# u = %s  v = %s\n# t xrange[%s] yrange[%s] pdf[%s,%s]\n",
             u.name, v.name, u.name, v.name, u.name, v.name);
    for (int i = 0; i < nx; i++){
      for (int j = 0; j < ny; j++)
        fprintf (fp, "%.12g %.12g %.12g %.12g \n", t, p_xrange[i], p_yrange[j], p_pdf[i*ny+j]);
      fputs ("\n", fp);
    }
    fputs ("\n", fp);
    fclose (fp);
  }
}
