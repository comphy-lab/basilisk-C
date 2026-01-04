#include <gsl/gsl_spline.h>
#include <gsl/gsl_histogram.h>
#pragma autolink -lgsl -lgslcblas

#define NINT(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))
void reference_height(scalar yref, scalar T, double thetamin, double thetamax, int nbin){
	/* Compute the volume weighted PDF of the temperature field */
	gsl_histogram * h = gsl_histogram_alloc (nbin);
	gsl_histogram_set_ranges_uniform (h, thetamin-0.05, thetamax+0.05);

	foreach()
		gsl_histogram_accumulate(h, T[], dv());

	gsl_histogram_pdf * p = gsl_histogram_pdf_alloc (nbin);
	gsl_histogram_pdf_init (p, h);

	double p_range[nbin], p_sum[nbin], p_sum_loc[nbin];
	for (int ii = 0; ii < nbin; ii++){
		p_range[ii] = p->range[ii];
		p_sum_loc[ii] = p->sum[ii];
	}
	gsl_histogram_free (h);
	gsl_histogram_pdf_free (p);

	/* Sum the values over the ensemble of sub-domains */
	@if _MPI
		MPI_Allreduce (&p_sum_loc[0], &p_sum[0], nbin, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
	@else
		for (int ii = 0; ii < nbin; ii++)
			p_sum[ii] = p_sum_loc[ii];
	@endif

	/* Use the CDF to compute the reference height Yref(x,y,z,t) */
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, nbin);
  gsl_spline_init(spline, p_range, p_sum, nbin);
	foreach()
  	yref[] = gsl_spline_eval (spline, T[], acc) + Y0;
	gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

	/* Save the PDF and CDF (optional)*/
	// FILE * fp = fopen("reference_state_gsl.asc", "w");
	// for (int i = 0; i < nbin; i++)
	// 	fprintf (fp, "%g %g \n", p_range[i], p_sum[i]);
	// fclose (fp);
}
