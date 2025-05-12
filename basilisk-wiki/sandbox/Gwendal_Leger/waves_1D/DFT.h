#include "dalloc.h"

// Custom type containing values of the DFT
typedef struct {
  double N; // Resolution
  double * f; // Frequencies
  double * r; // Real part
  double * i; // Imaginary part
  double * n; // Norm (euclidian)
} DFT;



void interpolate_signal (int M, int N, double *samples, double *signal, double *x, double *s) {
  int n;
  // Uniform sample of the signal (via linear interpolation), to the resolution wanted
  double dx = (samples[N-1] - samples[0])/(double)(M - 1);
  x[0] = samples[0];
  s[0] = signal[0];
  OMP_PARALLEL()
  {
    OMP(omp for schedule(static))
    for (n = 1; n < M-1; n++) {
      int k = 0;
      x[n] = n*dx;
      while (x[n] > samples[k+1])
	k++;
      double t = (x[n] - samples[k])/(samples[k+1] - samples[k]);
      s[n] = (1. - t)*signal[k] + t*signal[k+1]; // Linear interpolation
    }
  } // end parallel
  x[M-1] = x[M-2] + dx;
  s[M-1] = signal[N-1];
}


DFT discrete_Fourier_transform (int N, double sample_length, double *signal) {
  DFT S;
  int k;
  S.N = N;
  S.f = dalloc (S.N);
  S.r = dalloc (S.N);
  S.i = dalloc (S.N);
  S.n = dalloc (S.N);
  OMP_PARALLEL()
  {
    int SN = S.N;
    OMP(omp for schedule(static))
    for (k = 0; k < SN; k++) {
      S.f[k] = (double)k/sample_length;
      int n;
      for (n = 0; n < N; n++) {
	double theta = -2.*pi*(double)(n*k)/(double)N;
	S.r[k] += signal[n]*cos(theta);
	S.i[k] += signal[n]*sin(theta);
      }
      S.n[k] = sqrt (S.r[k]*S.r[k] + S.i[k]*S.i[k]);
    }
  } // end parallel
  return S;
}


double * inverse_discrete_Fourier_transform (int N, DFT S) {
  int n;
  double * signal = dalloc (N);
  OMP_PARALLEL()
  {
    OMP(omp for schedule(static))
    for (n = 0; n < N; n++) {
      int k;
      for (k = 0; k < S.N; k++) {
	double theta = 2.*pi*(double)(n*k)/(double)S.N;
	signal[n] += S.r[k]*cos(theta) - S.i[k]*sin(theta);
      }
      signal[n] /= (double)N;
    }
  } // end parallel  
  return signal;
}


DFT nonuniform_discrete_Fourier_transform (int N, double sample_length, double *samples, double *signal) {
  DFT S;
  int k;
  S.N = N;
  S.f = dalloc (S.N);
  S.r = dalloc (S.N);
  S.i = dalloc (S.N);
  S.n = dalloc (S.N);
  OMP_PARALLEL()
  {
    int SN = S.N;
    OMP(omp for schedule(static))
    for (k = 0; k < SN; k++) {
      S.f[k] = (double)k/sample_length;
      int n;
      for (n = 0; n < N; n++) {
	double theta = -2.*pi*samples[n]/sample_length*(double)k;
	S.r[k] += signal[n]*cos(theta);
	S.i[k] += signal[n]*sin(theta);
      }
      S.n[k] = sqrt (S.r[k]*S.r[k] + S.i[k]*S.i[k]);
    }
  } // end parallel
  return S;
}



// Function to read data from a file typically produced with "output_gauge" (2 columns and the first line is the title), Fourier transform said data and plot the results.
void discrete_Fourier_transform_from_file (char * filename, int M, char * title) {
  int N; // N is the signal resolution, M is the wanted resolution (<0 to interpolate with N points, 0 to keep the points, >0 to interpolate with M points)
  double * samples;
  double * signal;
  int interpolated = 0;
  int k;
  FILE * fp = fopen (filename, "r");
  // Count the number of lines in the file
  int ch = 0;
  k = 0;
  /* while (EOF != (fscanf (fp, "%*[^\n]"), fscanf (fp, "%*c"))) {k++;} */
  while (EOF != (ch = fgetc (fp))) {
    if (ch == '\n')
      k++;
  }
  fclose (fp);

  N = k - 1;
  if (M < 0)
    M = N;
  if (title == NULL)
    title = filename;
  
  fprintf (stdout, "Discrete Fourier transform of %s (%s) signal with : M = %i, N = %i\n", filename, title, M, N);
  
  samples = dalloc (N);
  signal = dalloc (N);

  fp = fopen (filename, "r");
  char buffer[256];
  fgets (buffer, sizeof(buffer), fp); // Jumping the first line
  for (k = 0; k < N; k++)
    fscanf (fp, "%lg %lg\n", &samples[k], &signal[k]);
  fclose (fp);

    
  // Interpolation of the signal
  double * samples_interpolated;
  double * signal_interpolated;
  if (M != 0) { // Interpolate
    interpolated = 1;
    samples_interpolated = dalloc (M);
    signal_interpolated = dalloc (M);
    interpolate_signal (M, N, samples, signal, samples_interpolated, signal_interpolated);
  } else { // Use the discrete signal as it is
    interpolated = 0;
    samples_interpolated = dalloc (N);
    signal_interpolated = dalloc (N);
    for (k = 0; k < N; k++) {
      samples_interpolated[k] = samples[k];
      signal_interpolated[k] = signal[k];
    }
  }
  
  // Fourier transform
  DFT signal_Fourier_transformed;
  if (interpolated)
    signal_Fourier_transformed = discrete_Fourier_transform (M, samples_interpolated[M-1], signal_interpolated);
  else
    signal_Fourier_transformed = nonuniform_discrete_Fourier_transform (N, samples[N-1], samples, signal);
  
  // Inverse Fourier transform
  double * inverse_signal_Fourier_transformed;
  if (interpolated)
    inverse_signal_Fourier_transformed = inverse_discrete_Fourier_transform (M, signal_Fourier_transformed);
  else {
    inverse_signal_Fourier_transformed = dalloc (N);
    M = N;
  }
  
  // Write data to files and a gnuplot script
  //#pragma omp parallel default (none) shared (filename, M, N, samples, signal, samples_interpolated, signal_interpolated, signal_Fourier_transformed, inverse_signal_Fourier_transformed, interpolated) private (k, fp)
  {
    //#pragma omp sections nowait
    {
      //#pragma omp section
      { // Writing signal data
	if (interpolated) {
	  char fn[128];
	  sprintf (fn, "%s_signal", filename);
	  fp = fopen (fn, "w");
	  for (k = 0; k < N; k++)
	    fprintf (fp, "%g %g\n", samples[k], signal[k]);
	  fclose (fp);
	}
      }
      //#pragma omp section
      { // Writing the interpolated signal and its Fourier transform
	char fn[128];
	sprintf (fn, "%s_out", filename);
	fp = fopen (fn, "w");
	for (k = 0; k < M; k++)
	  fprintf (fp, "%g %g %g %g %g\n",
		   samples_interpolated[k], signal_interpolated[k],
		   signal_Fourier_transformed.f[k], signal_Fourier_transformed.n[k],
		   inverse_signal_Fourier_transformed[k]);
	fclose (fp);	
      }
      //#pragma omp section
      { // Writing the gnuplot script
	fp = fopen ("plot.plot", "w");
	if (interpolated) {
	  fprintf (fp,
		   "reset session\n"
		   "set term pdfcairo\n"
		   "set output '%s_data.pdf'\n"
		   "set title '%s'\n"
		   "set tics nomirror\n"
		   "file_signal = '%s_signal'\n"
		   "file_out = '%s_out'\n"
		   "plot file_signal using 1:2 with linespoints lc rgb 'red' title 'Sample data', file_out using 1:2 with points lc rgb 'green' title 'interpolation'\n",
		   filename, title, filename, filename);
	}
	fprintf (fp,
		 "reset session\n"
		 "set term pdfcairo\n"
		 "set output '%s_DFT.pdf'\n"
		 "set title '%s'\n"
		 "set autoscale noextend\n"
		 "set tics nomirror\n"
		 "file_out = '%s_out'\n"
		 "set multiplot layout 2,1\n"
		 "set xlabel 'Time (s)'\n"
		 "set ylabel 'Amplitude (m)'\n"
		 "plot file_out using 1:2 with lines lc rgb 'red' title 's'\n"
		 "unset title\n"
		 "set xlabel 'Frequency (Hz)'\n"
		 "set logscale x\n"
		 //"set logscale y\n"
		 "set ylabel 'Magnitude'\n"
		 "unset y2label\n"
		 "unset y2tics\n"
		 "plot file_out every ::0::%i using 3:4 with lines lc rgb 'red' title '||S||_2'\n"
		 "unset multiplot\n",
		 filename, title, filename, (int)(M*0.6));
	fclose (fp);
      }
    } // end sections
  } // end parallel
  
  system ("gnuplot plot.plot");
  
  free (samples);
  free (signal);
  free (samples_interpolated);
  free (signal_interpolated);
  free (signal_Fourier_transformed.f);
  free (signal_Fourier_transformed.r);
  free (signal_Fourier_transformed.i);
  free (signal_Fourier_transformed.n);
  free (inverse_signal_Fourier_transformed);
}



