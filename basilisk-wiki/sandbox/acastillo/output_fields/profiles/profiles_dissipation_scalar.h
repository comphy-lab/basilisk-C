double average_dissipation_scalar_plane(scalar c, scalar w, double *mean_grad, double *mean_grad_sq, double *mean_chi, double hprof = 0., coord box[2] = {{X0, hprof}, {X0 + L0, hprof}}, coord nsamples = {N, 1})
{
  #if dimension == 2
    int len_grad = 2;
  #else
    int len_grad = 3;
  #endif

  double sample_count = 0;  // Counter for points within the y-range
  double total_weight = 0;  // Accumulator for weights

  coord p;

  // Main loop: Iterate over the specified region
  foreach_region(p, box, nsamples, reduction(+ : total_weight) reduction(+ : sample_count) reduction(+ : mean_grad[:len_grad]) reduction(+ : mean_grad_sq[:len_grad]) reduction(+ : mean_chi[0:1]) ){
    sample_count++;

    // Compute weight
    double w_factor = (w.i == unity.i) ? 1. : w[];
    total_weight += w_factor;

    // Compute centred-difference gradients of c
    double dcdx = (c[1]   - c[-1]  )/(2.*Delta);
    double dcdy = (c[0,1] - c[0,-1])/(2.*Delta);
    #if dimension == 3
      double dcdz = (c[0,0,1] - c[0,0,-1])/(2.*Delta);
    #endif

    // Scalar dissipation rate: chi = |grad c|^2
    double chi = sq(dcdx) + sq(dcdy);
    #if dimension == 3
      chi += sq(dcdz);
    #endif

    // Accumulate gradient components
    mean_grad[0] += dcdx * w_factor;
    mean_grad[1] += dcdy * w_factor;
    #if dimension == 3
      mean_grad[2] += dcdz * w_factor;
    #endif

    // Accumulate squared gradient components
    mean_grad_sq[0] += sq(dcdx) * w_factor;
    mean_grad_sq[1] += sq(dcdy) * w_factor;
    #if dimension == 3
      mean_grad_sq[2] += sq(dcdz) * w_factor;
    #endif

    // Accumulate scalar dissipation
    mean_chi[0] += chi * w_factor;
  }

  // Normalize accumulated values
  for (int g = 0; g < len_grad; g++) {
    mean_grad[g]    /= total_weight;
    mean_grad_sq[g] /= total_weight;
  }
  mean_chi[0] /= total_weight;

  // Return average cell size or its square root, depending on dimension
#if (dimension == 2)
  return total_weight / (double)sample_count;
#else
  return sqrt(total_weight / (double)sample_count);
#endif
}

#define _mindelta (L0 / (double)(1 << (grid->maxdepth + 1)))
void profile_dissipation_scalar_foreach_region(scalar c,
                            scalar w = unity,
                            char *filename = "profiles_chi.asc",
                            double hmin = Z0 + _mindelta,
                            double hmax = Z0 + L0 - _mindelta,
                            double xmin = X0,
                            double xmax = X0 + L0,
                            double ymin = Y0,
                            double ymax = Y0 + L0,
                            double rf = 1,
                            int n = N,
                            int m1 = N,
                            int m2 = N,
                            const char *mode = "a") {

  double deltahn = (hmax - hmin) / ((double)n - 0.99999999);

  #if dimension == 2
    int len_grad = 2;
  #else
    int len_grad = 3;
  #endif

  FILE *fp = NULL;
  // The primary worker (pid 0) handles file writing
  if (pid() == 0) {
    fp = fopen(filename, mode);
    if (fp == NULL) {
      perror(filename);
      exit(1);
    }

    // Write header
    fprintf(fp, "#y@%g\t", t);
    #if dimension == 2
      fprintf(fp, "dcdx\tdcdy\tdcdx2\tdcdy2\tchi\n");
    #else
      fprintf(fp, "dcdx\tdcdy\tdcdz\tdcdx2\tdcdy2\tdcdz2\tchi\n");
    #endif
  }

  double hprof = hmin;

  // Iterate over different y-coordinates (in 2D) or z-coordinates (in 3D)
  while (hprof <= hmax) {
    // Define the region of interest
    #if dimension == 2
      coord box[2] = {{xmin, hprof}, {xmax, hprof}};
      coord nsamples = {m1, 1};
    #else
      coord box[2] = {{xmin, ymin, hprof}, {xmax, ymax, hprof}};
      coord nsamples = {m1, m2, 1};
    #endif

    // Allocate and zero accumulation arrays
    double mean_grad[len_grad], mean_grad_sq[len_grad], mean_chi[1];
    memset(&mean_grad,    0, sizeof(mean_grad));
    memset(&mean_grad_sq, 0, sizeof(mean_grad_sq));
    mean_chi[0] = 0.;

    double deltah = average_dissipation_scalar_plane(c, w, mean_grad, mean_grad_sq, mean_chi, hprof, box, nsamples);

    // Write results to file (primary worker only)
    if (pid() == 0) {
      fprintf(fp, "%.17g", hprof);
      for (int g = 0; g < len_grad; g++)
        fprintf(fp, "\t%.17g", mean_grad[g]);
      for (int g = 0; g < len_grad; g++)
        fprintf(fp, "\t%.17g", mean_grad_sq[g]);
      fprintf(fp, "\t%.17g\n", mean_chi[0]);
    }

    // Calculate next y- (or z-) coordinate
    deltah = deltahn / rf;
    hprof += rf * deltah;
  }

  // Close file (primary worker only)
  if (pid() == 0) {
    fprintf(fp, "\n");
    fflush(fp);
    if (fp != stdout)
      fclose(fp);
  }
}
#undef _mindelta