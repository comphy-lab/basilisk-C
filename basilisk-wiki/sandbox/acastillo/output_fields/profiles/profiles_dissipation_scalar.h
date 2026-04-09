double average_dissipation_scalar_plane(scalar c, scalar w, double *mean_grad, double *mean_grad_sq, double mean_chi, double hprof, coord box[2], coord nsamples)
{
  #if dimension == 2
    int len_grad = 2;
  #else
    int len_grad = 3;
  #endif

  double sample_count = 0;
  double total_weight = 0;
  coord p;

  // Main loop: Iterate over the specified region
  foreach_region(p, box, nsamples, reduction(+ : total_weight) reduction(+ : sample_count) reduction(+ : mean_grad[:len_grad]) reduction(+ : mean_grad_sq[:len_grad]) reduction(+ : mean_chi) ){
    sample_count++;
    double weight = (w.i != unity.i) ? w[] : 1.;
    total_weight += weight;

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

    // Accumulate weighted gradients and scalar dissipation
    mean_grad[0] += dcdx * weight;
    mean_grad[1] += dcdy * weight;
    #if dimension == 3
      mean_grad[2] += dcdz * weight;
    #endif

    mean_grad_sq[0] += sq(dcdx) * weight;
    mean_grad_sq[1] += sq(dcdy) * weight;
    #if dimension == 3
      mean_grad_sq[2] += sq(dcdz) * weight;
    #endif

    mean_chi += chi * weight;
  }

  for (int g = 0; g < len_grad; g++) mean_grad[g]    /= total_weight;
  for (int g = 0; g < len_grad; g++) mean_grad_sq[g] /= total_weight;
  mean_chi /= total_weight;

  #if (dimension == 2)
    return total_weight / (double)sample_count;
  #else
    return sqrt(total_weight / (double)sample_count);
  #endif
}

void profile_dissipation_scalar_foreach_region(scalar c, PROFILE_PARAMS) {

  double deltahn = (hmax - hmin) / ((double)n - 0.99999999);
  #if dimension == 2
    int len_grad = 2;
  #else
    int len_grad = 3;
  #endif

  FILE *fp = NULL;
  if (pid() == 0) {
    fp = fopen(filename, mode);
    if (fp == NULL) { perror(filename); exit(1); }

    fprintf(fp, "# Profile Scalar Dissipation: t = %.10g, L0 = %g\n", t, L0);
    fprintf(fp, "# [0]iprof [1]y [2]delta");
    int k = 3;
    for (int g = 0; g < len_grad; g++) fprintf(fp, " [%d]grad", k++);
    for (int g = 0; g < len_grad; g++) fprintf(fp, " [%d]grad_sq", k++);
    fprintf(fp, " [%d]chi", k++);
    fputc('\n', fp);
  }

  int iprof = 0;
  double hprof = hmin;

  // Iterate over different y-coordinates (in 2D) or z-coordinates (in 3D)
  while (hprof <= hmax) {
    SETUP_PROFILE_PLANE(hprof)

    double mean_grad[len_grad], mean_grad_sq[len_grad], mean_chi;
    memset(mean_grad,    0, sizeof(mean_grad));
    memset(mean_grad_sq, 0, sizeof(mean_grad_sq));
    mean_chi = 0.;

    double deltah = average_dissipation_scalar_plane(c, w, mean_grad, mean_grad_sq, mean_chi, hprof, box, nsamples);

    // Write results to file (primary worker only)
    if (pid() == 0) {
      fprintf(fp, "%-6d %15.8e %15.8e", iprof, hprof, deltah);
      for (int k = 0; k < len_grad; k++) fprintf(fp, " %24.15e", mean_grad[k]);
      for (int k = 0; k < len_grad; k++) fprintf(fp, " %24.15e", mean_grad_sq[k]);
      fprintf(fp, " %24.15e", mean_chi);
      fputc('\n', fp);
    }

    // Calculate next y- (or z-) coordinate
    hprof += deltahn;
    iprof++;
  }

  if (pid() == 0) {
    fputc('\n', fp);
    fputc('\n', fp);
    fflush(fp);
    if (fp != stdout) fclose(fp);
  }
}