double average_dissipation_plane(vector u, scalar w, double *mean_grad, double *mean_grad_sq, double *mean_strain_sq, double hprof, coord box[2], coord nsamples)
{
  #if dimension == 2
    int len1 = 4, len2 = 4;
  #else
    int len1 = 9, len2 = 7;
  #endif

  double sample_count = 0;
  double total_weight = 0;
  coord p;

  foreach_region(p, box, nsamples, reduction(+ : total_weight) reduction(+ : sample_count) reduction(+ : mean_grad[:len1]) reduction(+ : mean_grad_sq[:len1]) reduction(+ : mean_strain_sq[:len2]) ){
    sample_count++;
    double weight = (w.i != unity.i) ? w[] : 1.;
    total_weight += weight;

    double dudx, dvdx, dudy, dvdy;
    dudx = (u.x[1]   - u.x[-1]  )/(2.*Delta);
    dvdx = (u.y[1]   - u.y[-1]  )/(2.*Delta);
    dudy = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    dvdy = (u.y[0,1] - u.y[0,-1])/(2.*Delta);

    #if dimension == 3
      double dwdx, dwdy, dudz, dvdz, dwdz;
      dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
      dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
      dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
      dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
      dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    #endif

    // Strain-rate tensor
    double Sxx = dudx;
    double Sxy = 0.5*(dudy + dvdx);
    double Syy = dvdy;
    double S2  = sq(Sxx) + 2.*sq(Sxy) + sq(Syy);
    #if dimension == 3
      double Szz = dwdz;
      double Sxz = 0.5*(dwdx + dudz);
      double Syz = 0.5*(dwdy + dvdz);
      S2 += sq(Szz) + 2.*sq(Sxz) + 2.*sq(Syz);
    #endif

    // Accumulate weighted gradients and strains
    mean_grad[0] += dudx * weight;
    mean_grad[1] += dvdx * weight;
    mean_grad[2] += dudy * weight;
    mean_grad[3] += dvdy * weight;
    #if dimension == 3
      mean_grad[4] += dwdx * weight;
      mean_grad[5] += dwdy * weight;
      mean_grad[6] += dudz * weight;
      mean_grad[7] += dvdz * weight;
      mean_grad[8] += dwdz * weight;
    #endif

    mean_grad_sq[0] += sq(dudx) * weight;
    mean_grad_sq[1] += sq(dvdx) * weight;
    mean_grad_sq[2] += sq(dudy) * weight;
    mean_grad_sq[3] += sq(dvdy) * weight;
    #if dimension == 3
      mean_grad_sq[4] += sq(dwdx) * weight;
      mean_grad_sq[5] += sq(dwdy) * weight;
      mean_grad_sq[6] += sq(dudz) * weight;
      mean_grad_sq[7] += sq(dvdz) * weight;
      mean_grad_sq[8] += sq(dwdz) * weight;
    #endif

    #ifndef mu
      mean_strain_sq[0] += S2 * weight;
    #else
      mean_strain_sq[0] += (mu(f[])/rhov[]) * S2 * weight;
    #endif 
    mean_strain_sq[1] += sq(Sxx) * weight;
    mean_strain_sq[2] += sq(Sxy) * weight;
    mean_strain_sq[3] += sq(Syy) * weight;
    #if dimension == 3
      mean_strain_sq[4] += sq(Szz) * weight;
      mean_strain_sq[5] += sq(Sxz) * weight;
      mean_strain_sq[6] += sq(Syz) * weight;
    #endif
  }

  for (int g = 0; g < len1; g++) mean_grad[g]      /= total_weight;
  for (int g = 0; g < len1; g++) mean_grad_sq[g]   /= total_weight;
  for (int g = 0; g < len2; g++) mean_strain_sq[g] /= total_weight;

  #if (dimension == 2)
    return total_weight / (double)sample_count;
  #else
    return sqrt(total_weight / (double)sample_count);
  #endif
}

void profile_dissipation_foreach_region(vector v, PROFILE_PARAMS) {

  double deltahn = (hmax - hmin) / ((double)n - 0.99999999);
  #if dimension == 2
    int len1 = 4, len2 = 4;
  #else
    int len1 = 9, len2 = 7;
  #endif

  FILE *fp = NULL;
  if (pid() == 0) {
    fp = fopen(filename, mode);
    if (fp == NULL) { perror(filename); exit(1); }

    fprintf(fp, "# Profile Dissipation: t = %.10g, L0 = %g\n", t, L0);
    fprintf(fp, "# [0]iprof [1]y [2]delta");
    int k = 3;
    for (int g = 0; g < len1; g++) fprintf(fp, " [%d]grad", k++);
    for (int g = 0; g < len1; g++) fprintf(fp, " [%d]grad_sq", k++);
    for (int g = 0; g < len2; g++) fprintf(fp, " [%d]strain_sq", k++);
    fputc('\n', fp);
  }

  int iprof = 0;
  double hprof = hmin;

  // Iterate over different y-coordinates (in 2D) or z-coordinates (in 3D)
  while (hprof <= hmax) {
    SETUP_PROFILE_PLANE(hprof)

    // Calculate averages for the current region

    double mean_grad[len1], mean_grad_sq[len1], mean_strain_sq[len2];
    memset(mean_grad,      0, sizeof(mean_grad));
    memset(mean_grad_sq,   0, sizeof(mean_grad_sq));
    memset(mean_strain_sq, 0, sizeof(mean_strain_sq));

    double deltah = average_dissipation_plane(v, w, mean_grad, mean_grad_sq, mean_strain_sq, hprof, box, nsamples);

    // Write results to file (primary worker only)
    if (pid() == 0) {
      fprintf(fp, "%-6d %15.8e %15.8e", iprof, hprof, deltah);
      for (int k = 0; k < len1; k++) fprintf(fp, " %24.15e", mean_grad[k]);
      for (int k = 0; k < len1; k++) fprintf(fp, " %24.15e", mean_grad_sq[k]);
      for (int k = 0; k < len2; k++) fprintf(fp, " %24.15e", mean_strain_sq[k]);
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
