double average_dissipation_plane(vector u, scalar w, double *mean_grad, double *mean_grad_sq, double *mean_strain_sq, double hprof = 0., coord box[2] = {{X0, hprof}, {X0 + L0, hprof}}, coord nsamples = {N, 1})
{
  #if dimension == 2
    int len1 = 4, len2 = 4, len3 = 4;
  #else
    int len1 = 9, len2 = 9, len3 = 7;
  #endif

  double sample_count = 0;  // Counter for points within the y-range
  double total_weight = 0;  // Accumulator for weights

  coord p;

  // Main loop: Iterate over the specified region
  foreach_region(p, box, nsamples, reduction(+ : total_weight) reduction(+ : sample_count) reduction(+ : mean_grad[:len1]) reduction(+ : mean_grad_sq[:len2]) reduction(+ : mean_strain_sq[:len3]) ){
    sample_count++;

    // Compute gradients: either of u (unweighted) or of phi = w*u (weighted momentum)
    double dudx, dvdx, dudy, dvdy;
    if (w.i == unity.i) {
      total_weight += 1.;
      dudx = (u.x[1]   - u.x[-1]  )/(2.*Delta);
      dvdx = (u.y[1]   - u.y[-1]  )/(2.*Delta);
      dudy = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
      dvdy = (u.y[0,1] - u.y[0,-1])/(2.*Delta);
    } else {
      total_weight += w[];
      dudx = (w[1,0]*u.x[1,0] - w[-1,0]*u.x[-1,0])/(2.*Delta);
      dvdx = (w[1,0]*u.y[1,0] - w[-1,0]*u.y[-1,0])/(2.*Delta);
      dudy = (w[0,1]*u.x[0,1] - w[0,-1]*u.x[0,-1])/(2.*Delta);
      dvdy = (w[0,1]*u.y[0,1] - w[0,-1]*u.y[0,-1])/(2.*Delta);
    }
    #if dimension == 3
      double dwdx, dwdy, dudz, dvdz, dwdz;
      if (w.i == unity.i) {
        dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
        dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
        dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
        dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
        dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
      } else {
        dwdx = (w[1,0,0]*u.z[1,0,0] - w[-1,0,0]*u.z[-1,0,0])/(2.*Delta);
        dwdy = (w[0,1,0]*u.z[0,1,0] - w[0,-1,0]*u.z[0,-1,0])/(2.*Delta);
        dudz = (w[0,0,1]*u.x[0,0,1] - w[0,0,-1]*u.x[0,0,-1])/(2.*Delta);
        dvdz = (w[0,0,1]*u.y[0,0,1] - w[0,0,-1]*u.y[0,0,-1])/(2.*Delta);
        dwdz = (w[0,0,1]*u.z[0,0,1] - w[0,0,-1]*u.z[0,0,-1])/(2.*Delta);
      }
    #endif

    // Strain-rate tensor (symmetric, so Syx=Sxy, Szx=Sxz, Szy=Syz)
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

    // Accumulate
    mean_grad[0] += dudx;
    mean_grad[1] += dvdx;
    mean_grad[2] += dudy;
    mean_grad[3] += dvdy;
    #if dimension == 3
      mean_grad[4] += dwdx;
      mean_grad[5] += dwdy;
      mean_grad[6] += dudz;
      mean_grad[7] += dvdz;
      mean_grad[8] += dwdz;
    #endif

    mean_grad_sq[0] += sq(dudx);
    mean_grad_sq[1] += sq(dvdx);
    mean_grad_sq[2] += sq(dudy);
    mean_grad_sq[3] += sq(dvdy);
    #if dimension == 3
      mean_grad_sq[4] += sq(dwdx);
      mean_grad_sq[5] += sq(dwdy);
      mean_grad_sq[6] += sq(dudz);
      mean_grad_sq[7] += sq(dvdz);
      mean_grad_sq[8] += sq(dwdz);
    #endif

    mean_strain_sq[0] += sq(S2);
    mean_strain_sq[1] += sq(Sxx);
    mean_strain_sq[2] += sq(Sxy);
    mean_strain_sq[3] += sq(Syy);
    #if dimension == 3
      mean_strain_sq[4] += sq(Szz);
      mean_strain_sq[5] += sq(Sxz);
      mean_strain_sq[6] += sq(Syz);
    #endif
  }

  // Normalize accumulated values
  for (int g = 0; g < len1; g++) mean_grad[g]      /= total_weight;
  for (int g = 0; g < len2; g++) mean_grad_sq[g]   /= total_weight;
  for (int g = 0; g < len3; g++) mean_strain_sq[g] /= total_weight;

  // Return average cell size or its square root, depending on dimension
#if (dimension == 2)
  return total_weight / (double)sample_count;
#else
  return sqrt(total_weight / (double)sample_count);
#endif
}

#define _mindelta (L0 / (double)(1 << (grid->maxdepth + 1)))
void profile_dissipation_foreach_region(vector *vlist,
                            scalar w = unity,
                            char *filename = "profiles.asc",                                  
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
    fprintf(fp, "\n");
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

    // Calculate averages for the current region
    #if dimension == 2
      int len1 = 4, len2 = 4, len3 = 4;
    #else
      int len1 = 9, len2 = 9, len3 = 7;
    #endif
    double mean_grad[len1], mean_grad_sq[len2], mean_strain_sq[len3];
    memset(&mean_grad,      0, sizeof(mean_grad));
    memset(&mean_grad_sq,   0, sizeof(mean_grad_sq));
    memset(&mean_strain_sq, 0, sizeof(mean_strain_sq));

    double deltah = average_dissipation_plane(vlist, w, mean_grad, mean_grad_sq, mean_strain_sq, hprof, box, nsamples);

    // Write results to file (primary worker only)
    if (pid() == 0) {
      for (int k = 0; k < len1; k++) {
        if (k == 0)
          fprintf(fp, "%.17g\t%.17g", hprof, mean_grad[k]);
        else
          fprintf(fp, "\t%.17g", mean_grad[k]);
      }
      for (int k = 0; k < len2; k++)
        fprintf(fp, "\t%.17g", mean_grad_sq[k]);
      for (int k = 0; k < len3; k++)
        fprintf(fp, "\t%.17g", mean_strain_sq[k]);
      fprintf(fp, "\n");
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