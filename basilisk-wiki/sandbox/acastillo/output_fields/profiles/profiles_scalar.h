double average_scalar_plane(scalar *list, scalar w, double *averages, double *averages_sq, double hprof, coord box[2], coord nsamples){
  
  int len = list_len(list);   // Get length of the scalar list
  double sample_count = 0;    // Counter for points within the y-range
  double total_weight = 0 ;   // Accumulator for weights
  coord p;
  NOT_UNUSED(len);

  // Main loop: Iterate over the specified region
  foreach_region(p, box, nsamples, reduction(+ : total_weight) reduction(+ : sample_count) reduction(+ : averages[:len]) reduction(+ : averages_sq[:len]) ){
    sample_count++;
    double weight = 1.;
    if (w.i != unity.i) 
      weight = w[];

    total_weight += weight;
    int k = 0;
    for (scalar s in list){
      double val = s[];
      averages[k] += val * weight;
      averages_sq[k] += sq(val) * weight;
      k++;
    }
  }

  // Normalize
  for (int g = 0; g < len; g++){
    averages[g] /= total_weight;
    averages_sq[g] /= total_weight;
  }

  // Return average cell size
  #if (dimension == 2)
    return total_weight / (double)sample_count;
  #else
    return sqrt(total_weight / (double)sample_count);
  #endif
}

void profile_foreach_region(scalar *list = all, PROFILE_PARAMS) {

  int len = list_len(list);
  double deltahn = (hmax - hmin) / ((double)n - 0.99999999);
  
  FILE * fp = NULL; 
  if (pid() == 0) { 
    fp = fopen(filename, mode); 
    if (fp == NULL) { perror(filename); exit(1); } 
  }

  // Write header 
  if (pid() == 0) {
    fprintf(fp, "# Profile: t = %.10g, L0 = %g\n", t, L0);
    fprintf(fp, "# [0]iprof [1]y [2]delta");
    int k = 3;
    for (scalar s in list){
      fprintf(fp, "[%d]mean(%s)\t[%d]mean(%s^2)\t", k, s.name, k+1, s.name);
      k += 2;
    }
    fputc('\n', fp);
  }

  // Iterate over different y-coordinates (in 2D) or z-coordinates (in 3D)
  int iprof = 0;
  double hprof = hmin;
  while (hprof <= hmax) {
    
    SETUP_PROFILE_PLANE(hprof)

    double aver[len], aver_sq[len];
    memset(aver, 0, sizeof(aver));
    memset(aver_sq, 0, sizeof(aver_sq));    
    double deltah = average_scalar_plane(list, w, aver, aver_sq, hprof, box, nsamples);
    
    if (pid() == 0) {
      fprintf(fp, "%-6d %15.8e %15.8e", iprof, hprof, deltah);
      
      for (int k = 0; k < len; k++) 
        fprintf(fp, " %24.15e %24.15e", aver[k], aver_sq[k]);
      
      fputc('\n', fp);
    }
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