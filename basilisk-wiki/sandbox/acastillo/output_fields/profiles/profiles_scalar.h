double average_scalar_plane(scalar *list, scalar w, double *averages, double *averages_sq, double *averages_cb, double hprof = 0., coord box[2] = {{X0, hprof}, {X0 + L0, hprof}}, coord nsamples = {N, 1})
{

  int len = list_len(list); // Get length of the scalar list
  double sample_count = 0;  // Counter for points within the y-range
  double total_weight = 0 ; // Accumulator for weights
  
  coord p;
  NOT_UNUSED(len);

  // Main loop: Iterate over the specified region
  foreach_region(p, box, nsamples, reduction(+ : total_weight) reduction(+ : sample_count) reduction(+ : averages[:len]) reduction(+ : averages_sq[:len]) reduction(+ : averages_cb[:len]) ){
    sample_count++;

    // Calculate weight factor 
    double w_factor = 1.;
    
    if (w.i == unity.i){ // If no weights are provided    
      total_weight += w_factor;
    }
    else { // If weights are provided  
      total_weight += w[] * w_factor;
    }
    
    // Accumulate weighted values for each scalar field      
    int k = 0;
    if (w.i == unity.i){ // If no weights are provided    
      for (scalar s in list){
        double val = s[];
        averages[k] += val * w_factor;
        averages_sq[k] += sq(val) * w_factor;
        averages_cb[k] += cube(val) * w_factor;
        k++;
      }
    }
    else {// If weights are provided  
      for (scalar s in list){
        double val = s[];
        double weight = w[] * w_factor;
        averages[k] += val * weight;
        averages_sq[k] += sq(val) * weight;
        averages_cb[k] += cube(val) * weight;
        k++;
      }
    }
  }

  // Normalize accumulated values
  for (int g = 0; g < len; g++){
    averages[g] /= total_weight;
    averages_sq[g] /= total_weight;
    averages_cb[g] /= total_weight;
  }

  // Return average cell size or its square root, depending on dimension
#if (dimension == 2)
  return total_weight / (double)sample_count;
#else
  return sqrt(total_weight / (double)sample_count);
#endif
}

#define _mindelta (L0 / (double)(1 << (grid->maxdepth + 1)))
void profile_foreach_region(scalar *list = all,                                  
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
  int len = list_len(list);
  
  // The primary worker (pid 0) handles file writing
  if (pid() == 0) {
    fp = fopen(filename, mode);
    if (fp == NULL) {
      perror(filename);
      exit(1);
    }

    // Write header
    fprintf(fp, "#y@%g\t", t);
    for (scalar s in list)
      fprintf(fp, "%s\t%s.m2\t", s.name, s.name);
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
    double aver[len], aver_sq[len], aver_cb[len];
    memset(&aver, 0, sizeof(aver));
    memset(&aver_sq, 0, sizeof(aver_sq));
    memset(&aver_cb, 0, sizeof(aver_cb));

    double deltah = average_scalar_plane(list, w, aver, aver_sq, aver_cb, hprof, box, nsamples);

    // Write results to file (primary worker only)
    if (pid() == 0) {
      for (int k = 0; k < len; k++) {
        double mean = aver[k];
        double m2 = aver_sq[k];
        double m3 = aver_cb[k];
        
        if (k == 0) {
          fprintf(fp, "%.17g\t%.17g\t%.17g\t%.17g", hprof, mean, m2, m3);
        }
        else {
          fprintf(fp, "\t%.17g\t%.17g\t%.17g", mean, m2, m3);
        }
      }
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