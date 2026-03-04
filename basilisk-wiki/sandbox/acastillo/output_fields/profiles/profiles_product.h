double average_product_plane(scalar *list1, scalar *list2, scalar w, double *averages, double hprof = 0., coord box[2] = {{X0, hprof}, {X0 + L0, hprof}}, coord nsamples = {N, 1})
{

  int len = list_len(list1); // Get length of the scalar list
  double sample_count = 0;  // Counter for points within the y-range
  double total_weight = 0 ; // Accumulator for weights
  
  coord p;
  NOT_UNUSED(len);


  // Main loop: Iterate over the specified region
  foreach_region(p, box, nsamples, reduction(+ : total_weight) reduction(+ : sample_count) reduction(+ : averages[:len]) ){
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
    scalar s1, s2;
    if (w.i == unity.i){ // If no weights are provided    
      for (s1, s2 in list1, list2){
        double val = s1[]*s2[];
        averages[k] += val * w_factor;
        k++;
      }
    }
    else {// If weights are provided  
      for (s1, s2 in list1, list2){
        double val = s1[]*s2[];
        double weight = w[] * w_factor;
        averages[k] += val * weight;
        k++;
      }
    }
  }

  // Normalize accumulated values
  for (int g = 0; g < len; g++){
    averages[g] /= total_weight;
  }

  // Return average cell size or its square root, depending on dimension
#if (dimension == 2)
  return total_weight / (double)sample_count;
#else
  return sqrt(total_weight / (double)sample_count);
#endif
}


void profile_product_foreach_region(scalar *list1 = all,
                                    scalar *list2 = all,
                                    scalar w = unity,
                                    char *filename = "profiles_prod.asc",
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

  int len1 = list_len(list1);
  int len2 = list_len(list2);
  if (len1 != len2) {
    fprintf(stderr, "profile_product_foreach_region: list lengths must match (%d vs %d)\n", len1, len2);
    return;
  }

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
    scalar s1, s2;
    for (s1, s2 in list1, list2)
      fprintf(fp, "%s_%s\t", s1.name, s2.name);
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
    double aver[len1];
    memset(&aver, 0, sizeof(aver));

    double deltah = average_product_plane(list1, list2, w, aver, hprof, box, nsamples);

    // Write results to file (primary worker only)
    if (pid() == 0) {
      for (int k = 0; k < len1; k++) {
        if (k == 0) {
          fprintf(fp, "%.17g\t%.17g", hprof, aver[k]);
        }
        else {
          fprintf(fp, "\t%.17g", aver[k]);
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