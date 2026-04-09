double average_product_plane(scalar *list1, scalar *list2, scalar w, double *averages, double hprof, coord box[2], coord nsamples)
{
  int len = list_len(list1);
  double sample_count = 0;
  double total_weight = 0;
  coord p;
  NOT_UNUSED(len);

  foreach_region(p, box, nsamples, reduction(+ : total_weight) reduction(+ : sample_count) reduction(+ : averages[:len]) ){
    sample_count++;
    double weight = (w.i != unity.i) ? w[] : 1.;
    total_weight += weight;

    int k = 0;
    scalar s1, s2;
    for (s1, s2 in list1, list2){
      averages[k++] += s1[] * s2[] * weight;
    }
  }

  for (int g = 0; g < len; g++)
    averages[g] /= total_weight;

  #if (dimension == 2)
    return total_weight / (double)sample_count;
  #else
    return sqrt(total_weight / (double)sample_count);
  #endif
}

void profile_product_foreach_region(scalar *list1 = all, scalar *list2 = all, PROFILE_PARAMS) {

  int len1 = list_len(list1);
  int len2 = list_len(list2);
  if (len1 != len2) {
    fprintf(stderr, "profile_product_foreach_region: list lengths must match (%d vs %d)\n", len1, len2);
    return;
  }

  double deltahn = (hmax - hmin) / ((double)n - 0.99999999);
  FILE *fp = NULL;

  if (pid() == 0) {
    fp = fopen(filename, mode);
    if (fp == NULL) { perror(filename); exit(1); }

    fprintf(fp, "# Profile Product: t = %.10g, L0 = %g\n", t, L0);
    fprintf(fp, "# [0]iprof [1]y [2]delta");
    int k = 3;
    scalar s1, s2;
    for (s1, s2 in list1, list2)
      fprintf(fp, " [%d]mean(%s*%s)", k++, s1.name, s2.name);
    fputc('\n', fp);
  }

  int iprof = 0;
  double hprof = hmin;
  while (hprof <= hmax) {
    SETUP_PROFILE_PLANE(hprof)

    double aver[len1];
    memset(aver, 0, sizeof(aver));
    double deltah = average_product_plane(list1, list2, w, aver, hprof, box, nsamples);

    if (pid() == 0) {
      fprintf(fp, "%-6d %15.8e %15.8e", iprof, hprof, deltah);
      for (int k = 0; k < len1; k++) {
        fprintf(fp, " %24.15e", aver[k]);
      }
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
