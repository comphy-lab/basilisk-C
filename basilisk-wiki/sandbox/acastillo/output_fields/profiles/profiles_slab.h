/**
The profiles of [profiles.h](profiles.h) -- `profile_foreach_region()`,
`profile_product_foreach_region()` and `profile_dissipation_foreach_region()`
-- computed in one traversal of the grid instead of one per slab. Those
routines call `foreach_region()` for each of the `n` slabs, each with its own
reduction, so a profile costs `n` MPI collectives; here every cell is binned
into its slab by coordinate in a single `foreach()`, and the profile is
reduced once.

The two agree only where every cell of the sampled region is at the same
resolution *and* `hmin`, `hmax` and `n` are chosen so that a slab is exactly
one row of those cells. `foreach_region()` point-samples at cell centres
without interpolation, so under that condition it visits the cells this
traversal visits. Over a region spanning several levels, or with slabs thicker
than a cell, the two weight cells differently and this is not a drop-in
replacement.

Include after `navier-stokes/centered.h`, as for `profiles.h`: both need `u`
and `unity`. The output files are in the format `profiles.h` writes. */

#ifndef PROFILES_SLAB_H
#define PROFILES_SLAB_H

#define PROFILE_SLAB_PARAMS \
  scalar w = unity, \
  const char * filename = "profiles.asc", \
  double hmin = -L0/2., \
  double hmax =  L0/2., \
  double xmin = X0, \
  double xmax = X0 + L0, \
  double ymin = Y0, \
  double ymax = Y0 + L0, \
  int n = N, \
  int m1 = N, \
  int m2 = N, \
  const char * mode = "a"

#if dimension == 2
  #define PROFILE_SLAB_COORD(point) y
#else
  #define PROFILE_SLAB_COORD(point) z
#endif

/**
## Scalar profiles -- mean and mean-square of each field in `list`
*/

void profile_scalar_slab (scalar * list = all, PROFILE_SLAB_PARAMS)
{
  int len = list_len (list);
  double deltahn = (hmax - hmin)/((double) n - 0.99999999);

  double aver[n*len], aver_sq[n*len], total_weight[n], sample_count[n];
  for (int i = 0; i < n*len; i++)
    aver[i] = aver_sq[i] = 0.;
  for (int i = 0; i < n; i++)
    total_weight[i] = sample_count[i] = 0.;

  foreach (reduction(+:aver[:n*len]) reduction(+:aver_sq[:n*len])
           reduction(+:total_weight[:n]) reduction(+:sample_count[:n])) {
    double hval = PROFILE_SLAB_COORD(point);
    if (hval < hmin - Delta/2. || hval > hmax + Delta/2.) continue;
    int iprof = (int)((hval - hmin)/deltahn + 0.5);
    if (iprof < 0) iprof = 0;
    if (iprof >= n) iprof = n - 1;

    double weight = (w.i != unity.i) ? w[] : 1.;
    total_weight[iprof] += weight;
    sample_count[iprof]++;

    int k = 0;
    for (scalar s in list) {
      double val = s[];
      aver[iprof*len + k]    += val*weight;
      aver_sq[iprof*len + k] += sq (val)*weight;
      k++;
    }
  }

  if (pid() == 0) {
    FILE * fp = fopen (filename, mode);
    if (fp == NULL) { perror (filename); exit (1); }

    fprintf (fp, "# Profile: t = %.10g, L0 = %g\n", t, L0);
    fprintf (fp, "# [0]iprof [1]y [2]delta ");
    int k = 3;
    for (scalar s in list) {
      fprintf (fp, "[%d]mean(%s)\t[%d]mean(%s^2)\t", k, s.name, k + 1, s.name);
      k += 2;
    }
    fputc ('\n', fp);

    for (int iprof = 0; iprof < n; iprof++) {
      double hprof = hmin + iprof*deltahn;
      double tw = total_weight[iprof];
      double sc = sample_count[iprof];
      #if (dimension == 2)
        double deltah = sc > 0. ? tw/sc : 0.;
      #else
        double deltah = sc > 0. ? sqrt (tw/sc) : 0.;
      #endif
      fprintf (fp, "%-6d %15.8e %15.8e", iprof, hprof, deltah);
      for (int k = 0; k < len; k++) {
        double m  = tw > 0. ? aver[iprof*len + k]/tw    : 0.;
        double m2 = tw > 0. ? aver_sq[iprof*len + k]/tw : 0.;
        fprintf (fp, " %24.15e %24.15e", m, m2);
      }
      fputc ('\n', fp);
    }
    fputc ('\n', fp);
    fputc ('\n', fp);
    fflush (fp);
    fclose (fp);
  }
}

/**
## Product profiles -- mean of `s1*s2` for matching pairs in `list1`, `list2`
*/

void profile_product_slab (scalar * list1 = all, scalar * list2 = all,
                            PROFILE_SLAB_PARAMS)
{
  int len = list_len (list1);
  if (list_len (list2) != len) {
    fprintf (stderr, "profile_product_slab: list lengths must match (%d vs %d)\n",
             len, list_len (list2));
    return;
  }
  double deltahn = (hmax - hmin)/((double) n - 0.99999999);

  double aver[n*len], total_weight[n], sample_count[n];
  for (int i = 0; i < n*len; i++)
    aver[i] = 0.;
  for (int i = 0; i < n; i++)
    total_weight[i] = sample_count[i] = 0.;

  foreach (reduction(+:aver[:n*len]) reduction(+:total_weight[:n])
           reduction(+:sample_count[:n])) {
    double hval = PROFILE_SLAB_COORD(point);
    if (hval < hmin - Delta/2. || hval > hmax + Delta/2.) continue;
    int iprof = (int)((hval - hmin)/deltahn + 0.5);
    if (iprof < 0) iprof = 0;
    if (iprof >= n) iprof = n - 1;

    double weight = (w.i != unity.i) ? w[] : 1.;
    total_weight[iprof] += weight;
    sample_count[iprof]++;

    int k = 0;
    scalar s1, s2;
    for (s1, s2 in list1, list2)
      aver[iprof*len + k++] += s1[]*s2[]*weight;
  }

  if (pid() == 0) {
    FILE * fp = fopen (filename, mode);
    if (fp == NULL) { perror (filename); exit (1); }

    fprintf (fp, "# Profile Product: t = %.10g, L0 = %g\n", t, L0);
    fprintf (fp, "# [0]iprof [1]y [2]delta");
    int k = 3;
    scalar s1, s2;
    for (s1, s2 in list1, list2)
      fprintf (fp, " [%d]mean(%s*%s)", k++, s1.name, s2.name);
    fputc ('\n', fp);

    for (int iprof = 0; iprof < n; iprof++) {
      double hprof = hmin + iprof*deltahn;
      double tw = total_weight[iprof];
      double sc = sample_count[iprof];
      #if (dimension == 2)
        double deltah = sc > 0. ? tw/sc : 0.;
      #else
        double deltah = sc > 0. ? sqrt (tw/sc) : 0.;
      #endif
      fprintf (fp, "%-6d %15.8e %15.8e", iprof, hprof, deltah);
      for (int k = 0; k < len; k++) {
        double m = tw > 0. ? aver[iprof*len + k]/tw : 0.;
        fprintf (fp, " %24.15e", m);
      }
      fputc ('\n', fp);
    }
    fputc ('\n', fp);
    fputc ('\n', fp);
    fflush (fp);
    fclose (fp);
  }
}

/**
## Dissipation profiles -- velocity-gradient and strain-rate statistics
*/

void profile_dissipation_slab (vector u, PROFILE_SLAB_PARAMS)
{
  #if dimension == 2
    int len1 = 4, len2 = 4;
  #else
    int len1 = 9, len2 = 7;
  #endif
  double deltahn = (hmax - hmin)/((double) n - 0.99999999);

  double mean_grad[n*len1], mean_grad_sq[n*len1], mean_strain_sq[n*len2];
  double total_weight[n], sample_count[n];
  for (int i = 0; i < n*len1; i++)
    mean_grad[i] = mean_grad_sq[i] = 0.;
  for (int i = 0; i < n*len2; i++)
    mean_strain_sq[i] = 0.;
  for (int i = 0; i < n; i++)
    total_weight[i] = sample_count[i] = 0.;

  foreach (reduction(+:mean_grad[:n*len1]) reduction(+:mean_grad_sq[:n*len1])
           reduction(+:mean_strain_sq[:n*len2]) reduction(+:total_weight[:n])
           reduction(+:sample_count[:n])) {
    double hval = PROFILE_SLAB_COORD(point);
    if (hval < hmin - Delta/2. || hval > hmax + Delta/2.) continue;
    int iprof = (int)((hval - hmin)/deltahn + 0.5);
    if (iprof < 0) iprof = 0;
    if (iprof >= n) iprof = n - 1;

    double weight = (w.i != unity.i) ? w[] : 1.;
    total_weight[iprof] += weight;
    sample_count[iprof]++;

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

    double Sxx = dudx;
    double Sxy = 0.5*(dudy + dvdx);
    double Syy = dvdy;
    double S2  = sq (Sxx) + 2.*sq (Sxy) + sq (Syy);
    #if dimension == 3
      double Szz = dwdz;
      double Sxz = 0.5*(dwdx + dudz);
      double Syz = 0.5*(dwdy + dvdz);
      S2 += sq (Szz) + 2.*sq (Sxz) + 2.*sq (Syz);
    #endif

    double * g  = &mean_grad[iprof*len1];
    double * gs = &mean_grad_sq[iprof*len1];
    g[0] += dudx*weight;  gs[0] += sq (dudx)*weight;
    g[1] += dvdx*weight;  gs[1] += sq (dvdx)*weight;
    g[2] += dudy*weight;  gs[2] += sq (dudy)*weight;
    g[3] += dvdy*weight;  gs[3] += sq (dvdy)*weight;
    #if dimension == 3
      g[4] += dwdx*weight;  gs[4] += sq (dwdx)*weight;
      g[5] += dwdy*weight;  gs[5] += sq (dwdy)*weight;
      g[6] += dudz*weight;  gs[6] += sq (dudz)*weight;
      g[7] += dvdz*weight;  gs[7] += sq (dvdz)*weight;
      g[8] += dwdz*weight;  gs[8] += sq (dwdz)*weight;
    #endif

    double nu = 1.0;
    #ifdef mu
      #if defined(FOUR_PHASE)
        nu = mu (f1[], f2[], f3[])/rhov[];
      #elif defined(THREE_PHASE)
        nu = mu (f1[], f2[])/rhov[];
      #else
        nu = mu (f[])/rhov[];
      #endif
    #endif
    double * ss = &mean_strain_sq[iprof*len2];
    ss[0] += nu*S2*weight;
    ss[1] += sq (Sxx)*weight;
    ss[2] += sq (Sxy)*weight;
    ss[3] += sq (Syy)*weight;
    #if dimension == 3
      ss[4] += sq (Szz)*weight;
      ss[5] += sq (Sxz)*weight;
      ss[6] += sq (Syz)*weight;
    #endif
  }

  if (pid() == 0) {
    FILE * fp = fopen (filename, mode);
    if (fp == NULL) { perror (filename); exit (1); }

    fprintf (fp, "# Profile Dissipation: t = %.10g, L0 = %g\n", t, L0);
    fprintf (fp, "# [0]iprof [1]y [2]delta");
    int k = 3;
    for (int g = 0; g < len1; g++) fprintf (fp, " [%d]grad", k++);
    for (int g = 0; g < len1; g++) fprintf (fp, " [%d]grad_sq", k++);
    for (int g = 0; g < len2; g++) fprintf (fp, " [%d]strain_sq", k++);
    fputc ('\n', fp);

    for (int iprof = 0; iprof < n; iprof++) {
      double hprof = hmin + iprof*deltahn;
      double tw = total_weight[iprof];
      double sc = sample_count[iprof];
      #if (dimension == 2)
        double deltah = sc > 0. ? tw/sc : 0.;
      #else
        double deltah = sc > 0. ? sqrt (tw/sc) : 0.;
      #endif
      fprintf (fp, "%-6d %15.8e %15.8e", iprof, hprof, deltah);
      for (int k = 0; k < len1; k++)
        fprintf (fp, " %24.15e", tw > 0. ? mean_grad[iprof*len1 + k]/tw : 0.);
      for (int k = 0; k < len1; k++)
        fprintf (fp, " %24.15e", tw > 0. ? mean_grad_sq[iprof*len1 + k]/tw : 0.);
      for (int k = 0; k < len2; k++)
        fprintf (fp, " %24.15e", tw > 0. ? mean_strain_sq[iprof*len2 + k]/tw : 0.);
      fputc ('\n', fp);
    }
    fputc ('\n', fp);
    fputc ('\n', fp);
    fflush (fp);
    fclose (fp);
  }
}

#undef PROFILE_SLAB_COORD
#undef PROFILE_SLAB_PARAMS

#endif
