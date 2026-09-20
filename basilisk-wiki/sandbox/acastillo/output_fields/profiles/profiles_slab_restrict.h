/**
# Slab profiles on a mixed-resolution grid

The profiles of [profiles.h](profiles.h) in one traversal instead of one per
slab: `n` slabs cost one reduction, not `n`.

Binning is `foreach_level (slablevel)`, not `foreach()`: one level is one cell
size, so a bin is a plain average and `Delta` is constant. `restriction()`
fills the non-leaf cells first, so nonlinear quantities must be formed at the
leaves before it -- restriction commutes with nothing nonlinear.

`slablevel` is the coarsest level present everywhere the profile spans, and has
no default: `foreach_level(l)` visits exactly depth `l`, and a slab the tree is
coarser than comes back empty, indistinguishable from a zero average. `depth()`
is the trap -- always valid, right only on a uniform grid.

Include after `navier-stokes/centered.h`: needs `u` and `unity`. Output is in
`profiles.h` format. */

#ifndef PROFILES_SLAB_RESTRICT_H
#define PROFILES_SLAB_RESTRICT_H

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
  int slablevel = -1, \
  const char * mode = "a"

#if dimension == 2
  #define PROFILE_SLAB_COORD(point) y
#else
  #define PROFILE_SLAB_COORD(point) z
#endif

/** Refuse an unset `slablevel` rather than write zeros. */

static bool profile_slab_level_set (int slablevel, const char * filename)
{
  if (slablevel >= 0)
    return true;
  if (pid() == 0)
    fprintf (stderr, "%s: slablevel not set (the coarsest level present "
             "everywhere the profile spans). Nothing written.\n", filename);
  return false;
}

/** Fill the non-leaf cells the binning reads, `w` alongside the data. */

static void profile_slab_restrict (scalar * list, scalar w)
{
  scalar * r = list_copy (list);
  if (w.i != unity.i) r = list_add (r, w);
  restriction (r);
  free (r);
}

/** Slabs that got no cell: `slablevel` is finer than the grid there. */

#define PROFILE_SLAB_EPILOGUE() do { \
  int _e = 0; \
  for (int i = 0; i < n; i++) if (sample_count[i] <= 0.) _e++; \
  if (_e && pid() == 0) \
    fprintf (stderr, "%s: %d/%d slabs empty at level %d; those rows are " \
             "zeros, not averages\n", filename, _e, n, slablevel); \
} while (0)

/**
## Scalar profiles -- mean and mean-square of each field in `list`
*/

void profile_scalar_slab (scalar * list = all, PROFILE_SLAB_PARAMS)
{
  if (!profile_slab_level_set (slablevel, filename)) return;
  int len = list_len (list);
  double deltahn = (hmax - hmin)/((double) n - 0.99999999);

  double aver[n*len], aver_sq[n*len], total_weight[n], sample_count[n];
  for (int i = 0; i < n*len; i++)
    aver[i] = aver_sq[i] = 0.;
  for (int i = 0; i < n; i++)
    total_weight[i] = sample_count[i] = 0.;

  profile_slab_restrict (list, w);
  foreach_level (slablevel,
           reduction(+:aver[:n*len]) reduction(+:aver_sq[:n*len])
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

  PROFILE_SLAB_EPILOGUE();

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
      double deltah = sc > 0. ? pow (tw/sc, 1./(dimension - 1)) : 0.;
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
  if (!profile_slab_level_set (slablevel, filename)) return;
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

  /** Formed at the leaves: restricting the operands and multiplying here
      would give the product of the means. */

  scalar * prods = NULL;
  for (int k = 0; k < len; k++) {
    scalar p = new scalar;
    prods = list_append (prods, p);
  }
  foreach() {
    scalar s1, s2, p;
    for (s1, s2, p in list1, list2, prods)
      p[] = s1[]*s2[];
  }
  profile_slab_restrict (prods, w);
  foreach_level (slablevel,
           reduction(+:aver[:n*len]) reduction(+:total_weight[:n])
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
    for (scalar p in prods)
      aver[iprof*len + k++] += p[]*weight;
  }

  PROFILE_SLAB_EPILOGUE();

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
      double deltah = sc > 0. ? pow (tw/sc, 1./(dimension - 1)) : 0.;
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

  delete (prods);
  free (prods);
}

/**
## Dissipation profiles -- strain, enstrophy, and gradients for a Taylor scale

Five fields formed at the leaves then restricted, so the average is
$\overline{S_{ij}S_{ij}}$, not $\bar S_{ij}\bar S_{ij}$. The nine
$\partial_i u_j$ stay local doubles.

* `S2`, `O2` -- $S_{ij}S_{ij}$ and $\Omega_{ij}\Omega_{ij}$; their sum is
  $\sum(\partial_i u_j)^2$ pointwise.
* `nuS2` -- $\varepsilon/2$. `nu` varies with phase, so it is evaluated at a
  leaf where `f[]` is a real volume fraction.
* `dxx`, `dzz` -- the writer gives mean and mean square of each, hence
  $\overline{(\partial_i u_i')^2}$ and a longitudinal Taylor scale
  $\lambda_i^2 = \overline{u_i'^2}/\overline{(\partial_i u_i')^2}$, its
  numerator from a `profile_scalar_slab` call on `u`. $x$ is homogeneous, $z$
  carries the mean gradient, so the pair brackets the anisotropy. Longitudinal,
  so not comparable to a literature $Re_\lambda$.

Gradients are differenced at the leaves: this is the resolved dissipation, not
that of the filtered field. */

void profile_dissipation_slab (vector u, PROFILE_SLAB_PARAMS)
{
  if (!profile_slab_level_set (slablevel, filename)) return;

  scalar S2[], O2[], nuS2[], dxx[], dzz[];
  scalar * list = {S2, O2, nuS2, dxx, dzz};

  foreach() {
    /** `JJ[i][j]` is $\partial_j u_i$, built as `lambda2()` does: the
        component is held in a `scalar`, the stencil rotates. */
    double JJ[dimension][dimension];
    scalar s = u.x; int i = 0;
    foreach_dimension() JJ[0][i++] = center_gradient (s);
    s = u.y; i = 0;
    foreach_dimension() JJ[1][i++] = center_gradient (s);
    #if dimension == 3
      s = u.z; i = 0;
      foreach_dimension() JJ[2][i++] = center_gradient (s);
    #endif

    double s2 = 0., o2 = 0.;
    for (int i = 0; i < dimension; i++) {
      for (int j = 0; j < dimension; j++) {
        s2 += sq (0.5*(JJ[i][j] + JJ[j][i]));
        o2 += sq (0.5*(JJ[i][j] - JJ[j][i]));
      }
    }

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

    S2[]   = s2;
    O2[]   = o2;
    nuS2[] = nu*s2;
    dxx[]  = JJ[0][0];
    dzz[]  = JJ[dimension - 1][dimension - 1];
  }

  profile_scalar_slab (list, w, filename, hmin, hmax, xmin, xmax, ymin, ymax,
                       n, m1, m2, slablevel, mode);
}

#undef PROFILE_SLAB_COORD
#undef PROFILE_SLAB_EPILOGUE
#undef PROFILE_SLAB_PARAMS

#endif
