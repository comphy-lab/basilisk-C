/**
# Main macros
*/

#ifndef SURFACTANT_TRANSPORT_H
#define SURFACTANT_TRANSPORT_H

/* Suggested by Z.Xue, option II is set as defualt */
#ifndef OPTION_II
#  define OPTION_II
#endif

/* Test, do not change */
#ifndef BYAOLD
#  define BYAOLD
#endif

#ifndef DEBUG
#  define DEBUG 0
#endif

/* Same value as F_ERR */
#ifndef S_ERR
#  define S_ERR 0.
#endif

/**
# Computation of the direction-split stretching term.
*/
double stretching_velocity_split(Point point, double A, coord n, face vector uf, int dim)
{
  if (A <= 0.) return 0.;

  normalize(&n);

  double dudx = (uf.x[1] / (fm.x[1] + SEPS) - uf.x[] / (fm.x[] + SEPS)) / Delta;
  double dudy = (0.5 * (uf.x[0, 1] / (fm.x[0, 1] + SEPS) + uf.x[1, 1] / (fm.x[1, 1] + SEPS)) -
                  0.5 * (uf.x[0, -1] / (fm.x[0, -1] + SEPS) + uf.x[1, -1] / (fm.x[1, -1] + SEPS))) /
    (2. * Delta);
  double dvdx = (0.5 * (uf.y[1, 0] / (fm.y[1, 0] + SEPS) + uf.y[1, 1] / (fm.y[1, 1] + SEPS)) -
                  0.5 * (uf.y[-1, 0] / (fm.y[-1, 0] + SEPS) + uf.y[-1, 1] / (fm.y[-1, 1] + SEPS))) /
    (2. * Delta);

  double dvdy = (uf.y[0, 1] / (fm.y[0, 1] + SEPS) - uf.y[] / (fm.y[] + SEPS)) / Delta;

  double S_split = 0.;
  if (dim == 0)
    S_split = (1. - sq(n.x)) * dudx - (n.x * n.y) * dvdx;
  else {
    S_split = (1. - sq(n.y)) * dvdy - (n.x * n.y) * dudy;
#if AXI
    S_split += (y > 0.) ? 0.5 * (uf.y[] / (fm.y[] + SEPS) + uf.y[0, 1] / (fm.y[0, 1] + SEPS)) / y : dvdy;
#endif
  }
  return A * S_split;
}

/**
## Function `facets()` with a less strict criterion

The original `facets()` function provided by Basilisk uses the criterion

`fabs(n.y) > 1.e-4`

which may cause problems for the species transport method of Xue et al.
Here, a less strict criterion is used.
*/

#if dimension <= 2
int facets_surfactant(coord n, double alpha, coord p[2])
{
  int i = 0;
  for (double s = -0.5; s <= 0.5; s += 1.)
    foreach_dimension()
    {
      if (fabs(n.y) > 1.e-30 && i < 2) {
        double a = (alpha - s * n.x) / n.y;
        if (a >= -0.5 && a <= 0.5) {
          p[i].x = s;
          p[i++].y = a;
        }
      }
    }
  return i;
}

double plane_area_center_surfactant(coord m, double alpha, coord *p)
{
  alpha += (m.x + m.y) / 2.;

  coord n = m;
  foreach_dimension(2) if (n.x < 0.)
  {
    alpha -= n.x;
    n.x = -n.x;
  }

  p->x = p->y = p->z = 0.;

  if (alpha <= 0. || alpha >= n.x + n.y) return 0.;

  foreach_dimension(2) if (n.x < 1.e-30)
  {
    p->x = 0.;
    p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;
    return 1.;
  }

  if (alpha >= n.x) {
    p->x += 1.;
    p->y += (alpha - n.x) / n.y;
  }
  else
    p->x += alpha / n.x;

  double ax = p->x, ay = p->y;
  if (alpha >= n.y) {
    p->y += 1.;
    ay -= 1.;
    p->x += (alpha - n.y) / n.x;
    ax -= (alpha - n.y) / n.x;
  }
  else {
    p->y += alpha / n.y;
    ay -= alpha / n.y;
  }

  foreach_dimension(2)
  {
    p->x /= 2.;
    p->x = clamp(p->x, 0., 1.);
    if (m.x < 0.) p->x = 1. - p->x;
    p->x -= 0.5;
  }

  return sqrt(ax * ax + ay * ay);
}
#endif

/**
## Data structure for interfacial fluxes

This structure stores the geometric information and the left/right fluxes
associated with one interfacial segment during the surfactant transport step.
*/

typedef struct {
  double x, y;
  double Aold, Astar;
  double Sdt;
  double Delta;
  int pid;
  double fluxl;
  double fluxr;
  double *pfluxl;
  double *pfluxr;
} FluxItem;

typedef struct {
  int n;
  scalar Aold, Astar;
  scalar segment;
  scalar Sdt;
  vector nold, nstar;
  vector p0, p1;
  vector pctr;
  vector *sgrads;
  face vector phi;
  face vector *psis;
} SurfTransInfo;

static inline void surftrans_init(SurfTransInfo *st, scalar *surfactants)
{
  st->n = 0;
  scalar Aold = new scalar;
  scalar Astar = new scalar;
  scalar segment = new scalar;
  scalar Sdt = new scalar;
  vector nold = new vector;
  vector nstar = new vector;
  vector p0 = new vector;
  vector p1 = new vector;
  vector pctr = new vector;
  face vector phi = new face vector;

  Aold.nodump = true;
  Astar.nodump = true;
  segment.nodump = true;
  Sdt.nodump = true;
  foreach_dimension()
  {
    nold.x.nodump = true;
    nstar.x.nodump = true;
    p0.x.nodump = true;
    p1.x.nodump = true;
    pctr.x.nodump = true;
    phi.x.nodump = true;
  }

  st->Aold = Aold;
  st->Astar = Astar;
  st->segment = segment;
  st->Sdt = Sdt;
  st->nold = nold;
  st->nstar = nstar;
  st->p0 = p0;
  st->p1 = p1;
  st->pctr = pctr;
  st->phi = phi;

  st->sgrads = NULL;
  st->psis = NULL;

  for (scalar surfactant in surfactants) {
    vector sg = new vector;
    foreach_dimension() sg.x.nodump = true;
    st->sgrads = vectors_append(st->sgrads, sg);

    face vector psi = new face vector;
    foreach_dimension() psi.x.nodump = true;
    st->psis = vectors_append(st->psis, psi);
  }
}

static inline void surftrans_destroy(SurfTransInfo *st)
{
  delete ((scalar *){st->Aold, st->Astar, st->segment, st->Sdt});
  delete ((scalar *){
    st->nold.x, st->nold.y, st->nstar.x, st->nstar.y, st->p0.x, st->p0.y, st->p1.x, st->p1.y, st->pctr.x, st->pctr.y});

  delete ((scalar *){st->phi.x, st->phi.y});

  if (st->sgrads) {
    for (vector sg in st->sgrads)
      delete ((scalar *){sg});
    free(st->sgrads);
    st->sgrads = NULL;
  }

  if (st->psis) {
    for (face vector psi in st->psis)
      delete ((scalar *){psi});
    free(st->psis);
    st->psis = NULL;
  }
}

void surfactant_clean_SurfTransInfo(Point point, SurfTransInfo *st)
{
  st->segment[] = 0.;
  st->Aold[] = 0.;
  st->Astar[] = 0.;
  st->Sdt[] = 0.;

  foreach_dimension()
  {
    st->nold.x[] = 0.;
    st->nstar.x[] = 0.;
    st->p0.x[] = 0.;
    st->p1.x[] = 0.;
    st->pctr.x[] = 0.;
  }
}

/**
## Store geometric information for surfactant transport

This function stores the geometric information associated with the old and
advected interfaces for the surfactant transport step.

For each interfacial cell, it computes and stores:

* the interfacial segment length/area,
* the interface normal,
* the interface center,
* the two segment endpoints,
* the direction-split stretching contribution.

The information is stored in the `SurfTransInfo` structure `st`.
*/

void surfactant_store_info(scalar cold, scalar cstar, face vector uf, int dim, SurfTransInfo *st)
{
  int localn = 0;
  foreach () {
    surfactant_clean_SurfTransInfo(point, st);

    if ((cold[] > S_ERR && cold[] < 1. - S_ERR) || (cstar[] > S_ERR && cstar[] < 1. - S_ERR)) localn++;

    if (cold[] > S_ERR && cold[] < 1. - S_ERR) {
      coord nold = interface_normal(point, cold), pold;
      double alpha_old = plane_alpha(cold[], nold);
      double seg = plane_area_center_surfactant(nold, alpha_old, &pold) * Delta;
      st->segment[] = seg;
      st->Aold[] = seg;
#if AXI
      st->Aold[] *= y + pold.y * Delta;
#endif

      foreach_dimension() { st->nold.x[] = nold.x; }

      st->pctr.x[] = x + pold.x * Delta;
      st->pctr.y[] = y + pold.y * Delta;

      coord p[2];
      p[0].x = p[0].y = p[1].x = p[1].y = 0.;
      if (facets_surfactant(nold, alpha_old, p) == 2) {
        foreach_dimension()
        {
          st->p0.x[] = p[0].x;
          st->p1.x[] = p[1].x;
        }
      }
      else
        fprintf(stderr, "WARNING : facets = %d, p[0].x = %.12e, p[0].y = %.12e, p[1].x = %.12e, p[1].y = %.12e\n",
          facets_surfactant(nold, alpha_old, p), p[0].x, p[0].y, p[1].x, p[1].y);

      st->Sdt[] = stretching_velocity_split(point, st->Aold[], nold, uf, dim) * dt;
    }

    if (cstar[] > S_ERR && cstar[] < 1. - S_ERR) {
      coord nstar = interface_normal(point, cstar), pstar;
      double alpha_star = plane_alpha(cstar[], nstar);

      st->Astar[] = plane_area_center_surfactant(nstar, alpha_star, &pstar) * Delta;
#if AXI
      st->Astar[] *= y + pstar.y * Delta;
#endif
      foreach_dimension() st->nstar.x[] = nstar.x;
    }
  }
  st->n = localn;
  boundary((scalar *){st->segment, st->Aold, st->Astar, st->Sdt});
  boundary((scalar *){
    st->nold.x, st->nold.y, st->nstar.x, st->nstar.y, st->p0.x, st->p0.y, st->p1.x, st->p1.y, st->pctr.x, st->pctr.y});
}

/**
## Comparison function for sorting flux items

This function is used to sort `FluxItem` entries first according to their
transverse coordinate and then according to their streamwise coordinate.

A tolerance proportional to the grid size is used when comparing the
transverse coordinate.
*/

foreach_dimension() 
int compare_x(const void *a, const void *b)
{
  const FluxItem *sa = a;
  const FluxItem *fb = b;
  if (fabs(sa->y - fb->y) > 1.e-3 * L0 / (1 << grid->maxdepth)) return (sa->y < fb->y) ? -1 : 1;
  return (sa->x < fb->x) ? -1 : 1;
}

/**
## Compute a consistent interfacial area flux

This function computes a consistent flux of interfacial area for the
direction-split surfactant transport step.

It first stores the geometric information associated with the old interface
`cold` and the advected interface `cstar`, then builds a list of
interfacial cells and computes the left/right area fluxes along the current
sweep direction.

The fluxes are determined so that the variation of interfacial area is
consistent with the stretching contribution stored in `st->Sdt`.

In MPI mode, the interfacial-cell information is gathered from all
processes, the fluxes are computed globally, and the resulting values are
scattered back to the local process. The final fluxes are stored in
`st->phi`.
*/
void compute_area_flux_consistent(scalar cold, scalar cstar, face vector uf, int dim, SurfTransInfo *st)
{
  boundary((scalar *){cold, cstar, uf});
  surfactant_store_info(cold, cstar, uf, dim, st);

  FluxItem *items = NULL;
  if (st->n > 0) {
    items = (FluxItem *)malloc(sizeof(FluxItem) * st->n);
  }

  foreach_face() st->phi.x[] = 0.;

  int idx = 0;
  foreach (serial) {
    if ((cold[] > S_ERR && cold[] < 1. - S_ERR) || (cstar[] > S_ERR && cstar[] < 1. - S_ERR)) {
      if (idx < st->n) {
        items[idx].x = x;
        items[idx].y = y;
        items[idx].Aold = st->Aold[];
        items[idx].Astar = st->Astar[];
        items[idx].Sdt = st->Sdt[];
        items[idx].Delta = Delta;
        items[idx].fluxl = 0.;
        items[idx].fluxr = 0.;

#if _MPI
        items[idx].pid = pid();
#else
        if (dim == 0) {
          items[idx].pfluxl = &st->phi.x[];
          items[idx].pfluxr = &st->phi.x[1];
        }
        else {
          items[idx].pfluxl = &st->phi.y[];
          items[idx].pfluxr = &st->phi.y[0, 1];
        }
#endif
        idx++;
      }
    }
  }

#if _MPI
  int n_procs = npe();
  int *counts = (int *)malloc(n_procs * sizeof(int));
  int *displs = (int *)malloc(n_procs * sizeof(int));

  int localn = (int)st->n;
  MPI_Allgather(&localn, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);

  int cellTot = 0;
  for (int i = 0; i < n_procs; i++) {
    displs[i] = cellTot * sizeof(FluxItem);
    cellTot += counts[i];
    counts[i] *= sizeof(FluxItem);
  }

  FluxItem *gitems = (FluxItem *)malloc(cellTot * sizeof(FluxItem));
  MPI_Allgatherv(items, localn * sizeof(FluxItem), MPI_BYTE, gitems, counts, displs, MPI_BYTE, MPI_COMM_WORLD);

  if (items) free(items);
  free(counts);
  free(displs);

  FluxItem *solve_list = gitems;
  int solve_count = cellTot;

#else
  FluxItem *solve_list = items;
  int solve_count = st->n;
#endif

  if (solve_count > 0) {
    if (dim == 0)
      qsort(solve_list, solve_count, sizeof(FluxItem), compare_x);
    else
      qsort(solve_list, solve_count, sizeof(FluxItem), compare_y);

    int i = 0;
    while (i < solve_count) {
      int start = i;
      int end = i;

      while (end + 1 < solve_count) {
        FluxItem *curr = &solve_list[end];
        FluxItem *next = &solve_list[end + 1];
        bool same_line = false, neighbors = false;

        if (dim == 0) {
          same_line = (fabs(curr->y - next->y) < 1.e-3 * L0 / (1 << grid->maxdepth));
          double dist = next->x - curr->x;
          neighbors = (dist > 0 && dist < 1.6 * curr->Delta);
        }
        else {
          same_line = (fabs(curr->x - next->x) < 1.e-3 * L0 / (1 << grid->maxdepth));
          double dist = next->y - curr->y;
          neighbors = (dist > 0 && dist < 1.6 * curr->Delta);
        }

        if (same_line && neighbors)
          end++;
        else
          break;
      }

      double dATot = 0., SdtTot = 0., WTot = 0.;
      for (int k = start; k <= end; k++) {
        dATot += (solve_list[k].Astar - solve_list[k].Aold);
        SdtTot += solve_list[k].Sdt;
#ifndef BYAOLD
        WTot += fabs(solve_list[k].Astar - solve_list[k].Aold);
#else
        WTot += solve_list[k].Aold;
#endif
      }

      double total_error = dATot - SdtTot;
      double current_flux = 0.0;

      for (int k = start; k <= end; k++) {
        solve_list[k].fluxl = current_flux;
        double local_error = 0.;
#ifndef BYAOLD
        local_error = total_error * (fabs(solve_list[k].Astar - solve_list[k].Aold) / WTot);
#else
        local_error = total_error * (solve_list[k].Aold / WTot);
#endif
        double flux_out = current_flux - (solve_list[k].Astar - solve_list[k].Aold) + solve_list[k].Sdt + local_error;
        solve_list[k].fluxr = flux_out;
        current_flux = flux_out;
      }
      i = end + 1;
    }
  }

#if _MPI
  for (int k = 0; k < solve_count; k++) {
    if (solve_list[k].pid == pid()) {
      Point point = locate(solve_list[k].x, solve_list[k].y);
      if (point.level >= 0) {
        if (dim == 0) {
          st->phi.x[] = solve_list[k].fluxl;
          st->phi.x[1] = solve_list[k].fluxr;
        }
        else {
          st->phi.y[] = solve_list[k].fluxl;
          st->phi.y[0, 1] = solve_list[k].fluxr;
        }
      }
    }
  }
  if (gitems) free(gitems);
#else
  for (int k = 0; k < solve_count; k++) {
    *(solve_list[k].pfluxl) = solve_list[k].fluxl;
    *(solve_list[k].pfluxr) = solve_list[k].fluxr;
  }
  if (items) free(items);
#endif
  boundary((scalar *){st->phi});
}

/**
## Reconstruct the interfacial surfactant gradient

This function reconstructs the tangential gradient of each surfactant field
on the interface.

For each interfacial cell, a local weighted least-squares reconstruction is
performed using neighboring interfacial cells. The full gradient is first
estimated from the neighboring interface-center values, and its normal
component is then removed so that only the tangential component is retained.

The reconstructed gradients are stored in `st->sgrads`.
*/
//#if 0
void reconstruct_surfactant_gradient(scalar *surfactants, SurfTransInfo *st)
{
  int ii = 0;
  for (scalar surfactant in surfactants) {
    vector sg = st->sgrads[ii++];

    foreach () {
      foreach_dimension() sg.x[] = 0.;

      if (st->segment[] > 0.) {
        coord nctr = {st->nold.x[], st->nold.y[]};
        double abs_pctr_x = st->pctr.x[];
        double abs_pctr_y = st->pctr.y[];
        double ctr_val = surfactant[];

        double Mxx = 0., Mxy = 0., Myy = 0.;
        double Rx = 0., Ry = 0.;
        double epsilon = 1.e-20;

        foreach_neighbor(1)
        {
          if (point.i == 0 && point.j == 0) continue;

          if (st->segment[] > 1.e-6 * L0 / (1 << grid->maxdepth)) {
            double abs_pnbr_x = st->pctr.x[];
            double abs_pnbr_y = st->pctr.y[];
#if AXI
            double w = st->Aold[] / Delta;
#else
            double w = st->segment[] / Delta;
#endif
            double dx = abs_pnbr_x - abs_pctr_x;
            double dy = abs_pnbr_y - abs_pctr_y;
            double dS = surfactant[] - ctr_val;

            Mxx += w * dx * dx;
            Mxy += w * dx * dy;
            Myy += w * dy * dy;
            Rx += w * dx * dS;
            Ry += w * dy * dS;
          }
        }

        Mxx += epsilon;
        Myy += epsilon;

        double det = Mxx * Myy - Mxy * Mxy;
        double Gx_full = 0., Gy_full = 0.;

        if (fabs(det) > 1.e-30) {
          Gx_full = (Myy * Rx - Mxy * Ry) / det;
          Gy_full = (Mxx * Ry - Mxy * Rx) / det;
        }

        normalize(&nctr);
        double n_dot_grad = nctr.x * Gx_full + nctr.y * Gy_full;

        sg.x[] = Gx_full - nctr.x * n_dot_grad;
        sg.y[] = Gy_full - nctr.y * n_dot_grad;
      }
    }
    boundary((scalar *){sg});
  }
}
//#endif

/**
## Advect surfactant consistently with the interfacial area flux

This function advects the surfactant fields consistently with the
direction-split interfacial area flux stored in `st->phi`.

It first reconstructs the tangential surfactant gradients on the interface.
For each surfactant field, the corresponding surfactant flux is then
computed on faces from the upwind interfacial segment, taking into account
the reconstructed tangential distribution along the interface.

Finally, the surfactant concentration is updated in each interfacial cell
from the old surfactant mass and the net surfactant flux, using the new
interfacial area.
*/
void advect_surfactant_consistent(scalar *surfactants, face vector uf, int dim, SurfTransInfo *st)
{
  // return unnormalized face nmoral n1, surface area segment, surfactant gradient sgrad
  reconstruct_surfactant_gradient(surfactants, st);

  // double lengthMin = 1.e-25 * L0 / (1 << grid->maxdepth);
  double lengthMin = 0.;

  int ii = 0;
  for (scalar surfactant in surfactants) {
    vector sg = st->sgrads[ii];
    face vector psi = st->psis[ii];
    ii++;

    foreach_face() psi.x[] = 0.;

    if (dim == 0) {
      foreach_face(x)
      {
        if (fabs(st->phi.x[]) > lengthMin) {
          int up = (st->phi.x[] > 0.) ? -1 : 0;

          if (st->segment[up] > lengthMin) {
            coord n_v = {st->nold.x[up], st->nold.y[up]};
            coord t_v = {-n_v.y, n_v.x};
            normalize(&t_v);

            double p0_x = st->p0.x[up];
            double p0_y = st->p0.y[up];
            double p1_x = st->p1.x[up];
            double p1_y = st->p1.y[up];

            double seg_x = p1_x - p0_x;
            double seg_y = p1_y - p0_y;

            if (seg_x * t_v.x + seg_y * t_v.y < 0.) {
              double tmp_x = p0_x, tmp_y = p0_y;
              p0_x = p1_x;
              p0_y = p1_y;
              p1_x = tmp_x;
              p1_y = tmp_y;
            }
            double area_tot = 0.;
#if AXI
            double r0 = y + p0_y * Delta;
            double r1 = y + p1_y * Delta;
            double r_mid = 0.5 * (r0 + r1);
            area_tot = st->segment[up] * r_mid;
#else
            area_tot = st->segment[up];
#endif

            double G_s = sg.x[up] * t_v.x + sg.y[up] * t_v.y;
#if AXI
            double dr = r1 - r0;
            double S_scal = st->segment[up] / fabs(dr);
            double geom = (2.0 / 3.0) * (sq(r1) + r1 * r0 + sq(r0)) / (r1 + r0);
            double Gamma_c = surfactant[up] - G_s * S_scal * (geom - r0);
#endif

#if AXI
            if (st->phi.x[] > area_tot)
              st->phi.x[] = area_tot;
            else if (st->phi.x[] < -area_tot)
              st->phi.x[] = -area_tot;
#else
            if (st->phi.x[] > st->segment[up])
              st->phi.x[] = st->segment[up];
            else if (st->phi.x[] < -st->segment[up])
              st->phi.x[] = -st->segment[up];
#endif
            double flux_vol = fabs(st->phi.x[]);
            double exit_rel = (st->phi.x[] > 0.) ? 0.5 : -0.5;

            double dist0 = fabs(p0_x - exit_rel);
            double dist1 = fabs(p1_x - exit_rel);

            double fraction = area_tot > SEPS ? clamp(flux_vol / area_tot, 0., 1.) : 0.;
            double psi_val = 0.;
#if AXI
            double y_out = (dist0 < dist1) ? r0 : r1;
            double y_in = (dist0 < dist1) ? r1 : r0;
            double term = sq(y_out) - fraction * (sq(y_out) - sq(y_in));
            double y_phi = (term > 0.) ? sqrt(term) : 0.;

            double K_common = S_scal;
            double term_lin = G_s * S_scal;
            double K1 = K_common * (Gamma_c - term_lin * r0);
            double K2 = K_common * term_lin;

            double P_out = K1 * (sq(y_out) / 2.0) + K2 * (pow(y_out, 3) / 3.0);
            double P_phi = K1 * (sq(y_phi) / 2.0) + K2 * (pow(y_phi, 3) / 3.0);

            psi_val = fabs(P_out - P_phi);
#else
            double len_flux = area_tot * fraction;
            double s_start = 0., s_end = 0.;

            if (dist0 < dist1) {
              s_start = 0.;
              s_end = len_flux;
            }
            else {
              s_start = st->segment[up] - len_flux;
              s_end = st->segment[up];
            }

            double intercept = surfactant[up] - G_s * st->segment[up] * 0.5;
            double P_end = G_s * sq(s_end) / 2.0 + intercept * s_end;
            double P_start = G_s * sq(s_start) / 2.0 + intercept * s_start;
            psi_val = fabs(P_end - P_start);
#endif
            psi.x[] = (st->phi.x[] > 0.) ? psi_val : -psi_val;
          }
        }
      }
    }
    else if (dim == 1) { // dim == 1 (Y-Direction)
      foreach_face(y)
      {
        if (fabs(st->phi.y[]) > lengthMin) {
          int up = (st->phi.y[] > 0.) ? -1 : 0;
          if (st->segment[0, up] > lengthMin) {
            coord n_v = {st->nold.x[0, up], st->nold.y[0, up]};
            coord t_v = {-n_v.y, n_v.x};
            normalize(&t_v);

            double p0_x = st->p0.x[0, up];
            double p0_y = st->p0.y[0, up];
            double p1_x = st->p1.x[0, up];
            double p1_y = st->p1.y[0, up];

            double seg_x = p1_x - p0_x;
            double seg_y = p1_y - p0_y;

            if (seg_x * t_v.x + seg_y * t_v.y < 0.) {
              double tmp_x = p0_x, tmp_y = p0_y;
              p0_x = p1_x;
              p0_y = p1_y;
              p1_x = tmp_x;
              p1_y = tmp_y;
            }

            double area_tot = 0.;
#if AXI
            double y_cent = y + up * Delta / 2.;
            double r0 = y_cent + p0_y * Delta;
            double r1 = y_cent + p1_y * Delta;
            double dr = r1 - r0;
            double r_mid = 0.5 * (r0 + r1);
            area_tot = st->segment[0, up] * r_mid;
#else
            area_tot = st->segment[0, up];
#endif
            double G_s = sg.x[0, up] * t_v.x + sg.y[0, up] * t_v.y;
#if AXI
            double S_scal = st->segment[0, up] / fabs(dr);
            double geom = (2.0 / 3.0) * (sq(r1) + r1 * r0 + sq(r0)) / (r1 + r0);
            double Gamma_c = surfactant[0, up] - G_s * S_scal * (geom - r0);
#endif

#if AXI
            if (st->phi.y[] > area_tot)
              st->phi.y[] = area_tot;
            else if (st->phi.y[] < -area_tot)
              st->phi.y[] = -area_tot;
#else
            if (st->phi.y[] > st->segment[0, up])
              st->phi.y[] = st->segment[0, up];
            else if (st->phi.y[] < -st->segment[0, up])
              st->phi.y[] = -st->segment[0, up];
#endif

            double flux_vol = fabs(st->phi.y[]);
            double exit_rel = (st->phi.y[] > 0.) ? 0.5 : -0.5;

            double dist0 = fabs(p0_y - exit_rel);
            double dist1 = fabs(p1_y - exit_rel);

            double fraction = area_tot > SEPS ? clamp(flux_vol / area_tot, 0., 1.) : 0.;
            double psi_val = 0.;

#if AXI
            double y_out = (dist0 < dist1) ? r0 : r1;
            double y_in = (dist0 < dist1) ? r1 : r0;

            double term = sq(y_out) - fraction * (sq(y_out) - sq(y_in));
            double y_phi = (term > 0.) ? sqrt(term) : 0.;

            double K_common = S_scal;
            double term_lin = G_s * S_scal;
            double K1 = K_common * (Gamma_c - term_lin * r0);
            double K2 = K_common * term_lin;

            double P_out = K1 * (sq(y_out) / 2.0) + K2 * (pow(y_out, 3) / 3.0);
            double P_phi = K1 * (sq(y_phi) / 2.0) + K2 * (pow(y_phi, 3) / 3.0);

            psi_val = fabs(P_out - P_phi);
#else
            double len_flux = area_tot * fraction;
            double s_start = 0., s_end = 0.;

            if (dist0 < dist1) {
              s_start = 0.;
              s_end = len_flux;
            }
            else {
              s_start = st->segment[0, up] - len_flux;
              s_end = st->segment[0, up];
            }
            double intercept = surfactant[0, up] - G_s * st->segment[0, up] * 0.5;
            double P_end = G_s * sq(s_end) / 2.0 + intercept * s_end;
            double P_start = G_s * sq(s_start) / 2.0 + intercept * s_start;
            psi_val = fabs(P_end - P_start);
#endif
            psi.y[] = (st->phi.y[] > 0.) ? psi_val : -psi_val;
          }
        }
      }
    }
    boundary((scalar *){psi});

    foreach () {
      if (st->Astar[] > lengthMin) {
        double flux_net = dim == 0 ? psi.x[] - psi.x[1] : psi.y[] - psi.y[0, 1];

        double mass_old = surfactant[] * st->Aold[];
        double mass_new = mass_old + flux_net;
        if (mass_new < 0.) mass_new = 0.;

        double Astar = 0.;

#ifdef OPTION_II
        double area_net = dim == 0 ? st->phi.x[] - st->phi.x[1] : st->phi.y[] - st->phi.y[0, 1];
        Astar = st->Aold[] + area_net + st->Sdt[];
#else
        Astar = st->Astar[];
#endif
        surfactant[] = (Astar > lengthMin) ? mass_new / Astar : 0.;
      }
      else
        surfactant[] = 0.;
    }
    boundary({surfactant});
  }
}
#endif

/**
## References
~~~bib
@misc{xue2025,
  title={A sharp and conservative method for modeling interfacial flows with insoluble surfactants in the framework of a geometric volume-of-fluid approach},
  author={Zhong-Han Xue and Jacques Magnaudet and Jie Zhang},
  year={2025},
  eprint={2507.09680},
  archivePrefix={arXiv},
  primaryClass={physics.flu-dyn},
  url={https://arxiv.org/abs/2507.09680}
}
~~~
*/