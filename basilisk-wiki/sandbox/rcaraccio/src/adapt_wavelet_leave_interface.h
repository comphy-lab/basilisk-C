/**
# Wavelet adaptation that keeps the interface at the finest level

This is a copy of Oystein Lande's function
([adapt_wavelet_leave_interface.h](/sandbox/oystelan/adapt_wavelet_leave_interface.h))
with some small (but relevant) tweaks.

The function takes a second list, 'vol_frac'. It refines every cell of that
list which holds an interface up to 'maxlevel', whatever the wavelet error is.
Use 'padding' to keep some layers of neighbour cells at the same level.

We identify the interface cells with the F_ERR macro. This keeps the threshold
equal to the other parts of this sandbox.

We also restrict 'vol_frac'. Basilisk's 'adapt_wavelet' restricts 'slist' only.
See the note in the code below.
*/
#if TREE

#ifndef F_ERR
# define F_ERR 1.e-10
#endif

astats adapt_wavelet_leave_interface(scalar *slist,      // list of scalars
                                     scalar *vol_frac,   // the volume fraction scalar
                                     double *max,        // tolerance for each scalar
                                     int maxlevel,       // maximum level of refinement
                                     int minlevel = 1,   // minimum level of refinement (default 1)
                                     int padding = 0,    // number of neighbor cells to padd on each side of the interface being preserved
                                     scalar *list = all) // list of fields to update
{
  scalar * ilist = list;

  /**
  We restrict 'vol_frac' together with 'slist'. Basilisk's 'adapt_wavelet'
  restricts 'slist' only, and the callers give the volume fraction in
  'vol_frac', never in 'slist'.

  Two steps read 'vol_frac' on the coarse levels. The 'vf[]' test below runs
  on children that can be internal cells. 'fraction_refine' builds the
  children of a cell from the coarse value of that cell. With MPI,
  'mpi_boundary_refine' applies 'fraction_refine' to cells that this process
  keeps only as coarse ghosts.

  Caution: keep this restriction. If the coarse value is not up to date, and
  it is between 0 and 1, 'fraction_refine' gives fractional children in the
  body of the particle. These artificial interface cells make the interface
  calculations incorrect.

  'restriction' skips the fields that are already up to date, thus the added
  cost is small.
  */

  scalar * listr = list_copy (slist);
  for (scalar vf in vol_frac)
    listr = list_add (listr, vf);

  if (is_constant(cm)) {
    if (list == NULL || list == all)
      list = list_copy (all);
  }
  else {
    if (list == NULL || list == all) {
      list = list_copy ({cm, fm});
      for (scalar s in all)
        list = list_add (list, s);
    }
    listr = list_add (listr, cm);
  }

  /**
  'boundary' acts on 'list', which holds every field of 'all', thus it covers
  'vol_frac' too. Do not add a second 'boundary' call for 'listr'. The two
  functions test different flags: 'boundary' tests 's_centered' and
  'restriction' tests 's_restriction'. 'boundary' does nothing for a field
  that is clean for 's_centered', and it leaves the coarse levels as they
  are. Only 'restriction' refreshes them.
  */

  boundary (list);
  restriction (listr);
  free (listr);

  astats st = {0, 0};
  scalar * listc = NULL;
  for (scalar s in list)
    listc = list_add_depend (listc, s);

  // refinement
  if (minlevel < 1)
    minlevel = 1;
  tree->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);
  foreach_cell() {
    if (is_active (cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
        if (cell.flags & too_coarse) {
          cell.flags &= ~too_coarse;
          refine_cell (point, listc, refined, &tree->refined);
          st.nf++;
        }
        continue;
      }
      else { // !is_leaf (cell)
        if (cell.flags & refined) {
          // cell has already been refined, skip its children
          cell.flags &= ~too_coarse;
          continue;
        }
        // check whether the cell or any of its children is local
        bool local = is_local (cell);
        if (!local)
          foreach_child() if (is_local (cell)) {
            local = true;
            break;
          }
        if (local) {
          int i = 0;
          static const int just_fine = 1 << (user + 3);
          for (scalar s in slist) {
            double emax = max[i++], sc[(1 << dimension) * s.block];
            double *b = sc;
            foreach_child()
                foreach_blockf(s)
                    *b++ = s[];
            s.prolongation(point, s);
            b = sc;
            foreach_child()
                foreach_blockf(s) {
              double e = fabs(*b - s[]);
              if (e > emax && level < maxlevel) {
                cell.flags &= ~too_fine;
                cell.flags |= too_coarse;
              }
              else if ((e <= emax / 1.5 || level > maxlevel) &&
                       !(cell.flags & (too_coarse | just_fine))) {
                if (level >= minlevel)
                  cell.flags |= too_fine;
              }
              else if (!(cell.flags & too_coarse)) {
                cell.flags &= ~too_fine;
                cell.flags |= just_fine;
              }
              // arnbo: always set interface cells to the finest level
              for (scalar vf in vol_frac) {
                if (vf[] > F_ERR && vf[] < 1. - F_ERR && level < maxlevel) {
                  cell.flags |= too_coarse;
                  cell.flags &= ~too_fine;
                  cell.flags &= ~just_fine;
                  if (padding > 0)
                  {
                    foreach_neighbor(padding)
                    {
                      cell.flags |= too_coarse;
                      cell.flags &= ~too_fine;
                      cell.flags &= ~just_fine;
                    }
                  }
                }
              }
              s[] = *b++;
            }
          }
          foreach_child() {
            cell.flags &= ~just_fine;
            if (!is_leaf(cell)) {
              cell.flags &= ~too_coarse;
              if (level >= maxlevel)
                cell.flags |= too_fine;
            }
            else if (!is_active(cell))
              cell.flags &= ~too_coarse;
          }
        }
      }
    }
    else // inactive cell
      continue;
  }
  mpi_boundary_refine(listc);

  // coarsening
  // the loop below is only necessary to ensure symmetry of 2:1 constraint
  for (int l = depth(); l >= 0; l--) {
    foreach_cell() if (!is_boundary(cell)) {
      if (level == l) {
        if (!is_leaf(cell)) {
          if (cell.flags & refined)
            // cell was refined previously, unset the flag
            cell.flags &= ~(refined | too_fine);
          else if (cell.flags & too_fine) {
            if (is_local(cell) && coarsen_cell(point, listc))
              st.nc++;
            cell.flags &= ~too_fine; // do not coarsen parent
          }
        }
        if (cell.flags & too_fine)
          cell.flags &= ~too_fine;
        else if (level > 0 && (aparent(0).flags & too_fine))
          aparent(0).flags &= ~too_fine;
        continue;
      }
      else if (is_leaf(cell))
        continue;
    }
    mpi_boundary_coarsen(l, too_fine);
  }
  free(listc);

  mpi_all_reduce(st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce(st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update(list);

  if (list != ilist)
    free(list);

  return st;
}

#endif

