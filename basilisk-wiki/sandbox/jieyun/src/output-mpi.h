/**
# Interface output with MPI

This is a modified version of the output_facets function in
[fractions.h](http://basilisk.fr/src/fractions.h#interface-output).
We can output facets to one single file during an MPI run, by performing the
output processor by processor.

Caution: The serial output implementation slows down the simulation.

## Facets from one volume fraction field
*/

trace
void output_facets_mpi (scalar c, char *file = "facets_vof.dat",
  face vector s = {{-1}}) {
  char name[strlen(file) + 2];
  strcpy (name, file);

  if (!s.x.i) s.x.i = -1;

  FILE *fp;

  int sign_mpi[2];
  vector nn[];

  // compute the normal first to avoid the deadlock.
  foreach()
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      foreach_dimension()
        nn.x[] = n.x;
    }
  boundary ((scalar *) {nn});

  // wait for the signal of previous rank.
  #if _MPI
  MPI_Status status;
  if (pid() != 0)
    MPI_Recv (&sign_mpi, 2, MPI_INT, pid() - 1, 0, MPI_COMM_WORLD, &status);
  #endif
  if (pid() == 0)
    fclose (fopen (name, "w"));

  fp = fopen (name, "a");
  foreach()
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord n;
      foreach_dimension()
        n.x = nn.x[];

      double alpha = plane_alpha (c[], n);
#if dimension == 2
      coord segment[2];
      if (facets (n, alpha, segment) == 2)
        fprintf (fp, "%g %g\n%g %g\n\n",
          x + segment[0].x*Delta, y + segment[0].y*Delta,
          x + segment[1].x*Delta, y + segment[1].y*Delta);
#else // dimension == 3
      coord v[12];
      int m = facets (n, alpha, v, 1.);
      for (int i = 0; i < m; i++)
        fprintf (fp, "%g %g %g\n",
          x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
      if (m > 0)
	      fputc ('\n', fp);
#endif
    }

  fflush (fp);
  fclose (fp);

  // send the signal to the next rank.
  #if _MPI
  if (pid() + 1 < npe())
    MPI_Send (&sign_mpi, 2, MPI_INT, pid() + 1, 0, MPI_COMM_WORLD);
  #endif
}

#if dimension <= 2
/**
## Facets from two volume fraction fields

Output the interface segments of two volume fraction field for symmetric
difference error computation. Reference phase is located on the left side of
a segment.
*/
void output_polygon_mpi (scalar f1, scalar f2, char *file,
  double epsf = 1.e-6) {
  char name[strlen(file) + 2];
  strcpy (name, file);

  face vector s;
  s.x.i = -1;

  FILE *fp;

  int sign_mpi[2];
  vector n1[], n2[];

  foreach() {
    if (f1[] > epsf && f1[] < 1. - epsf) {
      coord n = facet_normal (point, f1, s);
      foreach_dimension()
        n1.x[] = n.x;
    }
    if (f2[] > epsf && f2[] < 1. - epsf) {
      coord n = facet_normal (point, f2, s);
      foreach_dimension()
        n2.x[] = n.x;
    }
  }
  boundary ((scalar *) {n1, n2});

  #if _MPI
  MPI_Status status;
  if (pid() != 0)
    MPI_Recv (&sign_mpi, 2, MPI_INT, pid() - 1, 0, MPI_COMM_WORLD, &status);
  #endif
  if (pid() == 0)
    fclose (fopen (name, "w"));

  fp = fopen(name, "a");
  foreach(serial) {
    double c1 = clamp(f1[], 0., 1.), c2 = clamp(f2[], 0., 1.), cmin, cmax;
    cmin = min(c1, c2);
    cmax = max(c1, c2);
    if (cmax > epsf && cmin < 1. - epsf) {
      if (cmin < epsf)
        fprintf (fp, "%d %g %g %g %.12e %.12e\n\n", 0, x, y, Delta, c1, c2);
      else {
        coord segment1[2], segment2[2];
        foreach_dimension() {
          segment1[0].x = 0.;
          segment1[1].x = 0.;
          segment2[0].x = 0.;
          segment2[1].x = 0.;
        }

        fprintf (fp, "%d %g %g %g %.12e %.12e\n", 1, x, y, Delta, c1, c2);
        if (c1 > epsf && c1 < 1. - epsf) {
          coord n = {n1.x[], n1.y[]};
          double alpha = plane_alpha (f1[], n);
          if (facets (n, alpha, segment1) == 2) {
            coord dx;
            for (int is = 0; is < 2; is++) {
              segment1[is].x += 0.5;
              segment1[is].y += 0.5;
            }

            dx.x = segment1[1].x - segment1[0].x;
            dx.y = segment1[1].y - segment1[0].y;
            if (n.x*dx.y - n.y*dx.x < 0.) {
              swap(double, segment1[0].x, segment1[1].x);
              swap(double, segment1[0].y, segment1[1].y);
            }
          }
        }

        if (c2 > epsf && c2 < 1. - epsf) {
          coord n = {n2.x[], n2.y[]};
          double alpha = plane_alpha (f2[], n);
          if (facets (n, alpha, segment2) == 2) {
            coord dx;
            for (int is = 0; is < 2; is++) {
              segment2[is].x += 0.5;
              segment2[is].y += 0.5;
            }

            dx.x = segment2[1].x - segment2[0].x;
            dx.y = segment2[1].y - segment2[0].y;
            if (n.x*dx.y - n.y*dx.x < 0.) {
              swap(double, segment2[0].x, segment2[1].x);
              swap(double, segment2[0].y, segment2[1].y);
            }
          }
        }

        fprintf (fp, "%.12e %.12e %.12e %.12e\n%.12e %.12e %.12e %.12e\n\n",
          segment1[0].x, segment1[0].y,
          segment1[1].x, segment1[1].y,
          segment2[0].x, segment2[0].y,
          segment2[1].x, segment2[1].y);
      }
    }
  }

  fflush (fp);
  fclose (fp);

  #if _MPI
  if (pid() + 1 < npe())
    MPI_Send (&sign_mpi, 2, MPI_INT, pid() + 1, 0, MPI_COMM_WORLD);
  #endif
}
#endif

/**
# Cell number statistic
*/
event profile_vof (t = end) {
  int nc = grid -> n, tnc = grid -> tn;
  printf ("\n# Total number of (leaf) cells. This process (PID:%d): %d.\
    All processes: %d\n", pid(), nc, tnc);
}
