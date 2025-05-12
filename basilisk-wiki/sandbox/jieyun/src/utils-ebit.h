/**
# Various utility functions for the EBIT method

## Default parameters and variables
*/

@define exist(val) (val > 0. && val < 1.)
@define existNeg(val) (val > -1. && val < 0.)
#define within(v, vmin, vmax) (v <= vmax && v >= vmin)
// linear interpolation
#define DR_L(a) (max(0., 1. - fabs(a)))

#define bilinear_ebit(a1, a2, a3, a4, x, y) ((1. - y)*(1. - x)*a1 + (1. - y)*x*a2 \
  + y*x*a3 + y*(1. - x)*a4)

#define machine_zero 1.e-16

#ifndef DIT
#define DIT 1
#endif

#if _MYOUTPUT
#define OUTPUTPATH "output_tmp/"
#else
#define OUTPUTPATH ""
#endif

#ifndef TEST
#define TEST "output/"
#endif

#ifndef U_TEST
#ifndef BILINEAR
#define BILINEAR
#endif
#endif

#ifndef FIT_TEST
#ifndef CIRCLE_FIT
#define CIRCLE_FIT
#endif
#endif

double tTime = 0.;
event init (i = 0) {
  tTime = 0.;
}

/** color_pha[]: corner color vertex */
vertex scalar color_pha[], color_pha_new[];

#if dimension == 3
double volume, volume0, volume_int;
#else
double area, area0, area_int;
#endif

/**
## Functions for embedded boundary
*/

// metric for embedded boundary
(const) scalar cm_ebit = unity;
(const) face vector fm_ebit = unityf;

#if EMBED
event metric(i = 0) {
  fm_ebit = fs;
  cm_ebit = cs;
}

/** In order to simplify the marker advection near the embedded boundary, 
  we set the value of the ghost cell, but these values are not used in
  the NS solver. */
void boundary_ebit_embed (vector ue) {
  boundary ((scalar *) {ue});
  foreach() {
    if (cm_ebit[] <= 0.) {
      foreach_dimension() {
        if (cm_ebit[1] > 0.) {
          ue.x[] = -ue.x[1];
          ue.y[] = ue.y[1];
        }

        if (cm_ebit[-1] > 0.) {
          ue.x[] = -ue.x[-1];
          ue.y[] = ue.y[-1];
        }
      }
    }
  }
  boundary ((scalar *) {ue});
}

void boundary_ebit_embed_scalar (scalar se) {
  foreach() {
    if (cm_ebit[] <= 0.) {
      if (cm_ebit[1] > 0.)
        se[] = se[1];

      if (cm_ebit[-1] > 0.)
        se[] = se[-1];
    }
  }
}
#endif

/**
## Special iterators for the 3D EBIT method
*/

#if dimension == 3
# define foreach_edge_b()				\
    foreach_vertex()				\
      foreach_dimension()

macro foreach_edge_i(char flags = 0, Reduce reductions = None) {
  foreach_vertex(flags, reductions) {
    coord _a = {x - X0, y - Y0, z - Z0};
    foreach_dimension()
      if (_a.x + 0.5*Delta < L0)
        {...}
  }
}

macro foreach_edge_x(char flags = 0, Reduce reductions = None) {
  foreach_vertex(flags, reductions) {
    coord _a = {x - X0, y - Y0, z - Z0};
    if (_a.x + 0.5*Delta < L0)
      {...}
  }
}

macro foreach_edge_y(char flags = 0, Reduce reductions = None) {
  foreach_vertex(flags, reductions) {
    coord _a = {x - X0, y - Y0, z - Z0};
    if (_a.y + 0.5*Delta < L0)
      {...}
  }
}

macro foreach_edge_z(char flags = 0, Reduce reductions = None) {
  foreach_vertex(flags, reductions) {
    coord _a = {x - X0, y - Y0, z - Z0};
    if (_a.y + 0.5*Delta < L0)
      {...}
  }
}
#endif


/** Table for finding out the marker pair connected by the interface
 first index: dictionary index
 _[][i] = j: marker j is connected to marker i
*/
static int dict_to_markers[9][4] = {
  {-1, -1, -1, -1},
  {-1, 3, -1, 1},
  {2, -1, 0, -1},
  {1, 0, -1, -1},
  {-1, 2, 1, -1},
  {-1, -1, 3, 2},
  {3, -1, -1, 0},
  {3, 2, 1, 0},
  {1, 0, 3, 2}
};

/** Table for finding out all connected markers within a face (cell in 2D) 
 first index: dictionary index
 second index: the indices of marker
*/
static int dict_to_edge[9][4] = {
  {-1, -1, -1, -1},
  {3, 1, -1, -1},
  {0, 2, -1, -1},
  {0, 1, -1, -1},
  {1, 2, -1, -1},
  {2, 3, -1, -1},
  {3, 0, -1, -1},
  {3, 0, 1, 2},
  {0, 1, 2, 3}
};

static int edge_shift[12][3] = {
  {0, 0, 0},   /* 0 */
  {0, 1, 0}, /* 1 */
  {0, 1, 1},  /* 2 */
  {0, 0, 1},  /* 3 */
  {0, 0, 0},   /* 4 */
  {0, 0, 1},  /* 5 */
  {1, 0, 1},  /* 6 */
  {1, 0, 0},  /* 7 */
  {0, 0, 0}, /* 8 */
  {1, 0, 0}, /* 9 */
  {1, 1, 0},   /* 10 */
  {0, 1, 0}   /* 11 */
};

/** Transformation from local edge index (on each face) to global edge index (on cube) */
static int edge_face_to_cube[6][4] = {
  {8, 4, 11, 5}, /* x-, y-z plane, 0-4 */
  {9, 7, 10, 6}, /* x+, y-z plane, 0-4 */
  {0, 8, 3, 9},   /* y-, z-x plane, 0-4 */
  {1, 11, 2, 10}, /* y+, z-x plane, 0-4 */
  {4, 0, 7, 1}, /* z-, x-y plane, 0-4 */
  {5, 3, 6, 2}  /* z+, x-y plane, 0-4 */
};

/** Table for finding the adjacent faces sharing the same edge.
  First index i1: face index;
  Second index j: local edge index on face i1
  Third index: 0-2: index shift for face vector, 3: face index (i2) of adjacent face
  4: local edge index on face i2

  Using this table, we can find out the closed polygon within the each cell.
*/
static int face_to_face[6][4][5] = {
  {{0, 0, 0, 2, 1}, {0, 0, 0, 4, 0}, {0, 1, 0, 3, 1}, {0, 0, 1, 5, 0}}, /* x- plane, 0-4 edges */
  {{0, 0, 0, 2, 3}, {0, 0, 0, 4, 2}, {0, 1, 0, 3, 3}, {0, 0, 1, 5, 2}}, /* x+ plane, 0-4 edges */
  {{0, 0, 0, 4, 1}, {0, 0, 0, 0, 0}, {0, 0, 1, 5, 1}, {1, 0, 0, 1, 0}}, /* y- plane, 0-4 edges */
  {{0, 0, 0, 4, 3}, {0, 0, 0, 0, 2}, {0, 0, 1, 5, 3}, {1, 0, 0, 1, 2}}, /* y+ plane, 0-4 edges */
  {{0, 0, 0, 0, 1}, {0, 0, 0, 2, 0}, {1, 0, 0, 1, 1}, {0, 1, 0, 3, 0}}, /* z- plane, 0-4 edges */
  {{0, 0, 0, 0, 3}, {0, 0, 0, 2, 2}, {1, 0, 0, 1, 3}, {0, 1, 0, 3, 2}}  /* z+ plane, 0-4 edges */
};

/** Used for vertex vector */
void boundary_edge (vector v_e, vector v_buf) {
  foreach_dimension()
    v_buf.x.dirty = false;

  foreach_edge()
    v_buf.x[] = v_e.x[];

  boundary ((scalar *){v_e});
  foreach_edge()
    if (v_buf.x[] != v_e.x[])
      v_e.x[] = v_buf.x[];
}

/**
## Helper functions for retrieving the connectivity
*/
double get_dict_ind (int p1, int p2, int p3, int p4, int pc) {
  int ii = 0, ind, issame, col_ver[4] = {p1, p2, p3, p4};
  double conf;
  conf = 0.; // Without interface

  ind = p1 + p2 + p3 + p4;
  issame = (p1 == p2);

  if (pc == -1) { // Connect the opposite edges, one interface
    conf = (issame) ? 2. : 1.;
  }
  else if (ind == 3 || ind == 1) { // Connect the consecutive edges, one interface
    ii = 0;
    for (int iv = 0; iv < 4; iv++) {
      if (col_ver[iv] != pc) {
        ii = iv;
        break;
      }
    }
    conf = 3. + ii;
  }
  else if (ind == 2) { // Connect the consecutive edges, two interfaces
    conf = (col_ver[0] == pc) ? 7. : 8.;
  }
  return conf;
}

/** Given the dictionary index of cell and specific marker,
  output the marker connected to the input marker.
  This function is used to find out the markers for circle fit.*/
int get_end (int conf, int is) {
  int ie, dict[6] = {4, 2, 1, 3, 5, 3};
  if (conf == 8)
    conf = (is > 1) ? 5 : 3;
  else if (conf == 7)
    conf = (is == 1 || is == 2) ? 4 : 6;

  ie = dict[conf - 1] - is;
  return ie;
}


double get_scalar (Point point, scalar sc, int is, int js, int ks) {
  return sc[is,js,ks];
}

double get_vertex_scalar (Point point, vertex scalar sc, int is, int js, int ks) {
  return sc[is,js,ks];
}

/** Get the coodinate of the markers within specific cell. */
int get_segments_cell (Point point, scalar dict, face vector st, \
  coord *xm, int *nm, face vector sn = {{-1}}) {
  // return the two markers within the cell
  // st: tangential position of marker
  // sn: normal position of marker, default 0
  int ns = 0;
  // we need this because the value of zerof inside the function is not eaxctly zero
  double coef = 1.;
  face vector _sn;
  if (sn.x.i < 0) {
    coef = 0.;
    _sn = st;
  }
  else
    _sn = sn;

  int conf = (int) (dict[]);
  nm[0] = nm[1] = 0;

  if (conf == 0)
    return 0;

  double xy_edge[4][2];
  int ind_markers[4] = {-1, -1, -1, -1};

  xy_edge[0][0] = _sn.x[]*coef;
  xy_edge[0][1] = st.x[];

  xy_edge[1][0] = st.y[];
  xy_edge[1][1] = _sn.y[]*coef;

  xy_edge[2][0] = _sn.x[1]*coef + 1.;
  xy_edge[2][1] = st.x[1];

  xy_edge[3][0] = st.y[0,1];
  xy_edge[3][1] = _sn.y[0,1]*coef + 1.;

  for (int iind = 0; iind < 4; iind++)
    ind_markers[iind] = dict_to_edge[conf][iind];

  for (int iv = 0; iv < 4; iv += 2) {
    if (ind_markers[iv] != -1) {
      double x1, y1, x2, y2;
      int np = 0;
      int ie1 = ind_markers[iv], ie2 = ind_markers[iv + 1];

      x1 = xy_edge[ie1][0];
      y1 = xy_edge[ie1][1];

      x2 = xy_edge[ie2][0];
      y2 = xy_edge[ie2][1];

      xm[2*ns + np].x = x1; xm[2*ns + np].y = y1;
      xm[2*ns + np + 1].x = x2; xm[2*ns + np + 1].y = y2;
      np += 2;

      nm[ns] = np;
      ns++;
    }
  }
  return ns;
}

int get_segments (Point point, scalar dict, face vector st, \
  coord xm[8], int nm[2], face vector sn = {{-1}}) {
  // return the three/four markers used for fitting
  // st: tangential position of marker
  // sn: normal position of marker, default 0
  int ns = 0;
  double coef = 1.;
  face vector _sn;
  if (sn.x.i < 0) {
    coef = 0.;
    _sn = st;
  }
  else
    _sn = sn;

  int conf = (int) (dict[]);
  nm[0] = nm[1] = 0;

  if (conf == 0)
    return 0;

  double xy_edge[4][2];
  int ind_markers[4] = {-1, -1, -1, -1};

  xy_edge[0][0] = _sn.x[]*coef;
  xy_edge[0][1] = st.x[];

  xy_edge[1][0] = st.y[];
  xy_edge[1][1] = _sn.y[]*coef;

  xy_edge[2][0] = _sn.x[1]*coef + 1.;
  xy_edge[2][1] = st.x[1];

  xy_edge[3][0] = st.y[0,1];
  xy_edge[3][1] = _sn.y[0,1]*coef + 1.;

  for (int iind = 0; iind < 4; iind++)
    ind_markers[iind] = dict_to_edge[conf][iind];

  for (int iv = 0; iv < 4; iv += 2) {
    if (ind_markers[iv] != -1) {
      double x1, y1, x2, y2, x3, y3;
      int np = 0;
      int ie1 = ind_markers[iv], ie2 = ind_markers[iv + 1];

      x1 = xy_edge[ie1][0];
      y1 = xy_edge[ie1][1];

      x2 = xy_edge[ie2][0];
      y2 = xy_edge[ie2][1];

      xm[4*ns + np].x = x1; xm[4*ns + np].y = y1;
      xm[4*ns + np + 1].x = x2; xm[4*ns + np + 1].y = y2;
      np += 2;

      // find out the four markers used for fitting
      int iee[2] = {ie1, ie2}, ipx, ipy, conn, ipps, ippe;
      for (int ipe = 0; ipe < 2; ipe++) {
        if (iee[ipe] == 0) {
          ipx = -1;
          ipy = 0;
        }
        else if (iee[ipe] == 1) {
          ipx = 0;
          ipy = -1;
        }
        else if (iee[ipe] == 2) {
          ipx = 1;
          ipy = 0;
        }
        else {
          ipx = 0;
          ipy = 1;
        }
        conn = (int) dict[ipx,ipy];
        if (conn == 0)
          continue;
        ipps = 2*(iee[ipe] % 2 + 1) - iee[ipe];
        ippe = get_end (conn, ipps);

        if (ippe == 0) {
          x3 = _sn.x[ipx,ipy]*coef + ipx;
          y3 = st.x[ipx,ipy] + ipy;
        }
        else if (ippe == 1) {
          x3 = st.y[ipx,ipy] + ipx;
          y3 = _sn.y[ipx,ipy]*coef + ipy;
        }
        else if (ippe == 2) {
          x3 = _sn.x[ipx + 1,ipy]*coef + ipx + 1.;
          y3 = st.x[ipx + 1,ipy] + ipy;
        }
        else {
          x3 = st.y[ipx,ipy + 1] + ipx;
          y3 = _sn.y[ipx,ipy + 1]*coef + ipy + 1.;
        }

        xm[4*ns + np].x = x3;
        xm[4*ns + np].y = y3;

        np++;
      }

      nm[ns] = np;
      ns++;
    }
  }
  return ns;
}

// for test
const coord iadv = {0., 1., 2.};

/**
## Functions for IO
*/

void output_color_vertex (char *file, vertex scalar cv, \
  scalar cv_center = {{-1}}, face vector cv_fcenter = {{-1}},
  bool output_interfacial = true) {
  // file: compulsory, output file name
  // scalar cv: compulsory, corner color vertex
  // cv_center: optional, color vertex on cell center, for 2D (may be 3D later)
  // cv_fcenter: optional, color vertex on face center, for 3D

  char name[strlen(file) + 2];
  strcpy (name, file);

  FILE *fp1;

  #if _MPI
  int sign_mpi = 1;
  MPI_Status status;
  if (pid() != 0)
    MPI_Recv (&sign_mpi, 1, MPI_INT, pid() - 1, 0, MPI_COMM_WORLD, &status);
  #endif
  if (pid() == 0)
    fclose (fopen (name, "w"));
  fp1 = fopen (name, "a");

  // output the color vertices for test
  foreach_vertex(serial) {
    bool is_output = true;
    if (output_interfacial) {
      foreach_dimension()
        is_output = is_output && (cv[] == cv[1]) && (cv[] == cv[-1]);

      is_output = !is_output;
    }
    if (is_output) {
    #if dimension == 3
      fprintf (fp1, "%g %g %g %g %d %d\n", x, y, z, cv[], pid(), 1);
    #else
      fprintf (fp1, "%g %g %g %d %d\n", x, y, cv[], pid(), 1);
    #endif
    }
  }

  // output the central color point for test
#if dimension == 3
  foreach_face(serial) {
    fprintf (fp1, "%g %g %g %g %d %d\n", x, y, z, cv_fcenter.x[], pid(), -1);
  }
#else
  foreach(serial) {
    bool is_output = true;
    if (output_interfacial) {
      is_output = (cv_center[] != cv[]) || (cv_center[] != cv[1])\
        || (cv_center[] != cv[0,1]) || (cv_center[] != cv[1,1]);
    }
    if (is_output)
      fprintf (fp1, "%g %g %g %d %d\n", x, y, cv_center[], pid(), -1);
  }
#endif

  fflush (fp1);
  fclose (fp1);
  #if _MPI
  if (pid() + 1 < npe())
    MPI_Send (&sign_mpi, 1, MPI_INT, pid() + 1, 0, MPI_COMM_WORLD);
  #endif
}

void output_mesh (char *file) {
  char name[strlen(file) + 2];
  strcpy (name, file);

  FILE *fp1;

  #if _MPI
  int sign_mpi = 1;
  MPI_Status status;
  if (pid() != 0)
    MPI_Recv (&sign_mpi, 1, MPI_INT, pid() - 1, 0, MPI_COMM_WORLD, &status);
  #endif
  if (pid() == 0)
    fclose (fopen (name, "w"));
  fp1 = fopen (name, "a");

  // output the mesh for test
  foreach(serial) {
    double x1, y1, z1;
  #if dimension == 3
    for (int k = -1; k < 2; k += 2)
  #endif
    for (int j = -1; j < 2; j += 2)
      for (int i = -1; i < 2; i += 2) {
        x1 = x + 0.5*i*Delta;
        y1 = y + 0.5*j*Delta;
      #if dimension == 3
        z1 = z + 0.5*k*Delta;
        fprintf (fp1, "%g %g %g\n", x1, y1, z1);
      #else
        fprintf (fp1, "%g %g\n", x1, y1);
      #endif
      }
  }

  fflush (fp1);
  fclose (fp1);
  #if _MPI
  if (pid() + 1 < npe())
    MPI_Send (&sign_mpi, 1, MPI_INT, pid() + 1, 0, MPI_COMM_WORLD);
  #endif
}

/** Output the interface segments of two volume fraction field for symmetric
 difference error computation. Reference phase is located on the left side of
 a segment.
*/
void output_polygon_ebit_mpi (scalar f1, scalar f2, face vector s1, face vector s2, \
  scalar conf1, scalar conf2, vertex scalar col1, vertex scalar col2, \
  char *file, double epsf = 1.e-6)
{
  char name[strlen(file) + 2];
  strcpy (name, file);

  FILE *fp;

  int sign_mpi[2];

  face vector st;

  #if _MPI
  MPI_Status status;
  if (pid() != 0)
    MPI_Recv (&sign_mpi, 2, MPI_INT, pid() - 1, 0, MPI_COMM_WORLD, &status);
  #endif
  if (pid() == 0) {
    fclose (fopen (name, "w"));
  }

  fp = fopen (name, "a");
  foreach(serial) {
    double c1 = clamp(f1[], 0., 1.), c2 = clamp(f2[], 0., 1.), cmin, cmax;
    int confs[2] = {(int) (conf1[]), (int) (conf2[])};
    cmin = min(c1, c2);
    cmax = max(c1, c2);
    if (cmax > epsf && cmin < 1. - epsf) {
      if (cmin < epsf)
        fprintf (fp, "%d %g %g %g %.12e %.12e\n\n", 0, x, y, Delta, c1, c2);
      else {
        coord segs[4];
        for (int ii = 0; ii < 4; ii++) {
          foreach_dimension()
            segs[ii].x = 0.;
        }

        fprintf (fp, "%d %g %g %g %.12e %.12e\n", 1, x, y, Delta, c1, c2);

        for (int ii = 0; ii < 2; ii++) {
          double col;
          if (ii == 0) {
            st = s1;
            col = col1[];
          }
          else {
            st = s2;
            col = col2[];
          }

          coord xy_edge[4];
          int ind_markers[4] = {-1, -1, -1, -1};

          xy_edge[0].x = 0.;
          xy_edge[0].y = st.x[];

          xy_edge[1].x = st.y[];
          xy_edge[1].y = 0.;

          xy_edge[2].x = 1.;
          xy_edge[2].y = st.x[1];

          xy_edge[3].x = st.y[0,1];
          xy_edge[3].y = 1.;

          int conf = confs[ii];
          for (int iind = 0; iind < 4; iind++)
            ind_markers[iind] = dict_to_edge[conf][iind];

          for (int iv = 0; iv < 4; iv += 2) {
            if (ind_markers[iv] != -1) {
              coord p1, p2, dx, dm;
              int ie1 = ind_markers[iv], ie2 = ind_markers[iv + 1];

              foreach_dimension() {
                p1.x = xy_edge[ie1].x;
                p2.x = xy_edge[ie2].x;
                dx.x = p2.x - p1.x;
                dm.x = 0.5*(p2.x + p1.x);
              }

              double ds = sqrt(sq(dx.x) + sq(dx.y)) + 1.e-32;
              dx.x /= ds;
              dx.y /= ds;
              if ((col - 0.5)*(dm.x*dx.y - dm.y*dx.x) < 0.) {
                swap(double, p1.x, p2.x);
                swap(double, p1.y, p2.y);
              }

              foreach_dimension() {
                segs[2*ii].x = p1.x;
                segs[2*ii + 1].x = p2.x;
              }
            }
          }
        }

        fprintf (fp, "%.12e %.12e %.12e %.12e\n%.12e %.12e %.12e %.12e\n\n",
          segs[0].x, segs[0].y,
          segs[1].x, segs[1].y,
          segs[2].x, segs[2].y,
          segs[3].x, segs[3].y);
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

/**
## Functions for debugging
*/

#if _MYOUTPUT
FILE *fp_debug;

event defaults (i = 0) {
  color_pha.nodump = color_pha_new.nodump = true;
  char name[80];
  sprintf (name, "output_tmp/debug_log_ebit.dat");
  fp_debug = fopen (name, "w");
}

void debug_log (int i) {
  if (pid() == 0) {
  #if dimension == 3
    fprintf (fp_debug, "%d %g %.12e, %.12e\n", i, dt, volume, volume_int);
  #else
    fprintf (fp_debug, "%d %g %.12e, %.12e\n", i, dt, area, area_int);
  #endif
    fflush (fp_debug);
  }

}

event debug_end (t = end) {
  fclose (fp_debug);
}
#endif

event profile_ebit (t = end) {
  int nc = grid -> n, tnc = grid -> tn;
  printf ("\n# Total number of (leaf) cells. This process (PID:%d): %d.\
    All processes: %d\n", pid(), nc, tnc);
}
