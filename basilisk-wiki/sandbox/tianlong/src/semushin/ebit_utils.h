@define exist(val) (val > 0. && val < 1.)
@define existNeg(val) (val > -1. && val < 0.)

// linear interpolation
#define DR_L(a) (max(0., 1. - fabs(a)))
#define machine_zero 1.e-16

#ifndef DIT
#define DIT 1
#endif

#if dimension == 3
# define foreach_edge_b()				\
    foreach_vertex()				\
      foreach_dimension()

@def foreach_vertex_myaux()
foreach_vertex() {
  coord _a = {x, y, z};
@
@define end_foreach_vertex_myaux() } end_foreach_vertex()

# define foreach_edge_i() \
    foreach_vertex_myaux() \
      foreach_dimension() \
        if (_a.x + 0.5 * Delta < L0)
#endif

# define foreach_edge_x()				\
    foreach_vertex_myaux()				\
      if (_a.x + 0.5 * Delta < L0)
# define foreach_edge_y()				\
    foreach_vertex_myaux()				\
      if (_a.y + 0.5 * Delta < L0)
# define foreach_edge_z()				\
    foreach_vertex_myaux()				\
      if (_a.z + 0.5 * Delta < L0)

double get_dict_ind(int p1, int p2, int p3, int p4, int pc) {
  int ii=0, ind, issame, col_ver[4]={p1, p2, p3, p4};
  double conf;
  conf = 0.; // Without interface

  ind = p1 + p2 + p3 + p4;
  issame = (p1 == p2);

  if (pc == -1) { // Connect the opposite edges, one interface
    conf = (issame) ? 2. : 1.;
  }
  else if (ind == 3 || ind == 1) { // Connect the consecutive edges, one interface
    ii = 0;
    for (int iv=0; iv<4; iv++) {
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
int get_end(int conf, int is) {
  int ie, dict[6]={4, 2, 1, 3, 5, 3};
  if (conf == 8)
    conf = (is > 1) ? 5 : 3;
  else if (conf == 7)
    conf = (is == 1 || is == 2) ? 4 : 6;

  ie = dict[conf - 1] - is;
  return ie;
}


double get_scalar(Point point, scalar sc, int is, int js, int ks) {
  return sc[is, js, ks];
}

double get_vertex_scalar(Point point, vertex scalar sc, int is, int js, int ks) {
  return sc[is, js, ks];
}


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
  {1, 3, -1, -1},
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

// for test
coord iadv = {0., 1., 2.};

/** Functions for IO */
struct ColorVertex {
  char *file; // compulsory, output file name
  vertex scalar cv; // compulsory, corner color vertex
  scalar cv_center; // optional, color vertex on cell center, for 2D (may be 3D later)
  face vector cv_fcenter; // optional, color vertex on face center, for 3D
};

void output_color_vertex(struct ColorVertex a) {
  char *file = a.file;
  vertex scalar cv = a.cv;
  scalar cv_center = a.cv_center;
  face vector cv_fcenter = a.cv_fcenter;

  char name[strlen(file) + 2];
  strcpy (name, file);

  FILE *fp1;

  #if _MPI
  int sign_mpi=1;
  MPI_Status status;
  if (pid() != 0)
    MPI_Recv(&sign_mpi, 1, MPI_INT, pid() - 1, 0, MPI_COMM_WORLD, &status);
  #endif
  if (pid() == 0)
    fclose(fopen(name, "w"));
  fp1 =fopen(name, "a");

  // output the color vertices for test
  foreach_vertex(serial) {
  #if dimension == 3
    fprintf (fp1, "%g %g %g %g %d\n", x, y, z, cv[], pid());
  #else
    fprintf (fp1, "%g %g %g %d\n", x, y, cv[], pid());
  #endif
  }

  // output the central color point for test
#if dimension == 3
  foreach_face(serial) {
    fprintf (fp1, "%g %g %g %g %d\n", x, y, z, cv_fcenter.x[], pid());
  }
#else
  foreach(serial) {
    fprintf (fp1, "%g %g %g %d\n", x, y, cv_center[], pid());
  }
#endif

  fflush(fp1);
  fclose(fp1);
  #if _MPI
  if (pid() + 1 < npe())
    MPI_Send(&sign_mpi, 1, MPI_INT, pid() + 1, 0, MPI_COMM_WORLD);
  #endif
}

void output_mesh(char *file) {
  char name[strlen(file) + 2];
  strcpy (name, file);

  FILE *fp1;

  #if _MPI
  int sign_mpi=1;
  MPI_Status status;
  if (pid() != 0)
    MPI_Recv(&sign_mpi, 1, MPI_INT, pid() - 1, 0, MPI_COMM_WORLD, &status);
  #endif
  if (pid() == 0)
    fclose(fopen(name, "w"));
  fp1 =fopen(name, "a");

  // output the mesh for test
  foreach(serial) {
    double x1, y1, z1;
  #if dimension == 3
    for (int k=-1; k<2; k+=2)
  #endif
    for (int j=-1; j<2; j+=2)
      for (int i=-1; i<2; i+=2) {
        x1 = x + 0.5 * i * Delta;
        y1 = y + 0.5 * j * Delta;
      #if dimension == 3
        z1 = z + 0.5 * k * Delta;
        fprintf (fp1, "%g %g %g\n", x1, y1, z1);
      #else
        fprintf (fp1, "%g %g\n", x1, y1);
      #endif
      }
  }

  fflush(fp1);
  fclose(fp1);
  #if _MPI
  if (pid() + 1 < npe())
    MPI_Send(&sign_mpi, 1, MPI_INT, pid() + 1, 0, MPI_COMM_WORLD);
  #endif
}

#if dimension == 3
double volume, volume0;
#endif
FILE *fp_debug;

event defaults (i=0) {
#if MYDEBUG
  char name[80];
  sprintf (name, "output_tmp/debug_log_ebit.dat");
  fp_debug = fopen (name, "w");
#endif
}

void debug_log(int i) {
  if (pid() == 0) {
  #if dimension == 3
    fprintf (fp_debug, "i: %d, dt:%g vol:%.9f\n", i, dt, volume);
  #else
    fprintf (fp_debug, "i: %d, dt:%g\n", i, dt);
  #endif
    fflush (fp_debug);
  }

}

event profile_semushin (t = end) {
#if MYDEBUG
  int nc = grid -> n, tnc = grid -> tn;
  printf("\n# Total number of (leaf) cells. This process (PID:%d): %d.\
    All processes: %d\n", pid(), nc, tnc);
  fclose(fp_debug);
#endif
}

/** Used for vertex vector */
void boundary_edge(vector v_e, vector v_buf) {
  foreach_dimension()
    v_buf.x.dirty = false;

  foreach_edge() {
      v_buf.x[] = v_e.x[];
    }
    boundary ((scalar *){v_e});
    foreach_edge() {
      if (v_buf.x[] != v_e.x[])
        v_e.x[] = v_buf.x[];
    }
}