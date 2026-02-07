/**
# 3D Edge Based Interface Tracking solver

We wish to advect an interface identified by marker points that are constrained
to move only along the grid lines. There is at most one markers per edge. */

static scalar * _interface = NULL;

extern vector * u_ebit;

/* To avoid the advection of "f" */ 
event vof (i++) {
  interfaces = _interface;
}

#include "fitting-ebit.h"
#include "vof.h"
#include "utils-ebit.h"

#if TREE
#include "ebit-tree.h"
#endif

// To avoid the advection of "f"
event vof (i++) {
  _interface = interfaces;
  interfaces = NULL;
}

// Be careful !!! These become the vertex vector when we use them inside foreach_edge or foreach_vertex
// We should take care of resolution boundary condition
vector s[], snew[], ss_tmp[], with_marker[];

// color_pha_fcen[]: central color vertex,
// config_fdict[]: topology configuration dictionary
face vector color_pha_fcen[], config_fdict[];

/** Function for updating the central color point */
foreach_dimension()
  static void update_color_cen_x() {
    scalar sdir = ss_tmp.x;
    foreach_face() {
      int ind = (int) (color_pha_new[] + color_pha_new[0,1] \
        + color_pha_new[0,1,1] + color_pha_new[0,0,1]);
      int issame = ((int) color_pha_new[] == (int) color_pha_new[0, 1])\
        || ((int) color_pha_new[] == (int) color_pha_new[0,0,1]);

      if (ind > 2)
        color_pha_fcen.x[] = 1.;
      else if (ind < 2)
        color_pha_fcen.x[] = 0.;
      else if (ind == 2 && issame)
        color_pha_fcen.x[] = -1;
    }

    //face x, not fully correct for face_x
    foreach_face(x) {
      int ind, ind_diff = 0, issame, conf_old;
      ind = (int) (color_pha_new[] + color_pha_new[0,1] \
        + color_pha_new[0,1,1] + color_pha_new[0,0,1]);
      issame = ((int) color_pha_new[] == (int) color_pha_new[0,1])\
        || ((int) color_pha_new[] == (int) color_pha_new[0,0,1]);
      conf_old = (int) config_fdict.x[];

      if (ind == 2 && !issame) {
        if (ind_diff == 3) {
          // for new conf = 7 or 8, we should make sure that there is no cell with
          // conf = 7, 8 in the initial state, otherwise the code will fail
          // see the slide for all possible changes of the central color vertex
          // from 3, 4, 5, 6 to 7, 8, add this later
        }
        else if (ind_diff == 2 && conf_old == 1) {
          // from 1 to 7, 8, is it correct for x face?
          int is, ischange = 0, low_diff, conf_nei_1, conf_nei_2;
          double color1, color2;
          low_diff = (int) (fabs(color_pha_new[] - color_pha[]));
          is = (low_diff > 0) ? 0 : 1;
          conf_nei_1 = config_fdict.z[0,0,is];
          conf_nei_2 = config_fdict.z[1,0,is];
          ischange = (conf_nei_1 == 5) || (conf_nei_1 == 8) \
            || (conf_nei_2 == 6) || (conf_nei_2 == 7);
          color1 = color_pha_new[];
          color2 = color_pha_new[0,1];

          color_pha_fcen.x[] = ischange ? color2 : color1;
          printf ("complex rule for ind_diff %d: conf_old: %d, color_pha_cen: %g\n", ind_diff, conf_old, color_pha_fcen.x[]);
        }
        else if (ind_diff == 2 && conf_old == 2) {
          // from 2 to 7, 8, is it correct for x face?
          int is, ischange = 0, low_diff, conf_nei_1, conf_nei_2;
          double color1, color2;
          low_diff = (int) (fabs(color_pha_new[] - color_pha[]));
          is = (low_diff > 0) ? 0 : 1;
          conf_nei_1 = config_fdict.y[0,is,0];
          conf_nei_2 = config_fdict.y[1,is,0];
          ischange = (conf_nei_1 == 5) || (conf_nei_1 == 8) \
            || (conf_nei_2 == 4) || (conf_nei_2 == 7);

          color1 = color_pha_new[];
          color2 = color_pha_new[0,0,1];

          color_pha_fcen.x[] = ischange ? color2 : color1;
          printf ("complex rule for ind_diff %d: conf_old: %d, color_pha_cen: %g\n", ind_diff, conf_old, color_pha_fcen.x[]);
        }
        else if (ind_diff == 4) {
          // from 7, 8 to 7, 8, add this later
        }
      }
    }

    foreach_face(y) {
      int ind, ind_diff, issame, conf_old;
      ind = (int) (color_pha_new[] + color_pha_new[0,0,1] \
        + color_pha_new[1,0,1] + color_pha_new[1]);
      issame = ((int) color_pha_new[] == (int) color_pha_new[0,0,1])\
        || ((int) color_pha_new[] == (int) color_pha_new[1]);
      conf_old = (int) config_fdict.y[];

      ind_diff = 0; //test only

      if (ind == 2 && !issame) {
        if (ind_diff == 3) {
          // for new conf = 7 or 8, we should make sure that there is no cell with
          // conf = 7, 8 in the initial state, otherwise the code will fail
          // see the slide for all possible changes of the central color vertex
          // from 3, 4, 5, 6 to 7, 8, add this later
        }
        else if (ind_diff == 2 && conf_old == 2) {
          // from 1 to 7, 8, is it correct for x face?
          int ischange = 0, low_diff;
          low_diff = (int) (fabs(color_pha_new[] - color_pha[]) \
            + fabs(color_pha_new[1] - color_pha[1]));
          ischange = (sdir[] < 0. && low_diff != 0) || (sdir[0,0,1] > 1. && low_diff == 0);
          color_pha_fcen.y[] = ischange ? color_pha[] : color_pha[1];
          printf ("complex rule for ind_diff %d: conf_old: %d, color_pha_cen: %g\n", ind_diff, conf_old, color_pha_fcen.y[]);
        }
        else if (ind_diff == 2 && conf_old == 1) {
          // from 2 to 7, 8, is it correct for x face?
          int is, left_diff, conf_nei;
          left_diff = (int) (fabs(color_pha_new[] - color_pha[]) \
            + fabs(color_pha_new[0,0,1] - color_pha[0,0,1]));
          double color1 = color_pha_new[];
          double color2 = color_pha_new[1];

          is = (left_diff == 0) ? 1 : -1;
          conf_nei = config_fdict.y[is,0,0];
          color_pha_fcen.y[] = (conf_nei == 7) ? color2 : color1;
          printf ("complex rule for ind_diff %d: conf_old: %d, color_pha_cen: %g\n", ind_diff, conf_old, color_pha_fcen.y[]);
        }
        else if (ind_diff == 4) {
          // from 7, 8 to 7, 8, add this later
        }
      }
    }

    foreach_face(z) {
      int ind, ind_diff, issame, conf_old;
      ind = (int) (color_pha_new[] + color_pha_new[1] \
        + color_pha_new[1,1] + color_pha_new[0,1]);
      issame = ((int) color_pha_new[] == (int) color_pha_new[1])\
        || ((int) color_pha_new[] == (int) color_pha_new[0,1]);
      conf_old = (int) config_fdict.z[];

      ind_diff = 0; //test only

      if (ind == 2 && !issame) {
        if (ind_diff == 3) {
          // for new conf = 7 or 8, we should make sure that there is no cell with
          // conf = 7, 8 in the initial state, otherwise the code will fail
          // see the slide for all possible changes of the central color vertex
          // from 3, 4, 5, 6 to 7, 8, add this later
          
        }
        else if (ind_diff == 2 && conf_old == 1) {
          // from 1 to 7, 8, is it correct for x face?
          int ischange = 0, low_diff;
          low_diff = (int) (fabs(color_pha_new[] - color_pha[]) \
          + fabs(color_pha_new[1] - color_pha[1]));
          ischange = (sdir[] < 0. && low_diff != 0) || (sdir[0,1] > 1. && low_diff == 0);
          color_pha_fcen.y[] = ischange ? color_pha[] : color_pha[1];
          printf ("complex rule for ind_diff %d: conf_old: %d, color_pha_cen: %g\n", ind_diff, conf_old, color_pha_fcen.z[]);
        }
        else if (ind_diff == 2 && conf_old == 2) {
          // from 2 to 7, 8, is it correct for x face?
          int is, left_diff, conf_nei;
          double color1, color2;
          left_diff = (int) (fabs(color_pha_new[] - color_pha[]) \
            + fabs(color_pha_new[0,1] - color_pha[0,1]));
          color1 = color_pha_new[];
          color2 = color_pha_new[1];

          is = (left_diff == 0) ? 1 : -1;
          conf_nei = config_fdict.z[is,0,0];
          color_pha_fcen.z[] = (conf_nei == 8) ? color2 : color1;
          printf ("complex rule for ind_diff %d: conf_old: %d, color_pha_cen: %g\n", ind_diff, conf_old, color_pha_fcen.z[]);
        }
        else if (ind_diff == 4) {
          // from 7, 8 to 7, 8, add this later
        }
      }
    }
  }


static void update_color_fcen() {
  foreach_face() {
    int p1 = (int) color_pha_new[];
    int p2 = (int) color_pha_new[0,1];
    int p3 = (int) color_pha_new[0,1,1];
    int p4 = (int) color_pha_new[0,0,1];

    int ind = p1 + p2 + p3 + p4;
    int issame = (p1 == p2) || (p1 == p4);

    if (ind > 2)
      color_pha_fcen.x[] = 1.;
    else if (ind < 2)
      color_pha_fcen.x[] = 0.;
    else if (ind == 2) {
      if (issame)
        color_pha_fcen.x[] = -1;
    }
  }
  
}

/** Update the topology dictionary based on the color vertex.
 Be careful inside the foreach_dimension iterator.
 We need different results for different directions, so we use foreach_dimension*/
foreach_dimension()
  static double get_cdict_x (Point point) {
    int p1, p2, p3, p4, pc;
    double conf = 0.; // Without interface

    p1 = (int) color_pha[];
    p2 = (int) color_pha[0,1];
    p3 = (int) color_pha[0,1,1];
    p4 = (int) color_pha[0,0,1];
    pc = (int) (color_pha_fcen.x[]);

    conf = get_dict_ind (p1, p2, p3, p4, pc);
    return conf;
  }

// all the three functions produce the same result, which is different to that in 2D
foreach_dimension()
  static void update_dict_x() {
    foreach_face()
      config_fdict.x[] = get_cdict_x (point);

    boundary ((scalar *) {config_fdict});
  }


event defaults (i = 0) {
  // It's not consistent with the boundary condition of ss_tmp
  // This setting means that we only use one circle fit for makers
  // within the first layer of physical boundary
  config_fdict.t[left] = 0.;
  config_fdict.r[left] = 0.;
  config_fdict.t[right] = 0.;
  config_fdict.r[right] = 0.;

  config_fdict.t[top] = 0.;
  config_fdict.r[top] = 0.;
  config_fdict.t[bottom] = 0.;
  config_fdict.r[bottom] = 0.;

  config_fdict.t[front] = 0.;
  config_fdict.r[front] = 0.;
  config_fdict.t[back] = 0.;
  config_fdict.r[back] = 0.;

  foreach_dimension() {
    ss_tmp.x.nodump = true;
    snew.x.nodump = true;
    s.x.nodump = true;
    with_marker.x.nodump = true;
  }
}


int get_polygons (Point point, coord *xpoly, int *nv_poly) {
  int cv, flag_edge[12], polygon[4][12], config_f[6], color_v[8];

  cv = 0;
  for (int kk = 0; kk < 2; kk++)
    for (int jj = 0; jj < 2; jj++)
      for (int ii = 0; ii < 2; ii++) {
        int iv = ii + jj*2 + kk*4;
        color_v[iv] = (int) color_pha[ii,jj,kk];
        cv += color_v[iv];
      }
  if (cv == 8 || cv == 0)
    return 0;

  int iface = 0;
  foreach_dimension() {
    config_f[iface] = (int) config_fdict.x[];
    config_f[iface + 1] = (int) config_fdict.x[1];
    iface += 2;
  }

  for (int ii = 0; ii < 12; ii++) {
    flag_edge[ii] = 0;
    polygon[0][ii] = polygon[1][ii] = -1;
    polygon[2][ii] = polygon[3][ii] = -1;
  }

  int npoly = 0;
  for (int ii = 0; ii < 6; ii++) {
    int conf = config_f[ii];

    if (conf > 0) { // face with markers
      int ind_markers[4], ind_f;
      ind_f = ii; // initial face index

      for (int iind = 0; iind < 4; iind++)
        ind_markers[iind] = dict_to_edge[conf][iind];

      for (int iv = 0; iv < 4; iv += 2) {
        int with_new = 0;
        if (ind_markers[iv] != -1) {
          int ie1 = ind_markers[iv], ie2 = ind_markers[iv + 1];
          // we have found two markers, the following part will find out the rest
          for (int ie_poly = 0; ie_poly < 12; ie_poly++) {
            int ie_cube = edge_face_to_cube[ind_f][ie1];
            if (flag_edge[ie_cube] == 0) {

              int ind_fn, conf_n;
              with_new = 1; // find a new polygon

              polygon[npoly][ie_poly] = ie_cube;
              flag_edge[ie_cube] = 1;
              // edge on next cell face
              ind_fn = face_to_face[ind_f][ie2][3];
              ie1 = face_to_face[ind_f][ie2][4];

              conf_n = config_f[ind_fn];
              ie2 = dict_to_markers[conf_n][ie1];
              ind_f = ind_fn;
            }
            else {
              // we have found all the edges of a closed polygon
              break;
            }
          }
        }

        if (with_new)
          npoly++;
      }
    }
  } // finish searching six cell faces

  // local coordinates of vertex of the polygon

  for (int ip = 0; ip < npoly; ip++) {
    int nv = 12;

    for (int iv = 0; iv < 12; iv++) {
      int ie = polygon[ip][iv], ind_v = 12*ip + iv;
      if (ie >= 0) {
        int ii, jj, kk, idir, ds[3] = {0, 0, 0};
        idir = ie/4;
        ds[idir] = 1;

        ii = edge_shift[ie][0];
        jj = edge_shift[ie][1];
        kk = edge_shift[ie][2];
        xpoly[ind_v].x = ii + ds[0]*s.x[ii,jj,kk];
        xpoly[ind_v].y = jj + ds[1]*s.y[ii,jj,kk];
        xpoly[ind_v].z = kk + ds[2]*s.z[ii,jj,kk];
      }
      else {
        nv = iv;
        break;
      }
    }
    nv_poly[ip] = nv;
  }
  return npoly;
}


/**
 ## Output functions
 
 * Output interface facets. if tri == true, the polygon will be decomposed into
 * several triangles, and the connecticity is stored in the file_con.
 */
void output_facets_semushin (char *file, bool tri = false, \
  char *file_con = "connectivity.dat") {
  char name[strlen(file) + 2], name_con[strlen(file_con) + 2];
  strcpy (name, file);
  strcpy (name_con, file_con);
  update_dict_x();

  FILE *fp1, *fp2;
  int sign_mpi[2], nvertex, nelement;

  #if _MPI
  MPI_Status status;
  if (pid() != 0)
    MPI_Recv (&sign_mpi, 2, MPI_INT, pid() - 1, 0, MPI_COMM_WORLD, &status);
  #endif
  if (pid() == 0) {
    fclose (fopen (name, "w"));
    if (tri)
      fclose (fopen (name_con, "w"));
  }

  fp1 = fopen (name, "a");
  if (tri)
    fp2 = fopen (name_con, "a");

  // use three scalar to avoid index rotation.
  scalar sx = s.x, sy = s.y, sz = s.z;
  if (!tri) {
    int idim = -2;
    coord dx = {0., 1., 1.};

    foreach_dimension() {
      idim += 2;
      foreach_face(x, serial, noauto) {
        int ind_markers[4] = {-1, -1, -1, -1};
        double xo[3] = {x, y, z};

        // coordinate of origin, vertex [0, 0, 0]
        // we need to exchange the dx.y and dx.y, because that in foreach_dimension
        // x, y ,z -> y, z, x -> z, x, y
        xo[0] -= 0.5*Delta*dx.x;
        xo[1] -= 0.5*Delta*dx.z;
        xo[2] -= 0.5*Delta*dx.y;

        // Start calculating the connectivity of markers based on the
        // color vertex (central point)
        int conf = (int) config_fdict.x[];

        for (int iind = 0; iind < 4; iind++) {
          ind_markers[iind] = dict_to_edge[conf][iind];
        }

        for (int iv = 0; iv < 4; iv += 2) {
          if (ind_markers[iv] != -1) {
            double xx[2], yy[2], zz[2];
            int ie_locs[2] = {ind_markers[iv], ind_markers[iv + 1]};

            for (int iie = 0; iie < 2; iie++) {
              int ie_loc, ie_cube, is, js, ks, idir;
              double dd;
              ie_loc = ie_locs[iie];
              ie_cube = edge_face_to_cube[idim][ie_loc];
              
              idir = ie_cube/4;
              is = edge_shift[ie_cube][0];
              js = edge_shift[ie_cube][1];
              ks = edge_shift[ie_cube][2];

              xx[iie] = xo[0] + is*Delta;
              yy[iie] = xo[1] + js*Delta;
              zz[iie] = xo[2] + ks*Delta;

              if (idir == 0) {
                dd = get_scalar (point, sx, is, js, ks);
                xx[iie] += dd*Delta;
              }
              else if (idir == 1) {
                dd = get_scalar (point, sy, is, js, ks);
                yy[iie] += dd*Delta;
              }
              else {
                dd = get_scalar (point, sz, is, js, ks);
                zz[iie] += dd*Delta;
              }
            }

            fprintf (fp1, "%g %g %g %d\n%g %g %g %d\n\n", xx[0], yy[0], zz[0], pid(), xx[1], yy[1], zz[1], pid());
          }
        }
      }
    }
  }
  else {
    // nvertex, nelement: numbers of vertices and triangular elements from previous process
    if (pid() == 0) {
      nvertex = nelement = 0;
    }
    else {
      nvertex = sign_mpi[0];
      nelement = sign_mpi[1];
    }
    foreach(serial, noauto) {
      coord xo = {x - 0.5*Delta, y - 0.5*Delta, z - 0.5*Delta};
      coord xpoly[48];
      int npoly, nv_poly[4];
      npoly = get_polygons(point, xpoly, nv_poly);
      for (int ip = 0; ip < npoly; ip++) {
        int nv = nv_poly[ip];
        double xc = 0, yc = 0, zc = 0;
        for (int iv = 0; iv < nv; iv++) {
          double xx, yy, zz;
          xx = xo.x + xpoly[12*ip + iv].x*Delta;
          yy = xo.y + xpoly[12*ip + iv].y*Delta;
          zz = xo.z + xpoly[12*ip + iv].z*Delta;
          xc += xx; yc += yy; zc += zz;
          fprintf (fp1, "%g %g %g %d\n", xx, yy, zz, pid());
        }
        xc /= nv; yc /= nv; zc /= nv;

        if (nv > 3) {
          fprintf (fp1, "%g %g %g %d\n", xc, yc, zc, pid());
          for (int iv = 0; iv < nv; iv++) {
            fprintf (fp2, "%d %d %d\n", nvertex + iv, nvertex + (iv + 1) % nv, nvertex + nv);
          }
          nvertex += (nv + 1);
          nelement += nv;
        }
        else {
          fprintf (fp2, "%d %d %d\n", nvertex, nvertex + 1, nvertex + 2);
          nvertex += nv;
          nelement++;
        }

      }
    }
  }

  fflush (fp1);
  fclose (fp1);

  if (tri) {
    fflush (fp2);
    fclose (fp2);
    sign_mpi[0] = nvertex;
    sign_mpi[1] = nelement;
  }
  #if _MPI
  if (pid() + 1 < npe())
    MPI_Send (&sign_mpi, 2, MPI_INT, pid() + 1, 0, MPI_COMM_WORLD);
  #endif
}


/** Output interface, color vertex and mesh. Since it is very time consuming, 
 * it should be used for debugging purpose only.
*/
#if _MYOUTPUT
static void output_intf (int i) {
  char out[100], out_test[100], out_color[100];
  const char testName[50] = TEST;

  sprintf (out, "%s_%d_%d.dat", testName, N, i); 
  sprintf (out_color, "%s_%d_%d.dat", "output_debug/color_result", N, i);
  sprintf (out_test, "%s_%d.dat", "output_debug/grid", i);

  output_facets_semushin (out);
  // output_color_vertex (file=out_color, cv=color_pha, cv_fcenter=color_pha_fcen);
  // output_mesh (out_test);
}
#endif

/** Initialize the marker position based on the level-set function "phi". */
static void init_markers (vertex scalar phi) {
  boundary ({phi});
  fractions (phi, f); // we need this to for the stability events.
  // copy from fraction.h
  double val = 0.;
  foreach_edge_i() {
    if ((phi[] - val)*(phi[1] - val) < 0.) {
      snew.x[] = (phi[] - val)/(phi[] - phi[1]);
      if (phi[] < val)
	      snew.x[] = 1. - snew.x[];
    }
    else
      snew.x[] = (phi[] > val || phi[1] > val);
  }

  boundary ((scalar *) {snew});

  // Initialize the color vertex
  foreach_vertex() {
    if (phi[] >= 0.) // remember the equal sign
      color_pha[] = 1.; // Reference phase
    else
      color_pha[] = 0.;

    color_pha_new[] = color_pha[];
  }

  // For the color_pha_cen[] is undefined (-1) when only two consecutive color vertices
  // are the same color. When more then two vertices are the same color, we change the
  // color of central color point.
  foreach_face()
    color_pha_fcen.x[] = -1.;

  update_color_fcen();

  foreach_edge_b() {
    s.x[] = 0.;
    with_marker.x[] = 0.;
  }

  // Different values of color vertex at the two end of edge means there is a marker on the edge.
  // This method can deal with the corner case (interface pass through the cell vertex) correctly.
  // with_marker : number of markers on each edge
  boundary ((scalar *) {color_pha});
  foreach_edge_i() {
    int with_face = fabs(color_pha[] - color_pha[1]) > machine_zero;
    if (with_face) {
      if (color_pha[] > machine_zero)
        s.x[] = snew.x[];
      else
        s.x[] = 1. - snew.x[];

      with_marker.x[] = 1.;
    }
  }

  boundary ((scalar *) {s});
  // output the initial interface for test
  #if _MYDEBUG
  output_intf (-1);
  #endif

  #if TREE
  // Do not put this part into event defaults, there is something wrong with
  // the vertex vector field, the methods are not set correctly.
  // must set both prolongation, refine, coarsen, restriction manually
  foreach_dimension() {
    s.x.refine = my_no_restriction_r;
    s.x.prolongation = my_no_restriction_p;
    s.x.coarsen = my_no_restriction_r;
    s.x.restriction = my_restriction_vertex;

    color_pha_fcen.x.prolongation = color_pha_fcen.x.refine = no_restriction;
    color_pha_fcen.x.restriction = color_pha_fcen.x.coarsen = no_restriction;

    // we need the refine function to make MPI transfering the data
    // between the ghost cells. Do not use no_restriction for
    // the refine function
    config_fdict.x.prolongation = no_restriction;
    config_fdict.x.restriction = config_fdict.x.coarsen = no_restriction;
  }

  color_pha.restriction = my_restriction_vertex;
  color_pha.coarsen = my_restriction_vertex;
  color_pha.prolongation = prolongation_vertex;
  color_pha.refine = prolongation_vertex;

  // central color vertex is refine in color_pha.refine
  // we still need restriction
  color_pha_fcen.x.restriction = restriction_face;
  color_pha_fcen.x.coarsen = restriction_face;
  color_pha_fcen.x.refine = refine_face;
  foreach_dimension() {
    color_pha_fcen.x.prolongation = refine_face_injection_x;
  }

  #endif
}


void set_markers() {
  // determine the number of markers on edge based on color vertex, this is a robust method for corner case,
  // it can be generalized to double-Semushin easily in the furture.
  foreach_edge_i() {
    int with_face = fabs(color_pha[] - color_pha[1]) > machine_zero;
    if (with_face) {
      with_marker.x[] = 1.;
    }
    else { // without interface or with two interface
      with_marker.x[] = 0.;
      s.x[] = 0.;
    }
  }
}


/**
 ## Front2VOF algorithm

The representation of the interface with EBIT markers can be used to create
the associated VOF field. The connectivity of markers and the region of reference phase are 
idenfied by the "color vertex". 
The Front-to-VOF algorithm is employed to calculate the volume of reference*/

// How to determine the right sign of the normal direction
double f2vof (coord *xpoly, int npoly, int *nv_poly, int *color_v) {
  double frac = 0.;
  coord f_p = {0., 0., 0.}, fs_p = {0., 0., 0.};
  double ds_total_cell = 0.;
  coord dn_ave_cell = {0., 0., 0.};

  for (int ipoly = 0; ipoly < npoly; ipoly++) {
    coord xcen = {0., 0., 0.}, xx_poly[13];
    int np = nv_poly[ipoly];

    for (int ip = 0; ip < np; ip++) {
      foreach_dimension() {
        xx_poly[ip].x = xpoly[12*ipoly + ip].x;
        xcen.x += xx_poly[ip].x;
      }
    }

    foreach_dimension() {
      xx_poly[np].x = xx_poly[0].x;
      xcen.x /= np;
    }

    // calculate the average normal vector
    double ds_total = 0.;
    coord dn_ave = {0., 0., 0.};

    for (int ip = 0; ip < np; ip++) {
      coord dx[3], dn;
      double x1, x2, x3;

      foreach_dimension() {
        x1 = xx_poly[ip].x;
        x2 = xx_poly[ip + 1].x;
        x3 = xcen.x;

        dx[0].x = x2 - x1;
        dx[1].x = x3 - x2;
        dx[2].x = x1 - x3;
      }

      double ds_tri = 0.;
      foreach_dimension() {
        dn.x = 0.5*(dx[1].y*dx[2].z - dx[1].z*dx[2].y);
        dn_ave.x += dn.x;
        ds_tri += sq(dn.x);
      }
      ds_tri = sqrt(ds_tri);
      ds_total += ds_tri;
    }

    ds_total = max(ds_total, 1.e-32);
    foreach_dimension()
      dn_ave.x /= ds_total;

    if (ds_total > ds_total_cell) {
      ds_total_cell = ds_total;
      foreach_dimension()
        dn_ave_cell.x = dn_ave.x;
    }
    // identify the reference phase
    // not correct for corner case, fix this
    int iref[2] = {-1, -1}, ind_ref = 0;
    double sig_n = 1.;
    for (int kk = 0; kk < 2; kk++)
      for (int jj = 0; jj < 2; jj++)
        for (int ii = 0; ii < 2; ii++) {
          int iv, icol, iside;
          coord xv = {ii, jj, kk};
          double sig = 0.;
          iv = ii + jj*2 + kk*4;
          icol = color_v[iv];
          foreach_dimension()
            sig += (xv.x - xcen.x)*dn_ave.x;

          iside = sig > 0. ? 0 : 1;

          if (iref[iside] < 0)
            iref[iside] = icol;
          else if (iref[iside] != icol)
            ind_ref = 1 - iside;
        }

    sig_n = (ind_ref == iref[ind_ref]) ? 1. :-1.;

    // start front-to-vof
    for (int ip = 0; ip < np; ip++) {
      coord dx[3], p1, p2, dn, xc;
      double x1, x2, x3;

      foreach_dimension() {
        x1 = xx_poly[ip].x;
        x2 = xx_poly[ip + 1].x;
        x3 = xcen.x;
        if (sig_n < 0)
          swap(double, x1, x2);
        p1.x = x1;
        p2.x = x2;

        dx[0].x = x2 - x1;
        dx[1].x = x3 - x2;
        dx[2].x = x1 - x3;

        xc.x = (x1 + x2 + x3)/3.;
      }

      foreach_dimension() {
        if (fabs(p1.x - 1.) < 1.e-32 && fabs(p2.x - 1.) < 1.e-32) {
          double dlx = -dx[0].z;
          double dly = dx[0].y;
          double zc = 0.5*(p1.z + p2.z);
          double yc = 0.;
          if (fabs(p1.z - 1.) < 1.e-32)
            yc = (dlx >= 0) ? p1.y : 1. - p1.y;
          if (fabs(p2.z - 1.) < 1.e-32)
            yc = (dlx >= 0) ? p2.y : 1. - p2.y;

          fs_p.x += zc*dly + yc;
        }
      }

      foreach_dimension() {
        dn.x = 0.5*(dx[1].y*dx[2].z - dx[1].z*dx[2].y);
        f_p.x += dn.x*xc.x;
      }
    }
  }

  coord m2 = {fabs(dn_ave_cell.x), fabs(dn_ave_cell.y), fabs(dn_ave_cell.z)};
  int dir = (m2.x > m2.y) ? (m2.x > m2.z ? 0 : 2) : (m2.y > m2.z ? 1 : 2);

  int idir = 0;
  foreach_dimension() {
    if (idir == dir) {
      fs_p.x = fs_p.x < 0. ? 1 + fs_p.x : fs_p.x; // not correct for corner case
      f_p.x += fs_p.x;
      f_p.x = f_p.x < 0. ? 1 + f_p.x : f_p.x;
      frac = f_p.x;
    }
    idir++;
  }

  return frac;
}

/** 
 ## Volume fraction computation
*/
void semu2vof() {
  scalar with_faces[];
  update_dict_x();

  foreach() {
    int cv1 = (int) (color_pha[] + color_pha[1] + color_pha[1,1] + color_pha[0,1] \
      + color_pha[0,0,1] + color_pha[1,0,1] + color_pha[1,1,1] + color_pha[0,1,1]);

    with_faces[] = (cv1 == 8 || cv1 == 0) ? 0. : 1.;

    // cells without interface
    if (cv1 == 8)
      f[] = 1.;
    else
      f[] = 0.;
  }

  foreach(noauto) {
    int with_f = (int) (with_faces[]);

    if (with_f != 0) { // interfacial cell, cv is used for identifying the reference phase
      // find out the polygons within the cell (at most 2)
      int flag_edge[12], polygon[4][12], config_f[6], iface, color_v[8];

      for (int kk = 0; kk < 2; kk++)
        for (int jj = 0; jj < 2; jj++)
          for (int ii = 0; ii < 2; ii++) {
            int iv = ii + jj*2 + kk*4;
            color_v[iv] = (int) color_pha[ii,jj,kk];
          }

      iface = 0;
      foreach_dimension() {
        config_f[iface] = (int) config_fdict.x[];
        config_f[iface + 1] = (int) config_fdict.x[1];
        iface += 2;
      }

      for (int ii = 0; ii < 12; ii++) {
        flag_edge[ii] = 0;
        polygon[0][ii] = polygon[1][ii] = -1;
        polygon[2][ii] = polygon[3][ii] = -1;
      }

      int npoly = 0;
      for (int ii = 0; ii < 6; ii++) {
        int conf = config_f[ii];

        if (conf > 0) { // face with markers
          int ind_markers[4], ind_f;
          ind_f = ii; // initial face index

          for (int iind = 0; iind < 4; iind++)
            ind_markers[iind] = dict_to_edge[conf][iind];
          
          for (int iv = 0; iv < 4; iv += 2) {
            int with_new = 0;
            if (ind_markers[iv] != -1) {
              int ie1 = ind_markers[iv], ie2 = ind_markers[iv + 1];
              // we have found two markers, the following part will find out the rest
              for (int ie_poly = 0; ie_poly < 12; ie_poly++) {
                int ie_cube = edge_face_to_cube[ind_f][ie1];
                if (flag_edge[ie_cube] == 0) {

                  int ind_fn, conf_n;
                  with_new = 1; // find a new polygon

                  polygon[npoly][ie_poly] = ie_cube;
                  flag_edge[ie_cube] = 1;
                  // edge on next cell face
                  ind_fn = face_to_face[ind_f][ie2][3];
                  ie1 = face_to_face[ind_f][ie2][4];

                  conf_n = config_f[ind_fn];
                  ie2 = dict_to_markers[conf_n][ie1];
                  ind_f = ind_fn;
                }
                else {
                  // we have found all the edges of a closed polygon
                  break;
                }
              }
            }

            if (with_new)
              npoly++;
          }
        }
        
      } // finish searching six cell faces

      // calculate the volume fractions
      // local coordinates of vertex of the polygon
      double fpoly = 0.;
      // 0-11 for the first polygon, the rest for the second one, at most 12 vertices for each polygon
      coord xpoly[48];
      // number of polygons within the cell
      int nv_poly[4];

      for (int ip = 0; ip < npoly; ip++) {
        int nv = 12;

        for (int iv = 0; iv < 12; iv++) {
          int ie = polygon[ip][iv], ind_v = 12*ip + iv;
          if (ie >= 0) {
            int ii, jj, kk, idir, ds[3] = {0, 0, 0};
            idir = ie/4;
            ds[idir] = 1;

            ii = edge_shift[ie][0];
            jj = edge_shift[ie][1];
            kk = edge_shift[ie][2];
            xpoly[ind_v].x = ii + ds[0]*s.x[ii,jj,kk];
            xpoly[ind_v].y = jj + ds[1]*s.y[ii,jj,kk];
            xpoly[ind_v].z = kk + ds[2]*s.z[ii,jj,kk];
          }
          else {
            nv = iv;
            break;
          }
        }
        nv_poly[ip] = nv;
      }

      fpoly = f2vof (xpoly, npoly, nv_poly, color_v);

      f[] = fpoly;
    }
  }


  boundary ({f});
  volume = 0.;
  volume_int = 0.;
  foreach(reduction(+:volume) reduction(+:volume_int)) {
    int with_f = (int) (with_faces[]);
    volume += f[]*cube(Delta);
    if (with_f != 0)
      volume_int += f[]*cube(Delta);
  }
}

/**
  This event perform the 1-D advection of the EBIT marker points. */
foreach_dimension()
  void advect_x (vector u, int ind = 0) {
    double idir = iadv.x;
    boundary ((scalar *) {u});

    scalar ux = u.x; // use this to get rid of rotation in foreach_dimension(), check this

    set_markers();

    update_dict_x(); // check the output of this function for _y and _z

    foreach_edge_b() {
      snew.x[] = 0.;
      ss_tmp.x[] = 0.;
    }

    // it's needed in AMR, because color_pha_new[] is not updated after the adapt_wavelet,
    foreach_vertex()
      color_pha_new[] = color_pha[];

     // This performs the horizontal advection, for both aligned and unaligned markers
    foreach_dimension()
      foreach_vertex() {
        double um, u1, xx;

        um = 0.;
        #ifdef UTEST
        um = (ux[] + ux[0,-1] + ux[0,-1,-1] + ux[0,0,-1])/4.;
        #else
        for (int ii = -1; ii < 2; ii++) {
          u1 = (ux[ii] + ux[ii,-1] + ux[ii,-1,-1] + ux[ii,0,-1])/4.;
          xx = s.x[] - (0.5 + ii*1.);
          xx = DR_L(xx);
          um += xx*u1;
        }
        #endif

        ss_tmp.x[] = um*dt/Delta;
      }

    // final position of aligned markers
    foreach_vertex() {
      if ((int) with_marker.x[] > 0) {
        double xx = ss_tmp.x[];
        ss_tmp.x[] = s.x[] + xx;
        snew.x[] = s.x[] + xx;
      }
      else {
        ss_tmp.x[] = 0.;
      }
    }

    foreach_vertex() {
      if (snew.x[] > 1. || snew.x[] < 0.) {  // The point is crossing a vertical line
        snew.x[] = 0.;
        with_marker.x[] -= 1;
      }
    }
    boundary ((scalar *) {ss_tmp});

    // When there are two markers on the same edge, snew for the old one, s_tmp for the new one
    // For the simple Semushin method, We will remove these two markers on the
    // final stage based on the color vertex.
    foreach_vertex() {
      int withi = (int) with_marker.x[], withr, withl;
      withr = ss_tmp.x[1] < 0.;
      withl = ss_tmp.x[-1] > 1.;
      if (withi == 0) {
        if (ss_tmp.x[-1] > 1.)
          snew.x[] = ss_tmp.x[-1] - 1.;
        else if (ss_tmp.x[1] < 0.)
          snew.x[] = 1. + ss_tmp.x[1];

        #if _MYDEBUG
        if (withl && withr)
          printf ("istep:%d, idir: %g, pid: %d, Double markers on edge: x:%g, y:%g, z:%g, \n\
           ii:%g, jj:%g, kk:%g \n|current:%g, left:%g, right:%g\n\n", ind, \
           idir, pid(), x, y, z, x/Delta, y/Delta, z/Delta, snew.x[], ss_tmp.x[-1], ss_tmp.x[1]);
        #endif
      }
      else {
        #if _MYDEBUG
        if (withl || withr)
          printf ("istep:%d, idir: %g, pid: %d, Double markers on edge: x:%g, y:%g, z:%g,\n\
          ii:%g, jj:%g, kk:%g \n|current:%g, left:%g, right:%g\n\n", ind, \
           idir, pid(), x, y, z, x/Delta, y/Delta, z/Delta, snew.x[], ss_tmp.x[-1], ss_tmp.x[1]);
        #endif
      }

      with_marker.x[] += (withl + withr);
    }

    // Change the color of vertex when marker moves across the grid line
    boundary ((scalar *) {color_pha});
    foreach_vertex() {
      if (ss_tmp.x[-1] > 1.)
        color_pha_new[] = color_pha[-1];
      else if (ss_tmp.x[] < 0.) {
        color_pha_new[] = color_pha[1];
      }
    }

    boundary ((scalar *) {color_pha_new});

    double kr_max = 5.;
    // unaligned markers on y direction
    s.y.dirty = false;
    snew.y.dirty = false;
    foreach_vertex() {
      double y3, y2, y1, x3, x2, x1;
      int withi = (int) with_marker.y[];
      with_marker.y[] = 0.; // should be removed

      if (fabs(ss_tmp.y[]) < machine_zero && withi) {
        // For markers at the physical boundary (for symmetric boundary)
        snew.y[] = s.y[];
        with_marker.y[] = 1.; // should be removed
      }
      else {
        for (int is = 0; is < 2; is++) {
          // is=0 for the cell on the right, is=1 for the cell on the left
          // reconnect the markers and calculate the intersection basing on
          // the topology of cells on the both side of the cell face (last time step)
          double xy_edge[4][2], y0, y0c;
          int ind_markers[4] = {-1, -1, -1, -1};

          int conf = (int) config_fdict.z[-is];
          if (conf != 0) {
            xy_edge[0][0] = ss_tmp.y[-is] - is;
            xy_edge[0][1] = s.y[-is];

            xy_edge[1][0] = ss_tmp.x[-is] - is;
            xy_edge[1][1] = 0.;

            xy_edge[2][0] = 1. + ss_tmp.y[1 - is] - is;
            xy_edge[2][1] = s.y[1 - is];

            xy_edge[3][0] = ss_tmp.x[-is, 1] - is;
            xy_edge[3][1] = 1.;

            for (int iind = 0; iind < 4; iind++)
              ind_markers[iind] = dict_to_edge[conf][iind];

            // following part should be replaced by a function when we introduce circle fit
            for (int iv = 0; iv < 4; iv += 2) {
              if (ind_markers[iv] != -1) {
                int ie1=ind_markers[iv], ie2=ind_markers[iv + 1];
                x1 = xy_edge[ie1][0];
                y1 = xy_edge[ie1][1];

                x2 = xy_edge[ie2][0];
                y2 = xy_edge[ie2][1];

                y0 = my_intersect (x1, y1, x2, y2, 0.); // y0 == -1 means no intersection

                if (y0 >= 0.) {
                  snew.y[] = y0;
                  with_marker.y[] += 1.; // for checking the numner of markers

                  #ifdef CIRCLE_FIT
                  // find out the four markers used for circle fit
                  // Current version of implementation: the intersection point is the average of two circle fits
                  int iee[2] = {ie1, ie2}, ipx, ipy, conn, ipps, ippe, nsec = 0;
                  double yave = 0., kp[2] = {0., 0.}, ycyc_int[2] = {HUGE, HUGE};
                  for (int ipe = 0; ipe < 2; ipe++) {
                    if (iee[ipe] == 0) {
                      ipx = -1; ipy = 0;
                    }
                    else if (iee[ipe] == 1) {
                      ipx = 0; ipy = -1;
                    }
                    else if (iee[ipe] == 2) {
                      ipx = 1; ipy = 0;
                    }
                    else{
                      ipx = 0; ipy = 1;
                    }
                    conn = (int) config_fdict.z[-is + ipx,ipy];
                    if (conn == 0)
                      continue;
                    ipps = 2*(iee[ipe] % 2 + 1) - iee[ipe];
                    ippe = get_end (conn, ipps);

                    if (ippe == 0) {
                      x3 = ss_tmp.y[-is + ipx,ipy] - is + ipx;
                      y3 = s.y[-is + ipx,ipy] + ipy;
                    }
                    else if (ippe == 1) {
                      x3 = ss_tmp.x[-is + ipx,ipy] - is + ipx;
                      y3 = ipy;
                    }
                    else if (ippe == 2) {
                      x3 = ss_tmp.y[-is + ipx + 1,ipy] - is + ipx + 1.;
                      y3 = s.y[-is + ipx + 1,ipy] + ipy;
                    }
                    else {
                      x3 = ss_tmp.x[-is + ipx,ipy + 1] - is + ipx;
                      y3 = 1. + ipy;
                    }

                    double xrc, yrc, rc, xmin, xmax;
                    get_circle (x1, y1, x2, y2, x3, y3, &xrc, &yrc, &rc);
                    y0c = intersection_circle (xrc, yrc, rc, 0.);

                    if (rc > 0.) {
                      // sign is not correct
                      kp[ipe] = 1./rc;
                    }
                    // If the intersection point is outside [0., 1.], we abandon the result
                    ycyc_int[ipe] = y0;
                    if (rc > 0. && y0c > 0.) {
                      xmin = min(min(x1, x2), x3);
                      xmax = max(max(x1, x2), x3);
                      // we don't use the result of extrapolation
                      if (xmin <= 0. && xmax >= 0.) {
                        ycyc_int[ipe] = y0c;
                        yave += y0c;
                        nsec++;
                      }
                    }
                  }

                  // curvature ratio of two circle fittings
                  double kr = max(fabs(kp[0]), fabs(kp[1]))/(min(fabs(kp[0]), fabs(kp[1])) + 1.e-32);
                  if (nsec == 2 && kr > kr_max) {
                    // for the region with large gradient of curvature, choose the fitting
                    // resulting in smaller curvature (or linear fitting)
                    yave = fabs(kp[0]) < fabs(kp[1]) ? nsec*ycyc_int[0]:  nsec*ycyc_int[1];
                  }

                  if (nsec > 0) {
                    // We revert to straight line fit if all the circle fits fail
                    yave /= nsec;
                    snew.y[] = yave;
                  }
                  #endif
                }
              }
            }
          }

        }
      }
    } // end foreach_edge_y, advection of unaligned markers

    // unaligned markers on z direction
    s.z.dirty = false;
    snew.z.dirty = false;
    foreach_vertex() {
      double y3, y2, y1, x3, x2, x1;
      int withi = (int) with_marker.z[];
      with_marker.z[] = 0.; // should be removed

      if (fabs(ss_tmp.z[]) < machine_zero && withi) {
        // For markers at the physical boundary (for symmetric boundary)
        snew.z[] = s.z[];
        with_marker.z[] = 1.; // should be removed
      }
      else {
        for (int is = 0; is < 2; is++) {
          double xy_edge[4][2], y0, y0c;
          int ind_markers[4] = {-1, -1, -1, -1};

          // this is why we can not use foreach_dimension,
          int conf = (int) config_fdict.y[-is];
          if (conf != 0) {
            // exchange first index 1 <-> 0, 2 <-> 3.
            xy_edge[1][0] = ss_tmp.z[-is] - is;
            xy_edge[1][1] = s.z[-is];

            xy_edge[0][0] = ss_tmp.x[-is] - is;
            xy_edge[0][1] = 0.;

            xy_edge[3][0] = 1. + ss_tmp.z[1 - is] - is;
            xy_edge[3][1] = s.z[1 - is];

            xy_edge[2][0] = ss_tmp.x[-is, 0, 1] - is;
            xy_edge[2][1] = 1.;

            for (int iind = 0; iind < 4; iind++) {
              ind_markers[iind] = dict_to_edge[conf][iind];
            }

            for (int iv = 0; iv < 4; iv += 2) {
              if (ind_markers[iv] != -1) {
                int ie1 = ind_markers[iv], ie2 = ind_markers[iv + 1];
                x1 = xy_edge[ie1][0];
                y1 = xy_edge[ie1][1];

                x2 = xy_edge[ie2][0];
                y2 = xy_edge[ie2][1];

                y0 = my_intersect(x1, y1, x2, y2, 0.);

                if (y0 >= 0.) {
                  snew.z[] = y0;
                  with_marker.z[] += 1.;

                  #ifdef CIRCLE_FIT
                  int iee[2] = {ie1, ie2}, ipx, ipy, conn, ipps, ippe, nsec = 0;
                  double yave = 0., kp[2] = {0., 0.}, ycyc_int[2] = {HUGE, HUGE};;
                  for (int ipe = 0; ipe < 2; ipe++) {
                    if (iee[ipe] == 0) {
                      ipx = 0; ipy = -1;
                    }
                    else if (iee[ipe] == 1) {
                      ipx = -1; ipy = 0;
                    }
                    else if (iee[ipe] == 2) {
                      ipx = 0; ipy = 1;
                    }
                    else {
                      ipx = 1; ipy = 0;
                    }
                    conn = (int) config_fdict.y[-is + ipx,0,ipy];
                    if (conn == 0)
                      continue;
                    ipps = 2*(iee[ipe] % 2 + 1) - iee[ipe];
                    ippe = get_end (conn, ipps);

                    if (ippe == 1) {
                      x3 = ss_tmp.z[-is + ipx,0,ipy] - is + ipx;
                      y3 = s.z[-is + ipx,0,ipy] + ipy;
                    }
                    else if (ippe == 0) {
                      x3 = ss_tmp.x[-is + ipx,0,ipy] - is + ipx;
                      y3 = ipy;
                    }
                    else if (ippe == 3) {
                      x3 = ss_tmp.z[-is + ipx + 1,0,ipy] - is + ipx + 1.;
                      y3 = s.z[-is + ipx + 1,0,ipy] + ipy;
                    }
                    else {
                      x3 = ss_tmp.x[-is + ipx,0,ipy + 1] - is + ipx;
                      y3 = 1. + ipy;
                    }

                    double xrc, yrc, rc, xmin, xmax;
                    get_circle (x1, y1, x2, y2, x3, y3, &xrc, &yrc, &rc);
                    y0c = intersection_circle (xrc, yrc, rc, 0.);

                    if (rc > 0.) {
                      // sign is not correct
                      kp[ipe] = 1./rc;
                    }
                    // If the intersection point is outside [0., 1.], we abandon the result
                    ycyc_int[ipe] = y0;
                    if (rc > 0. && y0c > 0.) {
                      xmin = min(min(x1, x2), x3);
                      xmax = max(max(x1, x2), x3);
                      if (xmin <= 0. && xmax >= 0.) {
                        ycyc_int[ipe] = y0c;
                        yave += y0c;
                        nsec++;
                      }
                    }
                  }

                  // curvature ratio of two circle fittings
                  double kr = max(fabs(kp[0]), fabs(kp[1]))/(min(fabs(kp[0]), fabs(kp[1])) + 1.e-32);
                  if (nsec == 2 && kr > kr_max) {
                    // for the region with large gradient of curvature, choose the fitting
                    // resulting in smaller curvature (or linear fitting)
                    yave = fabs(kp[0]) < fabs(kp[1]) ? nsec*ycyc_int[0]:  nsec*ycyc_int[1];
                  }

                  if (nsec > 0) {
                    yave /= nsec;
                    snew.z[] = yave;
                  }
                  #endif
                }
              }
            }
          }

        }
      }
    } // end foreach_edge_z, advection of unaligned markers

    // Update the color vertex, don't change the function call order
    update_color_cen_x();
    foreach_vertex()
      color_pha[] = color_pha_new[];

    boundary ((scalar *) {color_pha});
    // determine the number of markers on edge based on color vertex, robust method for corner case,
    // it can be generalized to double-Semushin easily in the furture.
    // without noauto, it trigger the prolongation of color_pha, then onsameside give a floating
    // point exception, should check this later !!!
    foreach_edge_i() {
      int with_face = fabs(color_pha[] - color_pha[1]) > machine_zero;
      if (with_face) {
        with_marker.x[] = 1.;
      }
      else { // without interface or with two interface
        snew.x[] = 0.;
        with_marker.x[] = 0.;
      }
    }

    foreach_edge_b()
      s.x[] = snew.x[];

    // We should take care of the vertex vector filed on the resolution boundary manually?
    boundary ((scalar *) {s});

    // Output the information of cell face with odd number of markers, for debugging
    // we need noauto, the boundary() produces strange result
    foreach_face(z, x, y, noauto) {
      int ii = 0;
      ii += with_marker.x[] + with_marker.x[0,1] + with_marker.y[] + with_marker.y[1];

      if (ii % 2 != 0) {
        printf ("Illega number of markers (%d) on face:\n", ii);
        printf ("PID: %d x:%g, y:%g, z:%g, ii:%g, jj:%g, kk:%g\n", pid(), \
          x, y, z, x/Delta, y/Delta, z/Delta);
        printf ("color vertex: %g x:%g, y:%g, z:%g, \n", color_pha[], color_pha[1], \
          color_pha[0,1], color_pha[1,1]);
        printf ("with_markers: %g x:%g, y:%g, z:%g, \n", with_marker.x[], \
          with_marker.x[0,1], with_marker.y[], with_marker.y[1]);
      }
    }

    foreach_edge() {
      int nm = (int) with_marker.x[];
      if (nm > 1 || nm < 0) {
        printf ("Illega number of markers (%d) on edge: \n", nm);
        printf ("x:%g, y:%g, z:%g, ii:%g, jj:%g, kk:%g\n", x, y, z, x/Delta, y/Delta, z/Delta);
      }
    }
  }

/**
## Multi-dimensional EBIT marker advection */

void ebit_advection (vector u, int i) {
  void (* sweep[dimension]) ();
  int d = 0;
  foreach_dimension()
    sweep[d++] = advect_x;
  for (d = 0; d < dimension; d++) {
  #ifdef SWAP
    sweep[(i + d) % dimension] (u, i);
  #else
    sweep[(d) % dimension] (u, i);
  #endif
  }

  #ifdef ADAPT
  // marks the interfacial cell
  scalar with_intf[];

  foreach() {
    int cv = 0;
    for (int kk = 0; kk < 2; kk++)
      for (int jj = 0; jj < 2; jj++)
        for (int ii = 0; ii < 2; ii++)
          cv += (int) color_pha[ii,jj,kk];

    with_intf[] = (cv == 8 || cv == 0) ? 0. : 1.;
  }

  /* We only need to refine two adjacent layers of the interfacial cell
  to the maximum level, since vertex vertex is used to store with_markers,
  which is different to the 2D version. */
  set_mask (with_intf, 2);
  #endif

  semu2vof();
  
  #if _MYOUTPUT
  debug_log (i);
  #endif

  #ifdef OUT_GNUPLOT
  if (i % DIT == 0)
    output_intf (i);
  #endif
}

/** EBIT marker advection, we use the name "vof" to achieve consistency with VOF solver in Basilisk */

event vof (i++) {
  for (vector u in u_ebit)
    ebit_advection (u, i);
}

#if ADAPT
event adapt (i++) {
  set_markers(); // we need this to set the correct value for with_markers[] after adapt_wavelet

  #ifdef OUT_GNUPLOT
  // output interface again after AMR
  if (i % DIT == 0)
    output_intf (i);
  #endif
}
#endif