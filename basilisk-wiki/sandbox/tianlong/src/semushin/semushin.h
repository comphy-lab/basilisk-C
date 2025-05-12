/**
# Edge Based Interface Tracking solver

We wish to advect an interface identified by marker points that are constrained
to move only along the grid lines. We refer to it as EBIT method and
the following assumptions make it a Semushin's method:
The underlying grid is a 2D square grid, there are at most four markers per cell and at most one markers per edge.

*/
static scalar * _interface=NULL;

#if USE_MY_SOLID
extern scalar is_solid;
extern vertex scalar is_solid_vertex;
extern face vector is_solid_face;
void imposeBCforSemushinAdvection(face vector st, face vector sn);
#endif

// To avoid the advection of "f"
event vof(i++) {
  interfaces = _interface;
}

#include "fitting.h"
#include "vof.h"
#include "ebit_utils.h"
#include "semushin_tree.h"

// To avoid the advection of "f"
event vof(i++) {
  _interface = interfaces;
  interfaces = NULL;
}

#define bilinear(a1, a2, a3, a4, x, y) ((1. - y) * (1. - x) * a1 + (1. - y) * x * a2 \
  + y * x * a3 + y * (1. - x) * a4)

// make sure that CFL <= 0.5 when NTS == 2
#ifdef PRECOR
  #define NTS 2
#else
  #define NTS 1
#endif

face vector s[], snew[], s_tmp[], with_marker[], ss_tmp[];

/** color_pha[]: corner color vertex, color_pha_cen[]: central color vertex,
  config_dict[]: topology configuration dictionary
  there is something wrong with the automatic global vertex scalar, which causes segmentation fault.*/
vertex scalar color_pha, color_pha_new;
scalar color_pha_cen[], config_dict[];
double tTime = 0.;

static inline void restriction_empty(Point point, scalar s)
{  
  //do nothing
}

/** Function for updating the central color point */
foreach_dimension()
  static void update_color_cen_x(){
    foreach(){
#if USE_MY_SOLID
      if((int)is_solid[] == 1)
        continue;
#endif
      int ind, ind_diff, issame, conf_old;
      double cen_old = color_pha_cen[];
      ind = (int) (color_pha_new[] + color_pha_new[1, 0] + color_pha_new[1, 1] + color_pha_new[0, 1]);
      ind_diff = (int) (fabs(color_pha_new[] - color_pha[]) + fabs(color_pha_new[1, 0] - color_pha[1, 0]) \
        + fabs(color_pha_new[1, 1] - color_pha[1, 1]) + fabs(color_pha_new[0, 1] - color_pha[0, 1]));
      issame = ((int) color_pha_new[] == (int) color_pha_new[1, 0])\
        || ((int) color_pha_new[] == (int) color_pha_new[0, 1]);
      conf_old = (int) config_dict[];

      if (ind > 2){
        color_pha_cen[] = 1.;
      }
      else if (ind < 2){
        color_pha_cen[] = 0.;
      }
      else if (ind == 2){
        if (issame){
          color_pha_cen[] = -1;
        }
        else if (ind_diff == 3){
          int is, conf_nei, ie1, ie2, ie3, iy1, iy2, ischange;
          // for new conf = 7 or 8, we should make sure that there is no cell with
          // conf = 7, 8 in the initial state, otherwise the code will fail
          // see the slide for all possible changes of the central color vertex

          // from 3, 4, 5, 6 to 7, 8
          iy1 = (conf_old == 5 || conf_old == 6) ? 1 : 0;
          iy2 = 1 - iy1;
          if (conf_old == 3 || conf_old == 6){
            is = -1;
            ie2 = 0;
            // It's better to use a criteria based on the change of color vertex.
            ischange = (s_tmp.y[0, iy1] > 0.) && (s_tmp.y[is, iy2] > 0.);
          }
          else{
            is = 1;
            ie2 = 2;
            ischange = (s_tmp.y[0, iy1] < 0.) && (s_tmp.y[is, iy2] < 0.);
          }
          ie1 = get_end(conf_old, ie2);

          conf_nei = config_dict[is, 0];
          ie3 = 2 - ie2;
          ie3 = get_end(conf_nei, ie3);

          if (ie3 + ie1 == 4 && ischange){
            // two end points of the segments should on the oppsite edges,
            // and also, both of them should move cross the vertical edge on the same direction.
            color_pha_cen[] = 1. - cen_old;
          }
          printf("complex rule for ind_diff %d: conf_old: %d, color_pha_cen: %g\n", ind_diff, conf_old, color_pha_cen[]);
        }
        else if (ind_diff == 2 && conf_old == 1){
          // from 1 to 7, 8
          int ischange=0, low_diff;
          low_diff = (int) (fabs(color_pha_new[] - color_pha[]) + fabs(color_pha_new[1, 0] - color_pha[1, 0]));
          ischange = (s_tmp.y[0, 0] < 0. && low_diff != 0) || (s_tmp.y[0, 1] > 0. && low_diff == 0);
          color_pha_cen[] = ischange ? color_pha[] : color_pha[1, 0];
          printf("complex rule for ind_diff %d: conf_old: %d, color_pha_cen: %g\n", ind_diff, conf_old, color_pha_cen[]);
        }
        else if (ind_diff == 2 && conf_old == 2){
          // from 2 to 7, 8
          int is, ischange=0, left_diff, conf_nei;
          double color1, color2;
          left_diff = (int) (fabs(color_pha_new[] - color_pha[]) + fabs(color_pha_new[0, 1] - color_pha[0, 1]));
          color1 = color_pha_new[0, 0];
          color2 = color_pha_new[1, 0];

          is = (left_diff == 0) ? 1 : -1;
          conf_nei = config_dict[is, 0];
          color_pha_cen[] = (conf_nei == 8) ? color2 : color1;
          printf("complex rule for ind_diff %d: conf_old: %d, color_pha_cen: %g\n", ind_diff, conf_old, color_pha_cen[]);
        }
        else if (ind_diff == 4){
          // from 7, 8 to 7, 8
          int is, conf_nei, ie1, ie2, ie3, iy1, iy2, ischange;
          // right neighboor
          iy1 = (conf_old == 7) ? 0 : 1;
          iy2 = 1 - iy1;

          is = 1;
          ie2 = 2;
          ischange = (s_tmp.y[0, iy1] < 0.) && (s_tmp.y[is, iy2] < 0.);
          ie1 = get_end(conf_old, ie2);

          conf_nei = config_dict[is, 0];
          ie3 = 2 - ie2;
          ie3 = get_end(conf_nei, ie3);

          if (ie3 + ie1 == 4 && ischange){
            color_pha_cen[] = 1. - cen_old;
          }

          // left neighboor
          iy1 = (conf_old == 8) ? 0 : 1;
          iy2 = 1 - iy1;

          is = -1;
          ie2 = 0;
          ischange = (s_tmp.y[0, iy1] > 0.) && (s_tmp.y[is, iy2] > 0.);
          ie1 = get_end(conf_old, ie2);

          conf_nei = config_dict[is, 0];
          ie3 = 2 - ie2;
          ie3 = get_end(conf_nei, ie3);

          if (ie3 + ie1 == 4 && ischange){
            color_pha_cen[] = 1. - cen_old;
          }
          printf("complex rule for ind_diff %d: conf_old: %d, color_pha_cen: %g\n", ind_diff, conf_old, color_pha_cen[]);
        }
      }
    }
  }

/** Update the topology dictionary based on the color vertex.
 Be careful inside the foreach_dimension iterator.
 We need different results for different directions, so we use foreach_dimension*/
foreach_dimension()
  static double get_cdict_x(Point point){
#if USE_MY_SOLID
    if((int)is_solid[] == 1)
      return 0;
#endif
    int ii=0, ind, ind_cen, issame, col_ver[4];
    double conf;
    conf = 0.; // Without interface

    // Start calculating the connections of markers based on the color vertex
    // counter-clockeise for get_cdict_x, clockwise for get_cdict_y
    col_ver[0] = (int) color_pha[];
    col_ver[2] = (int) color_pha[1, 1];
    col_ver[1] = (int) color_pha[1, 0];
    col_ver[3] = (int) color_pha[0, 1];

    ind = col_ver[0] + col_ver[1] + col_ver[2] + col_ver[3];
    ind_cen = (int) (color_pha_cen[]);
    issame = (col_ver[0] == col_ver[1]);

    if (ind_cen == -1){ // Connect the opposite edges, one interface
      conf = (issame) ? 2. : 1.;
    }
    else if (ind == 3 || ind == 1){ // Connect the consecutive edges, one interface
      ii = 0;
      for (int iv=0; iv<4; iv++){
        if (col_ver[iv] != ind_cen){
          ii = iv;
          break;
        }
      }
      conf = 3. + ii;
    }
    else if (ind == 2){ // Connect the consecutive edges, two interfaces
      conf = (col_ver[0] == ind_cen) ? 7. : 8.;
    }
    return conf;
  }

foreach_dimension()
  static void update_dict_x(){
    foreach(){
      config_dict[] = get_cdict_x(point);
    }
    boundary({config_dict});
  }


/** Injection and restriction functions for markers and color vertex */

foreach_dimension()
  static void myrefine_face_x (Point point, scalar s)
  {
    /** In the current version, the 3 by 3 stencil of the interfacial cell are refined
      to maximum level after initialization, so refinement does not take place in the
      interfacial cell during the advection.
      For the refinement of an interfacial cell, I only tested this function for a simple case (circle).
      We may need this in the future.*/
    vector v = s.v;
    // we can't use the config_dict[] directly, so we recalcuate the dictionary
    int conf, ii=0, iip, ind_markers[4]={-1, -1, -1, -1};
    double xy_edge[4][2];

    conf = get_cdict_x(point);

    // if (conf != 0){
    //   printf("x: %g y: %g, conf: %d\n", x, y, conf);
    // }
    xy_edge[0][0] = 0.;
    xy_edge[0][1] = v.x[0, 0];

    xy_edge[1][0] = v.y[0, 0];
    xy_edge[1][1] = 0.;

    xy_edge[2][0] = 1.;
    xy_edge[2][1] = v.x[1, 0];

    xy_edge[3][0] = v.y[0, 1];
    xy_edge[3][1] = 1.;

    for (int iind=0; iind<4; iind++){
      ind_markers[iind] = dict_to_edge[conf][iind];
    }
    //
    if (!is_refined(neighbor(-1)) &&
      (is_local(cell) || is_local(neighbor(-1)))) {
        fine(v.x, 0, 1) = 0.;
        fine(v.x, 0, 0) = 1.; // This makes sure that we get a consistent representation for corner case

        if (v.x[0] >= 0.5){
          fine(v.x, 0, 1) = (v.x[0] - 0.5) * 2.;
        }
        else{
          fine(v.x, 0, 0) = v.x[0] * 2.;
        }
      }  
      
      if (!is_refined(neighbor(1)) && neighbor(1).neighbors &&
        (is_local(cell) || is_local(neighbor(1)))) {
          fine(v.x, 2, 1) = 0.;
          fine(v.x, 2, 0) = 1.;

          if (v.x[1] >= 0.5){
            fine(v.x, 2, 1) = (v.x[1] - 0.5) * 2.;
          }
          else {
            fine(v.x, 2, 0) = v.x[1] * 2.;
          }
        }
        
        if (is_local(cell)){     //Mid line
          double x1, y1, x2, y2, y0;
          fine(v.x, 1, 0) = 1.;
          fine(v.x, 1, 1) = 0.;
          for (int iv=0; iv<4; iv+=2){
            if (ind_markers[iv] != -1){
              int ie1=ind_markers[iv], ie2=ind_markers[iv + 1];
              x1 = xy_edge[ie1][0];
              y1 = xy_edge[ie1][1];

              x2 = xy_edge[ie2][0];
              y2 = xy_edge[ie2][1];

              if (x1 == x2 && x1 == 0.5)
                y0 = 0.5;
              else
                y0 = my_intersect(x1, y1, x2, y2, 0.5);
              if (y0 >= 0.5){
                fine(v.x, 1, 1) = (y0 - 0.5) * 2.;
              }
              else if (y0 >= 0){
                fine(v.x, 1, 0) = y0 * 2.;
              }
            }
          }
          
        }
  }

/** For color vertex*/
static void myrefine_vertex (Point point, scalar ss)
{ 
  // we can't use the config_dict, so we recalcuate the dictionary
  int conf, ii=0, iip, ind, ind_cen, issame;
  double xy_edge[4][2];

  conf = get_cdict_x(point);

  xy_edge[0][0] = 0.;
  xy_edge[0][1] = s.x[0, 0];

  xy_edge[1][0] = s.y[0, 0];
  xy_edge[1][1] = 0.;

  xy_edge[2][0] = 1.;
  xy_edge[2][1] = s.x[1, 0];

  xy_edge[3][0] = s.y[0, 1];
  xy_edge[3][1] = 1.;

  // consistent treatment for case where interface passes through the central color vertex
  double x1, y1, x2, y2, y0, xym;
  ind_cen = (int) color_pha_cen[];
  if (ind_cen != -1){
    fine(ss, 1, 1) = ind_cen;
  }
  else{
    int ie1, ie2;
    if (conf == 1){
      ie1 = 1;
      ie2 = 3;
    }
    else if (conf == 2){
      ie1 = 0;
      ie2 = 2;
    }
    x1 = xy_edge[ie1][0];
    y1 = xy_edge[ie1][1];

    x2 = xy_edge[ie2][0];
    y2 = xy_edge[ie2][1];
    xym = (conf == 1) ? 0.5 * (x1 + x2) : 0.5 * (y1 + y2);

    fine(ss, 1, 1) = (xym >= 0.5) ? ss[] : 1. - ss[];
  }

  // should check the meaning of allocated_child function
  int have_child=1;
  for (int i = 0; i <= 1; i++) {
    for (int j = 0; j <= 1; j++)
      if (allocated_child(2 * i, 2 * j))
        fine(ss, 2 * i, 2 * j) = ss[i, j];
      else{
        // printf("no_child: %d, i: %d, j: %d\n", have_child, i, j);
        have_child = 0;
      }
    
      foreach_dimension()
        if (neighbor(i).neighbors) {
          fine(ss, 2 * i, 1) = (s.x[i, 0] >= 0.5) ? ss[i, 0] : ss[i, 1];
        }
  }

  // refinement for the central color vertex
  // printf("x:%g, y:%g, cen:%g %d\n", x, y, ss[], is_leaf(cell));
  // fflush(stdout);
  if (have_child){
    // printf("havechild, x:%g, y:%g, cen:%g %d\n", x, y, ss[], is_leaf(cell));
    // fflush(stdout);
    ii = color_pha_cen[];
    for (int i = 0; i <= 1; i++) {
      for (int j = 0; j <= 1; j++){
        fine(color_pha_cen, i, j) = ii;
        int ind, c1, c2, c3, c4;
        c1 = (int) fine(ss, i, j);
        c2 = (int) fine(ss, i + 1, j);
        c3 = (int) fine(ss, i + 1, j + 1);
        c4 = (int) fine(ss, i, j + 1);
        ind = c1 + c2 + c3 + c4;
        issame = (c1 == c2) || (c1 == c4);
        if (ind > 2){
          fine(color_pha_cen, i, j) = 1.;
        }
        else if (ind < 2){
          fine(color_pha_cen, i, j) = 0.;
        }
        else if (ind == 2 && issame){
          fine(color_pha_cen, i, j) = -1;
        }
      }
    }
  }
  
}


void output_facets_semushin(char *file){
  char name[strlen(file) + 2];
  strcpy (name, file);
  update_dict_x();

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

  foreach(serial){
#if USE_MY_SOLID
    if((int)is_solid[] == 1)
      continue;
#endif
    int ii=0, iip, ind_markers[4]={-1, -1, -1, -1};
    double xy_edge[4][2];

    xy_edge[0][0] = x - Delta / 2.;
    xy_edge[0][1] = y - Delta / 2. + s.x[0, 0] * Delta;

    xy_edge[1][0] = x - Delta / 2. + s.y[0, 0] * Delta;
    xy_edge[1][1] = y - Delta / 2.;

    xy_edge[2][0] = x + Delta / 2.;
    xy_edge[2][1] = y - Delta / 2. + s.x[1, 0] * Delta;

    xy_edge[3][0] = x - Delta / 2. + s.y[0, 1] * Delta;
    xy_edge[3][1] = y + Delta / 2.;

    /** Start calculating the connections of markers basing on the 
      color vertex (central point)*/
    int conf=(int) config_dict[];

    for (int iind=0; iind<4; iind++){
      ind_markers[iind] = dict_to_edge[conf][iind];
    }

    for (int iv=0; iv<4; iv+=2){
      if (ind_markers[iv] != -1){
        double xx[2], yy[2];
        int ie1=ind_markers[iv], ie2=ind_markers[iv + 1];
        xx[0] = xy_edge[ie1][0];
        yy[0] = xy_edge[ie1][1];

        xx[1] = xy_edge[ie2][0];
        yy[1] = xy_edge[ie2][1];

        fprintf (fp1, "%g %g\n%g %g\n\n", xx[0], yy[0], xx[1], yy[1]);
      }
    }
  }

  fflush(fp1);
  fclose(fp1);
  #if _MPI
  if (pid() + 1 < npe())
    MPI_Send(&sign_mpi, 1, MPI_INT, pid() + 1, 0, MPI_COMM_WORLD);
  #endif
}


static void output_intf(int i) {
  char out[100], out_test[100], out_color[100];
  const char testName[50] = TEST;

  snprintf(out, sizeof(out), "%s_%d_%d.dat", testName, N, i); 
  snprintf(out_color, sizeof(out_color), "%s_%d_%d.dat", "output_debug/color_result", N, i);
  snprintf(out_test, sizeof(out_test), "%s_%d.dat", "output_debug/grid", i);

  output_facets_semushin(out);
  output_color_vertex(file=out_color, cv=color_pha, cv_center=color_pha_cen);
  output_mesh(out_test);
}


/** In `full` cells (volume fraction $f=1$) there is no interface and then the
   Semushin fraction is zero. To obtain the Semushin interface position we have
   to `rotate` the face fractions $s\_tmp$, since they are computed with respect
   to the normal, while in Semushin method they always refer to fixed vertexes
   (e.g. left and bottom). "phi" is levelset function used*/
static void init_markers(vertex scalar phi){
  // workaround, segmentation fault for global "vertex scalar color_pha[]"
  color_pha = new vertex scalar;
  color_pha_new = new vertex scalar;

  boundary ({phi});
  fractions (phi, f, s_tmp);
  boundary ((scalar *){f, s_tmp});

  // Initialize the color vertex
  foreach_vertex(){
    if (phi[] >= 0.){ // remember the equal sign
      color_pha[] = 1.; // Reference phase
    }
    else{
      color_pha[] = 0.;
    }
#if USE_MY_SOLID
    color_pha[] *= (1.0 - is_solid_vertex[]);
#endif
    color_pha_new[] = color_pha[];
  }

  // For the color_pha_cen[] is undefined (-1) when only two consecutive color vertices
  // are the same color. When more then two vertices are the same color, we change the
  // color of central color point.
  foreach(){
    color_pha_cen[] = -1.;
  }

  update_color_cen_x();

  boundary({color_pha_cen});

  foreach_face(){
    s.x[] = 0.;
    with_marker.x[] = 0.;
  }

  // Different values of color vertex at the two end of edge means there is a marker on the edge.
  // This method can deal with the corner case (interface pass through the cell vertex) correctly.
  // with_marker : number of markers on each edge
  foreach_face(){
#if USE_MY_SOLID
    if ((int)is_solid_face.x[] == 1)
    {
      s.x[] = 0.;
      with_marker.x[] = 0.;
      continue;
    }
#endif
    int with_face=fabs(color_pha[0, 0] - color_pha[0, 1]) > machine_zero;
    if (with_face){
      if (color_pha[] > machine_zero)
        s.x[] = s_tmp.x[];
      else
        s.x[] = 1. - s_tmp.x[];

      with_marker.x[] = 1.;
    }
  }

  // The boundary condtions for the scalar field used for the Semushin markers are set. It's not use now
  s_tmp.t[left] = 0.;
  s_tmp.t[right] = 0.;
  s_tmp.t[top] = 0.;
  s_tmp.t[bottom] = 0.;

  snew.t[left] = -min(0., s_tmp.y[0, 0]);
  snew.t[right] = max(0., s_tmp.y[0, 0]);

  snew.t[bottom] = -min(0., s_tmp.x[0, 0]);
  snew.t[top] = max(0., s_tmp.x[0, 0]);

  snew.t.depends = list_add(snew.t.depends, s_tmp.x);
  snew.t.depends = list_add(snew.t.depends, s_tmp.y);

  config_dict[left] = 0.;
  config_dict[right] = 0.;
  config_dict[top] = 0.;
  config_dict[bottom] = 0.;

  #if TREE
  // must set both prolongation, refine, coarsen, restriction manually
  // printf("Tree test\n");
  foreach_dimension(){
    s.x.restriction = no_restriction;
    s.x.coarsen = no_restriction;
    s.x.prolongation = myrefine_face_x;
  }
  s.x.refine = refine_face;
  s.x.restriction = myrestriction_face;
  s.x.coarsen = myrestriction_face;

  color_pha.restriction = restriction_vertex;
  color_pha.coarsen = restriction_vertex;
  color_pha.prolongation = myrefine_vertex;
  color_pha.refine = myrefine_vertex;

  // central color vertex is refine in color_pha.refine
  // we still need restriction
  color_pha_cen.prolongation = refine_injection;
  color_pha_cen.refine = restriction_empty; //for solid

  color_pha_cen.restriction = restriction_conf;
  color_pha_cen.coarsen = restriction_conf;

  // for color_pha_cen_new
  // 1. we need restriction to set vertex shared by two different levels
  // to the correct value (value at the finer level)
  // 2. we need a refine function to set vertex shared by two different levels
  // to the correct value (when the vertex is also shared by two different processors)
  // but why ???

  // we need the refine function to make MPI transfering the data
  // between the ghost cells. Do not use no_restriction for
  // the refine function
  config_dict.prolongation = no_restriction;
  config_dict.restriction = restriction_conf;
  // config_dict.refine = no_restriction;
  config_dict.coarsen = no_restriction;
  #endif

#if MYDEBUG
  // output the initial interface for test 
  output_intf(-1);
#endif
}


void set_markers(){
  // determine the number of markers on edge based on color vertex, this is a robust method for corner case,
  // it can be generalized to double-Semushin easily in the furture.
    foreach_face(){
#if USE_MY_SOLID
    if ((int)is_solid_face.x[] == 1)
    {
      with_marker.x[] = 0.;
      s.x[] = 0.;
      continue;
    }
#endif
      int with_face=fabs(color_pha[0, 0] - color_pha[0, 1]) > machine_zero;
      if (with_face){
        with_marker.x[] = 1.;
      }
      else{ // without interface or with two interface
        with_marker.x[] = 0.;
        s.x[] = 0.;
      }
    }
}

/**
The representation of the interface with Semushin markers can be used to create
the associated VOF field. The connection of makers and the region of reference phase are 
idenfied by the "color vertex". 
The VOF field is here used to easily compute the area of the reference phase. */

scalar f1[];
double area, area0;
void semu2vof() {
  update_dict_x();
  foreach(){
#if USE_MY_SOLID
    if ((int)is_solid[] == 1)
    {
      f1[] = 0.;
      continue;
    }
#endif
    int cv1=(int) (color_pha[] + color_pha[1, 0] + color_pha[1, 1] +color_pha[0, 1]);
    if (cv1 == 4){
      f1[] = 1.;
    }
    else{
      f1[] = 0.;
    }
  }

  foreach() {
#if USE_MY_SOLID
    if ((int)is_solid[] == 1)
    {
      continue;
    }
#endif
    int cv1, cv2, cv3, cv4, cv, conf;
    double ss;
    cv1 = (int) (color_pha[]);
    cv2 = (int) (color_pha[1, 0]);
    cv3 = (int) (color_pha[1, 1]);
    cv4 = (int) (color_pha[0, 1]);
    conf = (int) (config_dict[]);

    if (conf != 0){ // interfacial cell, cv is used for identifying the reference phase
      if (conf == 1){
        ss = (s.y[0, 0] + s.y[0, 1]) / 2.;
        cv = cv1;
      }
      else if (conf == 2){
        ss = (s.x[0, 0] + s.x[1, 0]) / 2.;
        cv = cv1;
      }
      else if (conf == 3){
        ss = (s.x[0, 0] * s.y[0, 0]) / 2.;
        cv = cv1;
      }
      else if (conf == 4){
        ss = ((1. - s.y[0, 0]) * s.x[1, 0]) / 2.;
        cv = cv2;
      }
      else if (conf == 5){
        ss = ((1. - s.x[1, 0]) * (1. - s.y[0, 1])) / 2.;
        cv = cv3;
      }
      else if (conf == 6){
        ss = ((1. - s.x[0, 0]) * s.y[0, 1]) / 2.;
        cv = cv4;
      }
      else if (conf == 7){ // Connect the consecutive edges, two interfaces
        ss = ((1. - s.x[0, 0]) * s.y[0, 1]) / 2. \
          + ((1. - s.y[0, 0]) * s.x[1, 0]) / 2.;
          cv = cv2;
      }
      else if (conf == 8){
        ss = (s.x[0, 0] * s.y[0, 0]) / 2. \
          + ((1. - s.x[1, 0]) * (1. - s.y[0, 1])) / 2.;
        cv = cv1;
      }

      // correct the volume fractions
      #ifdef CIRCLE_TEST
      double xy_edge[4][2];
      int ind_markers[4] = {-1, -1, -1, -1};

      xy_edge[0][0] = 0.;
      xy_edge[0][1] = s.x[];

      xy_edge[1][0] = s.y[];
      xy_edge[1][1] = 0.;

      xy_edge[2][0] = 1.;
      xy_edge[2][1] = s.x[1, 0];

      xy_edge[3][0] = s.y[0, 1];
      xy_edge[3][1] = 1.;

      for (int iind=0; iind<4; iind++){
        ind_markers[iind] = dict_to_edge[conf][iind];
      }
      // we need this beacause of the direction of segment
      if (conf == 1)
        swap(double, ind_markers[0], ind_markers[1]);

      for (int iv=0; iv<4; iv+=2){
        if (ind_markers[iv] != -1){
          double x1, y1, x2, y2, x3, y3, xrc, yrc, rc;
          double xp[4], yp[4];
          int np=2;
          int ie1=ind_markers[iv], ie2=ind_markers[iv + 1];

          x1 = xy_edge[ie1][0];
          y1 = xy_edge[ie1][1];

          x2 = xy_edge[ie2][0];
          y2 = xy_edge[ie2][1];

          xp[0] = x1; xp[1] = x2;
          yp[0] = y1; yp[1] = y2;
          // find out the four markers used for circle fit
          // not the same circle fit in the advect_x
          // The circle pass the first and second points, the last one or two point(s)
          // are used to define the loss function
          int iee[2]={ie1, ie2}, ipx, ipy, conn, ipps, ippe, nsec=0;
          double yave=0.,t1, t2, dtheta, area, sig;;
          for (int ipe=0; ipe<2; ipe++){
            if (iee[ipe] == 0){
              ipx = -1; ipy = 0;
            }
            else if (iee[ipe] == 1){
              ipx = 0; ipy = -1;
            }
            else if (iee[ipe] == 2){
              ipx = 1; ipy = 0;
            }
            else{
              ipx = 0; ipy = 1;
            }
            conn = (int) config_dict[ipx, ipy];
            if (conn == 0)
              continue;
            ipps = 2 * (iee[ipe] % 2 + 1) - iee[ipe];
            ippe = get_end(conn, ipps);

            if (ippe == 0){
              x3 = ipx;
              y3 = s.x[ipx, ipy] + ipy;
            }
            else if (ippe == 1){
              x3 = s.y[ipx, ipy] + ipx;
              y3 = ipy;
            }
            else if (ippe == 2){
              x3 = ipx + 1.;
              y3 = s.x[ipx + 1, ipy] + ipy;
            }
            else{
              x3 = s.y[ipx, ipy + 1] + ipx;
              y3 = 1. + ipy;
            }

            xp[np] = x3; yp[np] = y3;
            np++;
          }

          fit_circle(xp, yp, np, &xrc, &yrc, &rc);

          t1 = atan2(y1 - yrc, x1 - xrc);
          t2 = atan2(y2 - yrc, x2 - xrc);
          dtheta = fabs(t2 - t1);
          dtheta = (dtheta > pi) ? 2. * pi - dtheta : dtheta;
          sig = (((x1 - xrc) * (y2 - yrc) - (x2 - xrc) * (y1 - yrc)) > 0.) ? -1.: 1.;;
          if (rc > 0.)
            yave = 0.5 * sig * rc * rc * (dtheta - sin(dtheta));

          ss += yave;
        }
      }
      #endif
      f1[] = (cv == 1) ? ss : 1. - ss;
    }
  }

  foreach()
    f[] = f1[];
  boundary({f});
  area = 0.;
  foreach(reduction(+:area)) {
#if USE_MY_SOLID
    if ((int)is_solid[] == 1)
    {
      continue;
    }
#endif
    area += f[]*sq(Delta);
  }
}

/**
  This event perform the 1-D advection of the Semushin marker points. */
foreach_dimension()
  void advect_x(vector u) {
    int idim = (int) iadv.x + 1; // idim: for identifing the advection direction

    double tmp = 0;
    set_markers();

    update_dict_x();
      
    foreach_face() {
      snew.x[] = 0.;
      s_tmp.x[] = 0.;
      ss_tmp.x[] = 0.;
    }

    // it's needed in AMG
    foreach_vertex() {
      color_pha_new[] = color_pha[];
    }

    boundary ((scalar *){s_tmp});

    foreach_face(y) {
#if USE_MY_SOLID
      if ((int)is_solid_face.y[] == 1)
        continue;
#endif
       // This performs the horizontal advection
      double um = 0., u1, u2, u3, u4, xx, yy, sy;
      if ((int) with_marker.y[] > 0) {
        sy = s.y[];
        for (int its = 0; its < NTS; its++) {
          #ifdef BILINEAR
          if (sy >= 0.5) {
            u1 = u.x[0, -1];
            u2 = u.x[1, -1];
            u3 = u.x[1, 0];
            u4 = u.x[0, 0];
            xx = sy - 0.5;
          }
          else {
            u1 = u.x[-1, -1];
            u2 = u.x[0, -1];
            u3 = u.x[0, 0];
            u4 = u.x[-1, 0];
            xx = sy + 0.5;
          }
          yy = 0.5;

          if (idim == 2) {
            um += bilinear(u1, u4, u3, u2, yy, xx);
          }
          else {
            um += bilinear(u1, u2, u3, u4, xx, yy);
          }
          #else
          um += (u.x[0, 0] + u.x[0, -1]) / 2.;
          #endif
          sy += um * dt / Delta;
        }
        um /= NTS;

        ss_tmp.y[] = s.y[] + um * dt / Delta;
        snew.y[] = s.y[] + um * dt / Delta;

        #if AXI
        if (idim==2){
          double eps_axi=2.e-2;
          yy = y / Delta + (ss_tmp.y[] - 0.5);
          if (yy >= 0. && yy < eps_axi){
            ss_tmp.y[] = -yy;
            snew.y[] = -yy;
          }

        }
        #endif
      }
    }
    // horizontal displacement for the unaligned markers
    foreach_face(x) {
#if USE_MY_SOLID
      if ((int)is_solid_face.x[] == 1)
        continue;
#endif
      double um = 0., u1, u2, u3, u4, xx, yy;
      xx = 0.5;
      for (int its = 0; its < NTS; its++) {
        #ifdef BILINEAR
        if (s.x[] >= 0.5) {
          u1 = u.x[-1, 0];
          u2 = u.x[0, 0];
          u3 = u.x[0, 1];
          u4 = u.x[-1, 1];
          yy = s.x[] - 0.5;
        }
        else {
          u1 = u.x[-1, -1];
          u2 = u.x[0, -1];
          u3 = u.x[0, 0];
          u4 = u.x[-1, 0];
          yy = s.x[] + 0.5;
        }

        if (idim == 2) {
          um += bilinear(u1, u4, u3, u2, yy, xx);
        }
        else{
          um += bilinear(u1, u2, u3, u4, xx, yy);
        }
        #else
        um += (u.x[0, 0] + u.x[-1, 0]) / 2.;
        #endif

        xx += um * dt / Delta;
      }
      um /= NTS;

      ss_tmp.x[] =  um * dt / Delta;
    }

    foreach_face(y){
#if USE_MY_SOLID
      if ((int)is_solid_face.y[] == 1)
        continue;
#endif
      double tmpy;
      if (snew.y[] > 1.){  // The point is crossing a vertical line with $u>0$.
        tmpy = snew.y[] - 1.;
        snew.y[] = 0.;
        s_tmp.y[] = tmpy;
        with_marker.y[] -= 1;
      }
      else if (snew.y[] < 0.){  // The point is crossing a vertical line with $u<0$.
        tmpy = 1. + snew.y[];
        snew.y[] = 0.;
        s_tmp.y[] = -tmpy;
        with_marker.y[] -= 1;
      }
    }
    
    boundary ((scalar *){s_tmp, snew, ss_tmp});
#if USE_MY_SOLID
    imposeBCforSemushinAdvection(s_tmp, snew);
#endif
    // When there are two markers on the same edge, snew for the old one, s_tmp for the new one
    // For the simple Semushin method, We will remove these two markers on the
    // final stage based on the color vertex.
    foreach_face(y){
#if USE_MY_SOLID
      if ((int)is_solid_face.y[] == 1)
        continue;
#endif
      int withi=(int) with_marker.y[];
      if (withi == 0){
        if (exist(s_tmp.y[-1, 0]))
          snew.y[] = s_tmp.y[-1, 0];
        else if (existNeg(s_tmp.y[1, 0])){
          snew.y[] = -s_tmp.y[1, 0];
        }

        if (exist(s_tmp.y[-1, 0]) && existNeg(s_tmp.y[1, 0])) 
          printf("Double markers on edge: x:%g, y:%g, |current:%g, left:%g, right:%g\n",\
           x / Delta, y / Delta, snew.y[], s_tmp.y[-1,0], s_tmp.y[1, 0]);
      }
      else{
        if (exist(s_tmp.y[-1, 0]) || existNeg(s_tmp.y[1, 0])) 
          printf("Double markers on edge: x:%g, y:%g, |current:%g, left:%g, right:%g\n",\
           x / Delta, y / Delta, snew.y[], s_tmp.y[-1,0], s_tmp.y[1, 0]);
      }

      with_marker.y[] += (exist(s_tmp.y[-1, 0]) + existNeg(s_tmp.y[1, 0]));
    }

    // Change the color of vertex when marker moves across the grid line
    foreach_vertex(){
#if USE_MY_SOLID
      if ((int)is_solid_vertex[] == 1)
        continue;
#endif
      if (exist(s_tmp.y[-1, 0]))
        color_pha_new[] = color_pha[-1, 0];
      else if (existNeg(s_tmp.y[0, 0])){
        color_pha_new[] = color_pha[1, 0];
      }
    }

    boundary((scalar *){color_pha_new});

    foreach_face(x){
#if USE_MY_SOLID
      if ((int)is_solid_face.x[] == 1)
        continue;
#endif
      double y3, y2, y1, x3, x2, x1;
      int withi=(int) with_marker.x[];
      with_marker.x[] = 0.; // should be removed

      if (fabs(ss_tmp.x[]) < machine_zero){
        // For marker at the physical boundary (for symmetric boundary)
        if (withi){
          snew.x[] = s.x[];
          with_marker.x[] = 1.; // should be removed
        }
      }
      else{
        for (int is=0; is<2; is++){
          // is=0 for the cell on the right, is=1 for the cell on the left
          // reconnect the markers and calculate the intersection basing on
          // the topology of cells on the both side of the cell face (last time step)
          double xy_edge[4][2], y0, y0c;
          int conf, confl, confr, ind_markers[4] = {-1, -1, -1, -1};

          conf=(int) config_dict[-is, 0];
          if (conf != 0){
            xy_edge[0][0] = ss_tmp.x[-is, 0] - is;
            xy_edge[0][1] = s.x[-is, 0];

            xy_edge[1][0] = ss_tmp.y[-is, 0] - is;
            xy_edge[1][1] = 0.;

            xy_edge[2][0] = 1. + ss_tmp.x[1 - is, 0] - is;
            xy_edge[2][1] = s.x[1 - is, 0];

            xy_edge[3][0] = ss_tmp.y[-is, 1] - is;
            xy_edge[3][1] = 1.;

            for (int iind=0; iind<4; iind++){
              ind_markers[iind] = dict_to_edge[conf][iind];
            }

            for (int iv=0; iv<4; iv+=2){
              if (ind_markers[iv] != -1){
                int ie1=ind_markers[iv], ie2=ind_markers[iv + 1];
                x1 = xy_edge[ie1][0];
                y1 = xy_edge[ie1][1];

                x2 = xy_edge[ie2][0];
                y2 = xy_edge[ie2][1];

                y0 = my_intersect(x1, y1, x2, y2, 0.);
                // y0 == -1 means no intersection
                if (y0 >= 0.){
                  snew.x[] = y0;
                  with_marker.x[] += 1.; // should be removed

                  #ifdef CIRCLE_FIT
                  // find out the four markers used for circle fit
                  // Current version of implementation: the intersection point is the average of two circle fits
                  int iee[2]={ie1, ie2}, ipx, ipy, conn, ipps, ippe, nsec=0;
                  double yave=0.;
                  for (int ipe=0; ipe<2; ipe++){
                    if (iee[ipe] == 0){
                      ipx = -1; ipy = 0;
                    }
                    else if (iee[ipe] == 1){
                      ipx = 0; ipy = -1;
                    }
                    else if (iee[ipe] == 2){
                      ipx = 1; ipy = 0;
                    }
                    else{
                      ipx = 0; ipy = 1;
                    }
                    conn = (int) config_dict[-is + ipx, ipy];
                    if (conn == 0)
                      continue;
                    ipps = 2 * (iee[ipe] % 2 + 1) - iee[ipe];
                    ippe = get_end(conn, ipps);

                    if (ippe == 0){
                      x3 = ss_tmp.x[-is + ipx, ipy] - is + ipx;
                      y3 = s.x[-is + ipx, ipy] + ipy;
                    }
                    else if (ippe == 1){
                      x3 = ss_tmp.y[-is + ipx, ipy] - is + ipx;
                      y3 = ipy;
                    }
                    else if (ippe == 2){
                      x3 = ss_tmp.x[-is + ipx + 1, ipy] - is + ipx + 1.;
                      y3 = s.x[-is + ipx + 1, ipy] + ipy;
                    }
                    else{
                      x3 = ss_tmp.y[-is + ipx, ipy + 1] - is + ipx;
                      y3 = 1. + ipy;
                    }

                    double xrc, yrc, rc, xmin, xmax;
                    get_circle(x1, y1, x2, y2, x3, y3, &xrc, &yrc, &rc);
                    y0c = intersection_circle(xrc, yrc, rc, 0.);

                    // If the intersection point is outside [0., 1.], we abandon the result
                    if (rc > 0. && y0c > 0.){
                      xmin = min(min(x1, x2), x3);
                      xmax = max(max(x1, x2), x3);
                      // we don't use the result of extrapolation
                      if (xmin <= 0. && xmax >= 0.){
                        yave += y0c;
                        nsec++;
                      }
                    }

                    if (nsec > 0){
                      // We revert to straight line fit if all the circle fits fail
                      yave /= nsec;
                      snew.x[] = yave;
                    }
                  }
                  #endif
                }
              }
            }
          }

        }
      }
    }

    // Update the color vertex, don't change the function call order
    update_color_cen_x();
    foreach_vertex(){
      color_pha[] = color_pha_new[];
    }

    boundary({color_pha});
    // determine the number of markers on edge based on color vertex, robust method for corner case,
    // it can be generalized to double-Semushin easily in the furture.
    // without noauto, it trigger the prolongation of color_pha, then onsameside give a floating
    // point exception, should check this later !!!
    foreach_face(){
#if USE_MY_SOLID
      if ((int)is_solid_face.x[] == 1)
      {
        s.x[] = 0.;
        with_marker.x[] = 0.;
        continue;
      }
#endif
      int with_face=fabs(color_pha[0, 0] - color_pha[0, 1]) > machine_zero;
      if (with_face){
        with_marker.x[] = 1.;
      }
      else{ // without interface or with two interface
        snew.x[] = 0.;
        with_marker.x[] = 0.;
      }
    }

    foreach_face()
      s.x[] = snew.x[];
    
    boundary ((scalar *){s}); // It's needed in MPI

    // Output the information of cell with odd number of markers, for debugging
    foreach(noauto) {
#if USE_MY_SOLID
      if((int)is_solid[] == 1)
        continue;
#endif
      int ii = 0;
      ii += with_marker.x[0, 0] + with_marker.x[1, 0] + with_marker.y[0, 0] + with_marker.y[0, 1];
      if (ii % 2 != 0){
        printf("Illega number of markers on cell: %d\n", ii);
        printf("PID: %d idim: %d, x:%g, y:%g, ii:%g, jj:%g\n", pid(), idim, x, y, x / Delta, y / Delta);
        printf("|x00:%g, |x10:%g, |y00:%g, |y01:%g|\n\n", s.x[0,0], s.x[1,0], s.y[0,0], s.y[0,1]);
        printf("|x00:%g, |x10:%g, |y00:%g, |y01:%g|\n\n", with_marker.x[0,0], with_marker.x[1,0],\
         with_marker.y[0,0], with_marker.y[0,1]);
      }
    }
    foreach_face(noauto) {
#if USE_MY_SOLID
      if((int)is_solid_face.x[] == 1)
        continue;
#endif
      int nm = (int) with_marker.x[];
      if (nm > 1 || nm < 0){
        printf("Illega number of markers: %d\n", nm);
        printf("idim: %d, x:%g, y:%g, ii:%g, jj:%g\n", idim, x, y, x / Delta, y / Delta);
      }
    }
  }

//compute the facial fraction respect to fluid 1
void getFaceFraction(face vector ffraction)
{
    foreach_face()
    {
      bool with_face = (int)with_marker.x[] > 0;
      if (with_face)
      {
        ffraction.x[] = color_pha[] > machine_zero ? s.x[] : 1.0 - s.x[];
      }
      else
      {
        ffraction.x[] = color_pha[] > machine_zero ? 1.0 : 0.0;
      }
    }
}

/** EBIT marker advection, we use the name "vof" to make the solver consistent with VOF solver in Basilisk */
event vof (i++) {
#if !SEMUPC
  void (* sweep[dimension]) (vector);
  int d = 0;
#if MYDEBUG
  debug_log(i);
#endif
  foreach_dimension()
    sweep[d++] = advect_x;
  for (d = 0; d < dimension; d++) {
  #ifdef SWAP
    sweep[(i + d) % dimension] (u);
  #else
    sweep[(d) % dimension] (u);
  #endif
  }

  #ifdef ADAPT
  // marks the 3 * 3 stencil of the interfacial cell
  update_dict_x();

  foreach() {
    int with_itf= (config_dict[] > 0.5);
    foreach_neighbor(1){
      with_itf = with_itf || (config_dict[] > 0.5);
    }

    mask_intf[] = with_itf ? 1. : 0.;
  }
  boundary({mask_intf});
  #endif

  semu2vof();
  if (pid() == 0)
    printf("Step: %d, Area is %.9f \n", i, area);

  #ifdef OUT_GNUPLOT
  if (i % DIT == 0)
    output_intf(i);
  #endif

#endif //!SEMUPC
}

#if ADAPT
event adapt (i++) {
  set_markers(); // we need this to set the correct value for with_markers[] after adapt_wavelet
  #ifdef OUT_GNUPLOT
  // output interface again after AMR
  if (i % DIT == 0)
    output_intf(i);
  #endif
}
#endif
