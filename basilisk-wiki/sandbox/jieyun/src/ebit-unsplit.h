/**
# Unsplit advection scheme for the EBIT method
*/

/**
## Color Vertex and dictionary

Function for updating the central color vertex. We implement a simplifed
version for the case with two segments within one cell. We test if the two
cell vertices on the diagnol direction are located on two different sides
of one segment.
*/
static void update_color_cen_unsplit() {
  foreach() {
    int ind, issame, ind_cen_old;
    ind = (int) (color_pha_new[] + color_pha_new[1] + color_pha_new[1,1] + color_pha_new[0,1]);
    issame = ((int) color_pha_new[] == (int) color_pha_new[1])\
      || ((int) color_pha_new[] == (int) color_pha_new[0,1]);
    ind_cen_old = (int) color_pha_cen[];

    if (ind > 2)
      color_pha_cen[] = 1.;
    else if (ind < 2)
      color_pha_cen[] = 0.;
    else if (ind == 2) {
      if (issame)
        color_pha_cen[] = -1;
      else if (ind_cen_old == -1) {
        int nm[2];
        coord xm[8];
        Point pt = point;
        int ns = get_segments (pt, config_dict, st = s_tmp, xm = xm, nm = nm, sn = ss_tmp);
        coord p1, p2, p3 = {0., 0.}, p4 = {1., 1.};

        p1.x = xm[0].x;
        p1.y = xm[0].y;

        p2.x = xm[1].x;
        p2.y = xm[1].y;
        color_pha_cen[] = is_same_side (p1, p2, p3, p4) ? color_pha_new[] : 1. - color_pha_new[];
      }
    }
  }
}

// up: u(t + dt), urk1: u(t + dt/2), urk2: u(t + dt)
extern vector up, urk1, urk2;

/**
## Different time-integration schemes

Use the similar code structure as the
[Runge-Kutta time integrators](/src/runge-kutta.h) in Basilisk.

st, sn are the tangential and normal displacements of a marker, respectively.
*/

double ebit_update (face vector s, face vector with_marker,
  face vector kt, face vector kn,
  face vector st, face vector sn,
  vector uc, double dt, double w, double coef_t = 1.) {

  foreach_face() {
    if (with_marker.x[] > 0.) {
      double umx = 0., umy = 0., xx, yy;
      #if USE_ANA
      // use the analytical velocity for the single vortex test.
      coord xr = {kn.x[]*dt, kt.x[]*dt + (s.x[] - 0.5)*Delta};
      xx = x + xr.x;
      yy = y + xr.y;
      coord um = {sq(sin(xx*pi))*sin(2.*yy*pi), -sin(2.*xx*pi)*sq(sin(yy*pi))};
      umx = um.x*coef_t;
      umy = um.y*coef_t;
      #else
      // compute the marker velocity using the bilinear interpolation
      xx = kn.x[]*dt/Delta - 0.5;
      yy = kt.x[]*dt/Delta + s.x[] - 0.5;

      for (int i = -1; i < 2; i++)
        for (int j = -1; j < 2; j++) {
          double xi = xx - 1.*i, yi = yy - 1.*j;
          xi = DR_L(xi);
          yi = DR_L(yi);

          umx += xi*yi*uc.x[i,j];
          umy += xi*yi*uc.y[i,j];
        }
      #endif

      kn.x[] = umx;
      kt.x[] = umy;

      sn.x[] += w*umx;
      st.x[] += w*umy;
    }
  }
  return w;
}

void ebit_rk (face vector s, face vector with_marker,
  face vector st, face vector sn,
  vector *ul, double dt, int order, double *coef_ts) {
  double w = 0., coef_t;
  face vector kt[], kn[];
  vector uc;

  // ul is the velocity field at different time instants used in the integration
  foreach_face() {
    kt.x[] = kn.x[] = 0.;
    st.x[] = sn.x[] = 0.;
  }

  switch (order) {
  case 1:
    // first-order explicit Euler scheme
    uc = ul[0];
    coef_t = coef_ts[0];
    w += ebit_update (s, with_marker, kt, kn, st, sn, uc, dt, 1., coef_t);
    break;
  case 2:
    // predictor-corrector scheme
    uc = ul[0];
    coef_t = coef_ts[0];
    w += ebit_update (s, with_marker, kt, kn, st, sn, uc, dt, 1., coef_t);
    uc = ul[1];
    coef_t = coef_ts[2];
    w += ebit_update (s, with_marker, kt, kn, st, sn, uc, dt, 1., coef_t);
    break;
  case 4:
    // RK4 scheme
    uc = ul[0];
    coef_t = coef_ts[0];
    w += ebit_update (s, with_marker, kt, kn, st, sn, uc, dt, 1., coef_t);
    uc = ul[1];
    coef_t = coef_ts[1];
    w += ebit_update (s, with_marker, kt, kn, st, sn, uc, dt/2., 2., coef_t);
    w += ebit_update (s, with_marker, kt, kn, st, sn, uc, dt/2., 2., coef_t);
    uc = ul[2];
    coef_t = coef_ts[2];
    w += ebit_update (s, with_marker, kt, kn, st, sn, uc, dt, 1., coef_t);
    break;
  default:
    assert (false); // not implemented
  }

  foreach_face() {
    if (with_marker.x[] > 0.) {
      sn.x[] *= (dt/Delta/w);
      st.x[] = s.x[] + st.x[]*dt/Delta/w;
    }
  }
}

/**
## Unsplit advection scheme */

void advect_unsplit (int ind=0) {
  // ind for debugging purpose
  set_markers();

  update_dict_x();

  foreach_face()
    snew.x[] = 0.;

  // it's needed in AMG
  foreach_vertex()
    color_pha_new[] = color_pha[];

  boundary ({color_pha});

  // advection of markers, explict Euler,  predictor-corrector and RK4 methods
  boundary ((scalar *) {u, up, urk1, urk2});
  if (ebit_order == 1)
    ebit_rk (s, with_marker, s_tmp, ss_tmp, {u}, dt, 1, coef_ts);
  else if (ebit_order == 2)
    ebit_rk (s, with_marker, s_tmp, ss_tmp, {u, up}, dt, 2, coef_ts);
  else
    ebit_rk (s, with_marker, s_tmp, ss_tmp, {u, urk1, urk2}, dt, 4, coef_ts);

  boundary ((scalar *) {s_tmp, ss_tmp});

  /** compute the new position of markers on edge
    search the 2 * 3 (3 * 2) stencil for the x (y) face */

  int idim = -1;
  // maximum curvature ratio, used to select the appropriate fittging result
  // in the sharp corner region or tip region of ligament.
  double kr_max = 10.;
  foreach_dimension() {
    idim++;
    foreach_face(x) {
      for (int i = -1; i < 1; i++)
        for (int j = -1; j < 2; j++) {
          int conf = (int) config_dict[i,j];
          if (conf > 0) {
            Point pt = point;
            if (idim == 0) {
              pt.i += i;
              pt.j += j;
            }
            else {
              pt.i += j;
              pt.j += i;
            }

            int nm[2];
            coord xm[8];
            int ns = get_segments (pt, config_dict, st = s_tmp, xm = xm, nm = nm, sn = ss_tmp);
            for (int iseg = 0; iseg < ns; iseg++) {
              bool with_sl = false, with_circle = false;
              double y0 = HUGE, kp[2] = {0., 0.};
              coord pm[4];
              int np = nm[iseg];

              for (int ip = 0; ip < np; ip++) {
                pm[ip].x = xm[4*iseg + ip].x + i;
                pm[ip].y = xm[4*iseg + ip].y + j;
                if (ip >= 2)
                  kp[ip - 2] = get_kappa_circle (pm[0], pm[1], pm[ip]);
              }
              // curvature ratio of two circle fits
              double kr = max(fabs(kp[0]), fabs(kp[1]))/(min(fabs(kp[0]), fabs(kp[1])) + 1.e-32);

              // refactor this part
              if (pm[0].x*pm[1].x <= 0.) {
                // case with two segment endpoints located on two sides of the cell face
                // intersection point based on the linear fit
                y0 = my_intersect (pm[0].x, pm[0].y, pm[1].x, pm[1].y, 0., -HUGE, HUGE);
                with_sl = y0 >= 0. && y0 <= 1.;

                // circle fit
                int nsec = 0;
                double yave = 0., ymin = min(pm[0].y, pm[1].y),\
                  ymax = max(pm[0].y, pm[1].y), ycyc_int[2] = {HUGE, HUGE};
                double ycyc[4] = {-1., -1., -1, -1}; // for debugging

                for (int ip = 2; ip < np; ip++) {
                  double xrc, yrc, rc, yint[2];
                  get_circle (pm[0].x, pm[0].y, pm[1].x, pm[1].y, pm[ip].x, pm[ip].y, &xrc, &yrc, &rc);
                  if (rc > 0.) {
                    int nint = intersection_circle_two (xrc, yrc, rc, 0., yint);
                    if (nint > 0) {
                      double y0c;
                      int nwithin = 0;
                      // choose the correct intersection point
                      for (int iint = 0; iint < nint; iint++) {
                        if (yint[iint] >= ymin && yint[iint] <= ymax) {
                          nwithin++;
                          y0c = yint[iint];
                        }
                      }

                      // choose the point closer to that obtained with linear fit
                      if (nwithin == 0 || nwithin == 2)
                        y0c = fabs(yint[0] - y0) < fabs(yint[1] - y0) ?  yint[0] : yint[1];

                      ycyc_int[ip - 2] = y0c; // points based on two different circle fit
                      yave += y0c;
                      nsec++;

                      // for debugging
                      ycyc[2*ip - 4] = yint[0];
                      ycyc[2*ip - 3] = yint[1];
                    }
                  }
                  else
                    ycyc_int[ip - 2] = y0;
                }

                if (kr > kr_max) {
                  // for the region with large gradient of curvature, choose the fit
                  // resulting in smaller curvature (or linear fitting)
                  yave = fabs(kp[0]) < fabs(kp[1]) ? nsec*ycyc_int[0]:  nsec*ycyc_int[1];
                }
                else if (nsec == 2 && (within(ycyc_int[0], 0., 1.) != within(ycyc_int[1], 0., 1.))) {
                  // when two circle fits result in different results, i.e., different cutting results
                  // within [0., 1.], we choose the result obatined with the circle with smaller curvature.
                  yave = fabs(kp[0]) < fabs(kp[1]) ? 2.*ycyc_int[0]:  2.*ycyc_int[1];
                }

                if (nsec > 0)
                  yave /= nsec;

                if (nsec == 0 && y0 >= 0. && y0 <= 1.) {
                  snew.x[] = y0;
                  with_marker.x[] += 1.; // should be removed
                }
                else if (nsec > 0 && yave >= 0. && yave <= 1.) {
                  with_circle = true;
                  snew.x[] = yave;
                  with_marker.x[] += 1.; // should be removed
                }
              }
              else {
                // case with two segment endpoints located on the same side of the cell face
                // intersection can only result from circle fitting.
                int nsec = 0;
                double yave = 0., ymin = min(pm[0].y, pm[1].y), ymax = max(pm[0].y, pm[1].y);
                double ycyc[4] = {-1., -1., -1, -1};
                double ycyc_int[2] = {HUGE, HUGE};
                for (int ip = 2; ip < np; ip++) {
                  double xrc, yrc, rc, yint[2];
                  get_circle (pm[0].x, pm[0].y, pm[1].x, pm[1].y, pm[ip].x, pm[ip].y, &xrc, &yrc, &rc);
                  if (rc > 0.) {
                    int nint = intersection_circle_two (xrc, yrc, rc, 0., yint);
                    if (nint > 0) {
                      double y0c;
                      int nwithin = 0;
                      for (int iint = 0; iint < nint; iint++) {
                        if (within(yint[iint], ymin, ymax) && within(yint[iint], 0., 1.)) {
                          coord dx1 = {pm[0].x - xrc, pm[0].y - yrc}, dx2 = {pm[1].x - xrc, pm[1].y - yrc},
                            dx3 = {-xrc, yint[iint] - yrc};
                          double within_arc = (dx1.x*dx3.y - dx1.y*dx3.x)*(dx2.x*dx3.y - dx2.y*dx3.x);
                          if (within_arc < 0.) {
                            nwithin++;
                            y0c = yint[iint];
                          }
                        }
                      }

                      ycyc[2*ip - 4] = yint[0];
                      ycyc[2*ip - 3] = yint[1];

                      if (nwithin > 0) {
                        ycyc_int[ip - 2] = y0c;
                        yave += y0c;
                        nsec++;
                      }
                    }
                  }
                }

                // when two circles result in inconsistent intersections,
                // we still choose the circle with smaller curvature
                if (nsec == 1) {
                  yave = fabs(kp[0]) < fabs(kp[1]) ? ycyc_int[0] : ycyc_int[1];
                  if (within(yave, 0., 1.)) {
                    with_circle = true;
                    snew.x[] = yave;
                    with_marker.x[] += 1.;
                  }
                }
                else if (nsec == 2) {
                  yave = ycyc_int[0] + ycyc_int[1];
                  with_circle = true;
                  snew.x[] = yave/nsec;
                  with_marker.x[] += 1.;
                }
              }

            }
          }
        }
    }
  }

  /** Update the corner color vertex.

    Change the color of vertex:
    (1) when marker moves across the grid line, for case with zero velocity component;
    (2) cell vertex falls inside the polygon formed by the two segments before
    and after advection. */
  foreach_vertex() {
    if (is_leaf(cell)) {
      // leaf cell only, vertex field on resolution boundary is special
      if (fabs(ss_tmp.y[]) < machine_zero && s_tmp.y[] < 0.)
        color_pha_new[] = color_pha[1];
      else if (fabs(ss_tmp.y[-1]) < machine_zero && s_tmp.y[-1] > 1.)
        color_pha_new[] = color_pha[-1];
      else if (fabs(ss_tmp.x[]) < machine_zero && s_tmp.x[] < 0.)
        color_pha_new[] = color_pha[0,1];
      else if (fabs(ss_tmp.x[0,-1]) < machine_zero && s_tmp.x[0,-1] > 1.)
        color_pha_new[] = color_pha[0,-1];
      else {
        for (int i = -1; i < 1; i++)
          for (int j = -1; j < 1; j++) {
            int conf = (int) config_dict[i,j];
            if (conf > 0) {
              Point pt = point;
              pt.i += i;
              pt.j += j;

              int nm[2], nm_o[2];
              coord xm[8], xm_o[8];
              int ns = get_segments (pt, config_dict, st = s_tmp, xm = xm, nm = nm, sn = ss_tmp);
              int ns_o = get_segments (pt, config_dict, st = s, xm = xm_o, nm = nm_o);

              for (int iseg = 0; iseg < ns; iseg++) {
                coord xv[8], xo = {0., 0.}, xs = {i, j};
                double kp[2] = {0., 0.};

                coord pm[4];
                int np = nm[iseg];

                for (int ip = 0; ip < np; ip++) {
                  foreach_dimension()
                    pm[ip].x = xm[4*iseg + ip].x + xs.x;

                  if (ip >= 2)
                    kp[ip - 2] = get_kappa_circle (pm[0], pm[1], pm[ip]);
                }
                double kr = max(fabs(kp[0]), fabs(kp[1]))/(min(fabs(kp[0]), fabs(kp[1])) + 1.e-32);

                foreach_dimension() {
                  xv[0].x = xm[4*iseg].x + xs.x;
                  xv[1].x = xm_o[4*iseg].x + xs.x;
                  xv[2].x = xm_o[4*iseg + 1].x + xs.x;
                  xv[3].x = xm[4*iseg + 1].x + xs.x;
                }

                coord xy_ave = {0., 0.}, n_ave = {0., 0.}, xy_rc = {0., 0.};
                coord xy_cyc[2] = {{HUGE, HUGE}, {HUGE, HUGE}};
                for (int ip = 2; ip < np; ip++) {
                  double rc, yint[2];
                  get_circle (pm[0].x, pm[0].y, pm[1].x, pm[1].y, pm[ip].x, pm[ip].y, &xy_rc.x, &xy_rc.y, &rc);

                  /** we ignore the case with 'pm[0].x * pm[1].x > 0', actually we don't need this to
                  construct a correct polygon. The most important thing is obtaining the instersection
                  with a consistent scheme. */
                  foreach_dimension() {
                    if (pm[0].x*pm[1].x <= 0) {
                      double ymin = min(pm[0].y, pm[1].y), ymax = max(pm[0].y, pm[1].y);
                      double y0 = my_intersect (pm[0].x, pm[0].y, pm[1].x, pm[1].y, 0., -HUGE, HUGE);
                      if (rc > 0.) {
                        int nint = intersection_circle_two (xy_rc.x, xy_rc.y, rc, 0., yint);
                        if (nint > 0) {
                          double y0c;
                          int nwithin = 0;
                          for (int iint = 0; iint < nint; iint++) {
                            if (yint[iint] >= ymin && yint[iint] <= ymax) {
                              nwithin++;
                              y0c = yint[iint];
                            }
                          }

                          if (nwithin == 0 || nwithin == 2)
                            y0c = fabs(yint[0] - y0) < fabs(yint[1] - y0) ?  yint[0] : yint[1];

                          xy_cyc[ip - 2].x = y0c; // intersection point with x = 0 or y = 0 line
                          xy_ave.x += y0c;
                          n_ave.x += 1;
                        }
                      }
                      else {
                        xy_cyc[ip - 2].x = y0;
                      }
                    }
                  }
                }

                // add the intersection points resulted from the circle fit
                int nv_poly = 4;
                foreach_dimension() {
                  int nsec = (int) n_ave.x;
                  if (kr > kr_max)
                    xy_ave.x = fabs(kp[0]) < fabs(kp[1]) ? nsec*xy_cyc[0].x : nsec*xy_cyc[1].x;
                  else if (nsec == 2 && (within(xy_cyc[0].x, 0., 1.) != within(xy_cyc[1].x, 0., 1.)))
                    xy_ave.x = fabs(kp[0]) < fabs(kp[1]) ? 2.*xy_cyc[0].x : 2.*xy_cyc[1].x;
                  if (nsec > 0) {
                    xv[nv_poly].x = 0.;
                    xv[nv_poly].y = xy_ave.x/n_ave.x;
                    nv_poly++;
                  }
                }

                if (nv_poly >= 6) {
                  coord dx1, dx2;
                  double is_reverse = 0.;
                  foreach_dimension() {
                    dx1.x = xv[0].x - xv[3].x;
                    dx2.x = xv[5].x - xv[4].x;
                    is_reverse += dx1.x*dx2.x;
                  }
                  if (is_reverse < 0.) {
                    foreach_dimension()
                      swap(double, xv[4].x, xv[5].x);
                  }
                }

                bool is_in = is_inside (xv, 4, xo);
                bool is_in_circle = is_inside (xv, nv_poly, xo);
                if (is_in_circle)
                  color_pha_new[] = 1. - color_pha[];
              }
            }
          }
      }
    }
  }
  boundary ((scalar *) {color_pha_new});

  // update the central color vertex
  update_color_cen_unsplit();
  foreach_vertex()
    color_pha[] = color_pha_new[];

  // determine the number of markers on edge based on color vertex
  foreach_face() {
    int with_face = fabs(color_pha[] - color_pha[0,1]) > machine_zero;
    if (with_face)
      with_marker.x[] = 1.;
    else { // without interface or with two interface
      snew.x[] = 0.;
      with_marker.x[] = 0.;
    }
  }

  foreach_face()
    s.x[] = snew.x[];

  boundary ((scalar *) {s}); // It's needed in MPI
}
