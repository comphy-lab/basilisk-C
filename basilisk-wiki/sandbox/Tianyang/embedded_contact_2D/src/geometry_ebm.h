#define IATOL 1e-10

int facets_ebm (coord n, double alpha, coord p[2])
{
  int nx = 0;
  int ny = 0;
  coord px[2];
  coord py[2];

  for (double s = -0.5; s <= 0.5; s += 1.) {
    if (fabs (n.y) > 0.) {
      double a = (alpha - s*n.x)/n.y;
      if (a >= -0.5 && a <= 0.5) {
        px[nx].x = s;
	      px[nx++].y = a;
      }
    }
    if (fabs (n.x) > 0.) {
      double a = (alpha - s*n.y)/n.x;
      if (a >= -0.5 && a <= 0.5) {
        py[ny].y = s;
	      py[ny++].x = a;
      }
    }
  }

  if (nx == 2) {
    foreach_dimension() {
      p[0].x = px[0].x;
      p[1].x = px[1].x;
    }
    return 2;
  }
  else if (ny == 2) {
    foreach_dimension() {
      p[0].x = py[0].x;
      p[1].x = py[1].x;
    }
    return 2;
  }
  else {
    int i = 0;
    if (nx > 0) {
      assert (nx == 1);
      p[i].x = px[0].x;
      p[i].y = px[0].y;
      i++;
    }
    if (ny > 0) {
      assert (ny == 1);
      p[i].x = py[0].x;
      p[i].y = py[0].y;
      i++;
    }
    assert (i <= 2);
    return i;
  }
}

#if dimension == 2

double area_line_triangle (coord p[3], coord n, double alpha)
{
  assert (n.x != 0. || n.y != 0.);
  assert (p[1].x == 0.5 && p[1].y == 0.5); 
  assert (p[0].y == 0.5);
  assert (p[2].x == 0.5);

  double area = 0.5*(0.5 - p[0].x)*(0.5 - p[2].y);
  double d0 = n.x*p[0].x + n.y*p[0].y - alpha;
  double d1 = n.x*p[1].x + n.y*p[1].y - alpha;
  double d2 = n.x*p[2].x + n.y*p[2].y - alpha;

  if (d0 >= 0. && d1 >= 0. && d2 >= 0.)
    return 0.;
  else if (d0 <= 0. && d1 <= 0. && d2 <= 0.)
    return area;
  else if (d0 == 0. || d1 == 0. || d2 == 0.) {
    if (d0 == 0.) {
      assert (d1*d2 < 0.);
      assert (fabs(n.y) > 0.);
      double y0 = (alpha - 0.5*n.x)/n.y;
      assert (y0 > p[2].y && y0 < 0.5);

      return (d2 > 0. ? 0.5 - y0 : y0 - p[2].y)*(0.5 - p[0].x)/2.;
    }
    else if (d2 == 0.) {
      assert (d0*d1 < 0.);
      assert (fabs(n.x) > 0.);
      double x0 = (alpha - 0.5*n.y)/n.x;
      assert (x0 > p[0].x && x0 < 0.5);

      return (d0 > 0. ? 0.5 - x0 : x0 - p[0].x)*(0.5 - p[2].y)/2.;
    }
    else { //d1 == 0
      assert (d0*d2 < 0.);
      assert (fabs((p[2].x - p[0].x)*n.x + (p[2].y - p[0].y)*n.y) > 0.);
      double x0 = ((p[2].x - p[0].x)*alpha - n.y*(p[2].x*p[0].y - p[0].x*p[2].y))/((p[2].x - p[0].x)*n.x + (p[2].y - p[0].y)*n.y);
      assert (x0 > p[0].x && x0 < 0.5);

      return (d0 > 0. ? 0.5 - x0 : x0 - p[0].x)*(0.5 - p[2].y)/2.;
    }  
  }
  else {
    assert ((d0 != 0. && d1 != 0. && d2 != 0.));
    if (d0*d1 < 0. && d1*d2 < 0.) {
      assert (fabs(n.x) > 0 && fabs(n.y) > 0.);
      double x0 = (alpha - 0.5*n.y)/n.x;
      double y0 = (alpha - 0.5*n.x)/n.y;

      return (d1 < 0. ? (0.5 - x0)*(0.5 - y0)/2. : area - (0.5 - x0)*(0.5 - y0)/2.);
    }
    else if (d1*d2 < 0. && d2*d0 < 0.) {
      assert (fabs(n.y) > 0.);
      double y0 = (alpha - 0.5*n.x)/n.y;
      assert (fabs((p[2].x - p[0].x)*n.x + (p[2].y - p[0].y)*n.y) > 0.);
      double x0 = ((p[2].x - p[0].x)*alpha - n.y*(p[2].x*p[0].y - p[0].x*p[2].y))/((p[2].x - p[0].x)*n.x + (p[2].y - p[0].y)*n.y);

      return (d2 < 0. ? (y0 - p[2].y)*(0.5 - x0)/2. : area - (y0 - p[2].y)*(0.5 - x0)/2.);
    }
    else {
      assert (d2*d0 < 0. && d0*d1 < 0.);
      assert (fabs(n.x) > 0.);
      double x0 = (alpha - 0.5*n.y)/n.x;
      assert (fabs((p[2].x - p[0].x)*n.x + (p[2].y - p[0].y)*n.y) > 0.);
      double y0 = ((p[2].y - p[0].y)*alpha + n.x*(p[2].x*p[0].y - p[0].x*p[2].y))/((p[2].x - p[0].x)*n.x + (p[2].y - p[0].y)*n.y);

      return (d0 < 0. ? (x0 - p[0].x)*(0.5 - y0)/2. : area - (x0 - p[0].x)*(0.5 - y0)/2.);
    }
  }
}

double area_line_rectangle (coord p[4], coord n, double alpha)
{
  assert (p[2].x == 0.5 && p[2].y == 0.5);

  double area;
  if (p[1].y == 0.5) {
    assert (p[0].y >= p[3].y);
    assert (p[1].x == -0.5 && p[1].y == 0.5);
    assert (p[0].x == -0.5);
    assert (p[3].x ==  0.5);

    area = rectangle_fraction (n, alpha, p[0], p[2])*(p[1].y - p[0].y);

    if (p[0].y > p[3].y) {
      coord ptri[3];
      ptri[0] = (coord){-0.5, 0.5};
      ptri[1] = (coord){ 0.5, 0.5};
      ptri[2] = (coord){ 0.5, p[3].y - p[0].y + 0.5};
      alpha += (0.5 - p[0].y)*n.y;

      area += area_line_triangle (ptri, n, alpha);
    }
  }
  else {
    assert (p[0].x >= p[3].x);
    assert (p[1].x == 0.5 && p[1].y == -0.5);
    assert (p[3].y ==  0.5);
    assert (p[0].y == -0.5);

    area = rectangle_fraction (n, alpha, p[0], p[2])*(p[1].x - p[0].x);

    if (p[0].x > p[3].x) {
      coord ptri[3];
      ptri[0] = (coord){p[3].x - p[0].x + 0.5, 0.5};
      ptri[1] = (coord){0.5,  0.5};
      ptri[2] = (coord){0.5, -0.5};
      alpha += (0.5 - p[0].x)*n.x;

      area += area_line_triangle (ptri, n, alpha);
    }
  }

  return area;
}

double area_line_pentagon (coord p[5], coord n, double alpha)
{
  assert (p[1].x == -0.5 && p[1].y ==  0.5);
  assert (p[2].x ==  0.5 && p[2].y ==  0.5);
  assert (p[3].x ==  0.5 && p[3].y == -0.5);
  assert (p[0].x == -0.5 && p[4].y == -0.5);

  double area = plane_volume (n, alpha);

  foreach_dimension()
    n.x = - n.x;

  coord ptri[3];
  foreach_dimension() {
    ptri[0].x = - p[4].x;
    ptri[2].x = - p[0].x;
  }
  ptri[1] = (coord){0.5, 0.5};
  
  area -= area_line_triangle (ptri, n, alpha);

  return area;
}

double line_alpha_ebm (double c, double cs, coord ncs, double alphacs, coord n)
{
  assert (c > 0 && c < cs && cs < 1.);

  foreach_dimension()
    if (ncs.x > 0.) {
      ncs.x = - ncs.x;
      n.x = - n.x;
    }

  coord pcs[2];
  assert (facets_ebm (ncs, alphacs, pcs) == 2);
  if (pcs[0].x > pcs[1].x)
    swap (coord, pcs[0], pcs[1]);
  assert (pcs[0].x <= pcs[1].x);

  int vnum = 0;
  for (double sx = -0.5; sx <= 0.5; sx += 1.)
    for (double sy = -0.5; sy <= 0.5; sy += 1.)
      if (ncs.x*sx + ncs.y*sy - alphacs < 0.)
        vnum++;
  assert (vnum > 0 && vnum < 4);
  assert (ncs.x*0.5 + ncs.y*0.5 - alphacs < 0.);

  coord polygon[5];
  double (* area_line_polygon) (coord *, coord, double);
  if (vnum == 1) {
    polygon[0] = pcs[0];
    polygon[1] = (coord){0.5, 0.5};
    polygon[2] = pcs[1];

    area_line_polygon = area_line_triangle;
  }
  else if (vnum == 2) {
    if (ncs.x > ncs.y) {
      assert (ncs.x*(-0.5) + ncs.y*0.5 - alphacs < 0.);
      polygon[0] = pcs[0];
      polygon[1] = (coord){-0.5, 0.5};
      polygon[2] = (coord){ 0.5, 0.5};
      polygon[3] = pcs[1];
    }
    else {
      assert (ncs.x < ncs.y);
      assert (ncs.x*0.5 + ncs.y*(-0.5) - alphacs < 0.);
      polygon[0] = pcs[1];
      polygon[1] = (coord){0.5, -0.5};
      polygon[2] = (coord){0.5,  0.5};
      polygon[3] = pcs[0];
      if (polygon[0].y > polygon[3].y)
        swap (coord, polygon[0], polygon[3]);
    }

    area_line_polygon = area_line_rectangle;
  }
  else { //vnum ==3
    assert (ncs.x*(-0.5) + ncs.y*0.5 - alphacs < 0.);
    assert (ncs.x*0.5 + ncs.y*(-0.5) - alphacs < 0.);

    polygon[0] = pcs[0];
    polygon[1] = (coord){-0.5,  0.5};
    polygon[2] = (coord){ 0.5,  0.5};
    polygon[3] = (coord){ 0.5, -0.5};
    polygon[4] = pcs[1];
    
    area_line_polygon = area_line_pentagon;
  }

  double alpha_polygon[5];
  for (int k = 0; k < vnum+2; k++)
    alpha_polygon[k] = n.x*polygon[k].x + n.y*polygon[k].y;
  
  double alpha_min = nodata;
  for (int k = 0; k < vnum+2; k++)
    if (area_line_polygon (polygon, n, alpha_polygon[k]) < c) {
      alpha_min = alpha_polygon[k];
      break;
    }

  double alpha_max = nodata;
  for (int k = 0; k < vnum+2; k++)
    if (area_line_polygon (polygon, n, alpha_polygon[k]) > c) {
      alpha_max = alpha_polygon[k];
      break;
    }

  assert (alpha_min != nodata && alpha_max != nodata);

  int iternum = 0;
  double alpha_iter = 0.5*(alpha_min + alpha_max);
  double area = area_line_polygon (polygon, n, alpha_iter);
  while (fabs (area - c) > IATOL && iternum <= 100) {
    if (area > c)
      alpha_max = alpha_iter;
    else
      alpha_min = alpha_iter;

    alpha_iter = 0.5*(alpha_min + alpha_max);
    area = area_line_polygon (polygon, n, alpha_iter);
    iternum++;
  }
  
  if (iternum > 100)
    fprintf (stdout, "warning: iternum > 100\n");

  if (fabs (area - c) <= IATOL)
    return alpha_iter;
  else
    return nodata;
}

double line_area_ebm (coord pcs[2], double cs, coord ncs, double alphacs, coord n, double alpha)
{
  assert (cs > 0. && cs < 1.);

  foreach_dimension()
    if (ncs.x > 0.) {
      ncs.x = - ncs.x;
      n.x = - n.x;

      pcs[0].x = - pcs[0].x;
      pcs[1].x = - pcs[1].x;
    }

  if (pcs[0].x > pcs[1].x)
    swap (coord, pcs[0], pcs[1]);
  assert (pcs[0].x <= pcs[1].x);

  int vnum = 0;
  for (double sx = -0.5; sx <= 0.5; sx += 1.)
    for (double sy = -0.5; sy <= 0.5; sy += 1.)
      if (ncs.x*sx + ncs.y*sy - alphacs < 0.)
        vnum++;
  assert (vnum > 0 && vnum < 4);
  assert (ncs.x*0.5 + ncs.y*0.5 - alphacs < 0.);

  coord polygon[5];
  double (* area_line_polygon) (coord *, coord, double);
  if (vnum == 1) {
    polygon[0] = pcs[0];
    polygon[1] = (coord){0.5, 0.5};
    polygon[2] = pcs[1];

    area_line_polygon = area_line_triangle;
  }
  else if (vnum == 2) {
    if (ncs.x > ncs.y) {
      assert (ncs.x*(-0.5) + ncs.y*0.5 - alphacs < 0.);
      polygon[0] = pcs[0];
      polygon[1] = (coord){-0.5, 0.5};
      polygon[2] = (coord){ 0.5, 0.5};
      polygon[3] = pcs[1];
    }
    else {
      assert (ncs.x < ncs.y);
      assert (ncs.x*0.5 + ncs.y*(-0.5) - alphacs < 0.);
      polygon[0] = pcs[1];
      polygon[1] = (coord){0.5, -0.5};
      polygon[2] = (coord){0.5,  0.5};
      polygon[3] = pcs[0];
      if (polygon[0].y > polygon[3].y)
        swap (coord, polygon[0], polygon[3]);
    }

    area_line_polygon = area_line_rectangle;
  }
  else { //vnum ==3
    assert (ncs.x*(-0.5) + ncs.y*0.5 - alphacs < 0.);
    assert (ncs.x*0.5 + ncs.y*(-0.5) - alphacs < 0.);

    polygon[0] = pcs[0];
    polygon[1] = (coord){-0.5,  0.5};
    polygon[2] = (coord){ 0.5,  0.5};
    polygon[3] = (coord){ 0.5, -0.5};
    polygon[4] = pcs[1];
    
    area_line_polygon = area_line_pentagon;
  }

  return area_line_polygon (polygon, n, alpha);
}

double rectangle_fraction_cs (double cs, double fs, coord ncs, double alphacs, coord n, double alpha, coord a, coord b)
{
  coord n1, ncs1;
  foreach_dimension() {
    alpha -= n.x*(b.x + a.x)/2.;
    n1.x = n.x*(b.x - a.x);

    alphacs -= ncs.x*(b.x + a.x)/2.;
    ncs1.x = ncs.x*(b.x - a.x);
  }

  coord pcs[2];
  int num = facets_ebm (ncs1, alphacs, pcs);

  if (num < 2)
    return plane_volume (n1, alpha);

  return line_area_ebm (pcs, cs, ncs1, alphacs, n1, alpha);
}

#endif

