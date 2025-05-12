/**
# Redistancing a levelset function

This function reconstructs the zero contour of the distance function
*d* and recomputes *d* as the signed distance to this contour, using
the *distance()* function in [/src/distance.h](). */

#include "distance.h"

trace
void redistance (scalar d)
{
  assert (dimension == 2);

  Array * a = array_new();
  foreach_vertex() {
    coord p[4];
    int n = 0;
    if (d[-1,-1]*d[0,-1] < 0.)
      p[n].x = x - Delta/2. + Delta*(- d[-1,-1]/(d[0,-1] - d[-1,-1])),
	p[n].y = y - Delta/2., p[n++].z = 0.;
    if (d[0,-1]*d[0,0] < 0.)
      p[n].y = y - Delta/2. + Delta*(- d[0,-1]/(d[0,0] - d[0,-1])),
	p[n].x = x + Delta/2., p[n++].z = 0.;
    if (d[-1,0]*d[0,0] < 0.)
      p[n].x = x - Delta/2. + Delta*(- d[-1,0]/(d[0,0] - d[-1,0])),
	p[n].y = y + Delta/2., p[n++].z = 0.;
    if (d[-1,-1]*d[-1,0] < 0.)
      p[n].y = y - Delta/2. + Delta*(- d[-1,-1]/(d[-1,0] - d[-1,-1])),
	p[n].x = x - Delta/2., p[n++].z = 0.;
    if (n == 2 || n == 4)
      for (int i = 0; i < n/2; i++) {
	if (d[-1,-1] < 0.) {
	  array_append (a, &p[2*i], sizeof (coord));
	  array_append (a, &p[2*i+1], sizeof (coord));
	}
	else {
	  array_append (a, &p[2*i+1], sizeof (coord));
	  array_append (a, &p[2*i], sizeof (coord));
	}
      }
  }
  coord last = {nodata};
  array_append (a, &last, sizeof(coord));
  scalar dn[];

#if 0
  coord * p = a->p;
  while (p->x != nodata)
    fprintf (stderr, "%g %g %g %g\n", p->x, p->y, (p + 1)->x, (p + 1)->y),
      p += 2;
  exit (0);
#endif
  
  distance (dn, a->p);
  free (a);

  foreach()
    //    if (fabs(dn[]) > Delta)
      d[] = dn[];
  boundary ({d});
}
