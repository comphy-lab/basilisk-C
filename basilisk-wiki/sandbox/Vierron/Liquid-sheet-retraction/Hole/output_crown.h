/*
This code is the fonction Output_facets without 2D segments
**/
struct OutputCrown {
  scalar c;
  FILE * fp;     // optional: default is stdout
  face vector s; // optional: default is none
};

trace
void output_crown (struct OutputCrown p)
{
  scalar c = p.c;
  face vector s = p.s;
  if (!p.fp) p.fp = stdout;
  if (!s.x.i) s.x.i = -1;

  foreach()
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha = plane_alpha (c[], n);
      coord v[12];
      int m = facets (n, alpha, v, 1.);
      for (int i = 0; i < m; i++){
      	  if(z + v[i].z*Delta==0) // Find (x,y) interface coordinate at z=0
	  	fprintf (p.fp, "%g %g %g\n",
			 x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
	}
     if (m > 0)
	fputc ('\n', p.fp);
    }

  fflush (p.fp);
}