/**
# Boundary conditions for layers 

Note that the boundary conditions for f.y and v.x in the bottom-left
and bottom-right corners do not seem to be consistent, for quadtrees
only. This could be due to the complicated and obsolete treatment of
masked cells and should be checked when this is
removed. [boundaries2.c]() should be checked similarly. */

#define BGHOSTS 2
#define LAYERS 1

int main()
{
  nl = 2;
  
  init_grid (4);
  output_cells (stdout);

  vector v = new vector[nl];
  scalar s = new scalar[nl];
  face vector f = new face vector[nl];

  s[left] = dirichlet(x + y);
  s[top] = dirichlet(x + y);
  s[right] = dirichlet(x + y);
  s[bottom] = dirichlet(x + y);

  v.t[bottom] = dirichlet(0);
  f.n[bottom] = 0.;

  foreach()
    foreach_layer() {
      s[] = x + y;
      v.x[] = v.y[] = 1.;
    }
  foreach_face()
    foreach_layer()
      f.x[] = 1.;

  foreach()
    foreach_layer()
      foreach_neighbor()
        fprintf (stderr, "s %d %g %g %g %g %g\n",
		 point.l, x, y, s[], v.x[], v.y[]);  
  foreach_face(x)
    foreach_layer()
      for (int i = -2; i <= 2; i++)
        fprintf (stderr, "f.x %d %g %g %g\n", point.l, x, y + i*Delta, f.x[0,i]);
  foreach_face(y)
    foreach_layer()
      for (int i = -2; i <= 2; i++)
	fprintf (stderr, "f.y %d %g %g %g\n", point.l, x + i*Delta, y, f.y[i]);
}
