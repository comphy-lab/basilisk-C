#define BGHOSTS 2

int main() {
  init_grid (4);

  foreach_dimension()
    periodic (right);
  
  int l = 0;
  
  scalar s[];

  foreach_level(l)
    s[] = 0;
  boundary_level ({s}, l);

  foreach_boundary_level (l)
    printf ("%g %g %g\n", x, y, s[]);
  
  foreach_level (l + 1) {
    for (int i = -2; i <= 2; i++)
      for (int j = -2; j <= 2; j++)
	fprintf (stderr, "%g ", coarse(s,i,j));
    fprintf (stderr, "\n");
    exit (0);
  }
}
