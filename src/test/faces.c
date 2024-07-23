int main()
{
  origin (-0.5, -0.5);
  init_grid(4);
  refine (level == 2 && x*x + y*y < sq(0.25));
  foreach_face(x)
    fprintf (stderr, "%g %g\n", x, y);
  foreach_face(y)
    fprintf (stderr, "%g %g\n", x, y);
  foreach_face()
    fprintf (stderr, "%g %g\n", x, y);
  output_cells (stdout);
}
