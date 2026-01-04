int main()
{
  origin (-0.5, -0.5);
  init_grid(8);
  refine (level == 3 && x*x + y*y < sq(0.25));
  output_cells (stdout);
  foreach_vertex()
    fprintf (stderr, "%g %g %d\n", x, y, level);
}
