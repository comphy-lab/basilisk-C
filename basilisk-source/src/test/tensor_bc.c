int main()
{
  init_grid(8);
  L0 = 2.*pi;
  vector u[];
  foreach() {
    u.x[] = sin(x)*sin(y);
    u.y[] = sin(x)*sin(y);
  }

  // symmetric: fixme: boundary conditions don't work for symmetric tensors!
    tensor D[];
  foreach() {
    foreach_dimension()
      D.x.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
    D.x.y[] = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/(2.*Delta);
  }
  output_cells (stdout);
  D.y.x.dirty = true; // fixme
  boundary ({D.y.x}); // fixme
  foreach()
    foreach_neighbor(1)
      fprintf (stderr, "%g %g %g %g %g %g %g\n",
	       x, y, D.x.x[], D.y.y[], D.x.y[],
	       cos(x)*sin(y), cos(y)*sin(x));
}
