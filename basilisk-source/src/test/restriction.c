// checks that restriction works for face vectors on a non-trivial tree grid

int main(int argc, char * argv[])
{
  X0 = Y0 = -0.7;
  init_grid (1);
  int depth = argc > 1 ? atoi(argv[1]) : 4;
  refine (level < depth - 2 ||
	  level <= depth*(1. - sqrt(sq(x - 0.1) + sq(y - 0.1))));

  face vector u[];
  foreach_face()
    u.x[] = 1.;
  boundary ((scalar *){u});

  output_cells (stdout);

  for (int l = 0; l <= depth(); l++)
    foreach_level(l) {
      fprintf (stderr, "%g %g %g %d\n", x - Delta/2., y, u.x[], l);
      assert ((u.x[] == 1. && u.y[] == 1. && u.x[1,0] == 1. && u.y[0,1] == 1.));
    }
}
