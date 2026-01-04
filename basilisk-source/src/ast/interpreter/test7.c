/**
A complex function with undefined conditions and loops. */

#define radius 1./12.
#define length 0.025
int maxlevel = 10;

int main()
{
  init_grid (64);
  origin (0, -1.5, -1.5);
  size (3.);
#if 1
  {
    interpreter_verbosity (2);
    int refined;
    do {
      boundary (all);
      refined = 0;
      tree->refined.n = 0;
      foreach_leaf()
	if (x < 1.2*length && sq(y) + sq(z) < 2.*sq(radius) && level < maxlevel) {
	  refine_cell (point, all, 0, &tree->refined);
	  refined++;
	  continue;
	}
    } while (refined);
  }
#else
  refine (x < 1.2*length && sq(y) + sq(z) < 2.*sq(radius) && level < maxlevel);
#endif
}
