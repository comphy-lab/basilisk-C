/**
Checks that dimensions of conditional expressions are properly set. */

int main()
{
  init_grid (1);
  scalar f[];
  foreach() {
    double a;
    if (a)
      f[] = 1. [1];
    else
      f[] = 2.;
  }
  foreach()
    f[] = f[] < 3. ? 4. : f[];

  /**
  Multiplicative conditional expressions with a constant must be
  dimensionless. So [d] should be [c] since [(b > 0 ? 1 : -2)] = [0]
  since it is a multiplicative constant. */

  double b, c = 7, d;
  d = c*(b > 0 ? 1 : -2);
  display_value (d);

  /**
  Same here, but [e] must also be zero and so [4] is zero. */
  
  double e = 4;
  d = c*(b > 0 ? 3 : e);
  display_value (d);

  /**
  Values of undefined conditional expressions must have the same
  dimensions (equivalent to [test1.c]() for undefined branching). */

  double d = b ? 1 [2] : 2;
  3. + d; // [3] should be [2]
}
