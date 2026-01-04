/**
# "Undefined branched expressions" have the same dimensions */

int main()
{
  init_grid (1);
  {
    double a, b;

    if (a)
      b = 1. [1];
    else
      b = 2.;
    display_value (b);

    scalar s[];
    foreach() {
      if (a)
	s[] = 1. [1];
      else
      	s[] = 0.;
      display_value (s[]);
    }

    if (a)
      b = 1. [1];
    else
      b = 1. [2]; // should give an error since [2] != [1]

  }
}
