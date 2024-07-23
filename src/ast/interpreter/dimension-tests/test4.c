/**
Constants in functions can change dimensions between calls. */

double sum (double x)
{
  return x + 1.; // 1. has the dimension of x, so will change between
		 // the calls below
}

int main()
{
  double a = 2;
  a += sum (1. [1]);
  sum (2. [2]);
#if 0
  double b = 1.;
  sum (sq(b));
#endif
}
