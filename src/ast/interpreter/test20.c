/**
Checks that freed values modified by undefined conditionals are
properly handled. This should be run with valgrind. */

int main()
{
  double a;
  if (a) {
    double * b = malloc (2*sizeof (double));
    b[0] = 1.;
    b[1] = 2.;
    free (b);
  }
}
