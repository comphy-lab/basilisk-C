/**
# Checks that a field used in a boundary condition is properly indexed */

attribute { int dummy; }

scalar a[], b[];
a[left] = b[]*b.dummy;

int main()
{
  init_grid (1);

  foreach()
    b[] = 2;
  b.dummy = -1;

  foreach()
    a[] = 3;
  foreach (serial)
    fprintf (stderr, "%g %g\n", a[], a[-1]);
}
