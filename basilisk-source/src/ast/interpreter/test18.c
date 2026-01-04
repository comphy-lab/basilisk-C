/**
Checks some math functions. **/

int main()
{
  double c;
  {
    interpreter_verbosity (4);
    double a = 1, b = 2;

    double c = a*b;
    double d = sqrt (c);
    double e = pow (d, 2);
    double f = fabs (e);
    double g = sq(f);
    double a = 1.;
    a *= g;
    double b = sin (a/g);
  }
  
  init_grid (1);

  scalar s[], f[];
  foreach() {
    interpreter_verbosity (4);
    double a = L0*sin (x/L0);
    s[] = a*a;
    f[] = 2;

  }

  foreach() {
    interpreter_verbosity (4);
    f[] + L0*L0;
    s[] + L0*L0;
  }
}
