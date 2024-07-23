/**
Tests undefined conditions. */

int main()
{{ interpreter_verbosity (4);
  double a, b = 2, c = 1;
  if (a)
    c++;
  else
    c--;
  c;
  if (c)
    b;
  else
    b*b;
}}
