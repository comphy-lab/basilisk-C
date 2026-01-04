/**
Tests undefined values of structure members. */

int main()
{{ interpreter_verbosity (4);
  coord p;
  p.x = 2.;
  double a, b;
  a = p.x;
  b = p.y;
}}
