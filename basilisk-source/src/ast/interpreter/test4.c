/**
Tests array initialization. */

int main()
{{ interpreter_verbosity (4);
  double b[2] = {1, 2};
  double a = b[0]*3;
  a;
  coord p[1];
  (*p).x = b[1];
  a = p->x*2;
}}
