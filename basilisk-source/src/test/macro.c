/**
# Macros */

macro myothermacro()
{
  for (int j = 1; j < 3; j++) {
    fprintf (stderr, "other %d ", j);
    {...}
  }
}

macro mymacro (int len = 0, int end = 0)
{
  int i;
  for (i = len; i < end; i++)
    myothermacro()
      {...}
  {...}
}

macro mymacro2 (coord * view, float color[3])
{
  view->x = 1;
  view->y = 2;
  color[0] = 2;
  {...}
}

macro mymacro3 (double field, scalar s1, double expr)
{
  foreach() {
    field = expr;
    s1[] = expr;
    fprintf (stderr, "%g %g %g\n", x, y, expr);
  }
}

macro mymacro4() {}

macro mymacro4 (break = (i = j = 2)) // must overload the previous definition
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      {...}
}

auto macro mymacro4() {} // must not overload the previous definition

macro int myreturn (double x)
{
  return x*x;
}

macro
double myreturn2 (Point point, double expr)
{
#if 0 // simple macro
  return expr;
#else // complex macro
  double b = expr;
  return b;
#endif
}

macro
void withreturn (int i)
{
  if (i) {
    fprintf (stderr, "i != 0\n");
    return;
  }
  fprintf (stderr, "i == 0\n");  
}

macro m1() {
  fprintf (stderr, "m1\n");
}

macro m2() {
  m1();
}

void print_m2()
{
  m2(); // should print "m1"
}

macro m1() {
  fprintf (stderr, "m1bis\n");
}

void print_m2bis()
{
  m2(); // should print "m1bis"
}

/**
## Macro "inheritance" */

macro rmacro()
{
  fprintf (stderr, "rmacro\n");
  {...}
}

macro rmacro()
{
  rmacro() {{...}}
  fprintf (stderr, "rmacro1\n");
  {...}
}

int main()
{
  mymacro (end = 2)
    fprintf (stderr, "%d\n", i);
  coord p;
  float machin[3];
  mymacro2 (&p, machin)
    machin[1] = 3;
  fprintf (stderr, "%g %g %g %g\n", p.x, p.y, machin[0], machin[1]);

  init_grid (2);
  scalar s[], s1[];
  mymacro3 (s[], s1, x*x + y*y);
  foreach()
    fprintf (stderr, "%g %g %g %g\n", x, y, s[], s1[]);
  
  mymacro4() {
    for (int k = 0; k < 2; k++) {
      fprintf (stderr, "%d %d %d\n", i, j, k);
      break;
    }
    if (j == 1)
      break;
  }

  @define DEF M_PI
  mymacro3 (s[], s1, DEF);

  fprintf (stderr, "%d\n", myreturn(12));

  foreach_face()
    fprintf (stderr, "%g %g %g %g\n", x, y,
	     myreturn2 (point, x + y),
	     myreturn2 (point, 2.*(x + y)));

  withreturn (0);
  withreturn (1);

  print_m2();
  print_m2bis();

  rmacro()
    fprintf (stderr, "rmacro test\n");
}
