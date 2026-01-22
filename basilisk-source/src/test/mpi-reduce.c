#include "utils.h"

int main ()
{
  init_grid (64);

  scalar s[];
  foreach()
    s[] = x + y;

  // statsf() uses reduction operations
  stats stat = statsf (s);
  fprintf (qerr, "%g %g %g\n", stat.min, stat.sum, stat.max);

  // Array reduction
  #define arr_size 10
  int cells[arr_size] = {0};
  foreach (reduction(+:cells[:arr_size])) 
    cells[(int)(10*fabs(x))]++;

  for (int i = 0; i < arr_size; i++) 
    fprintf (qerr, "%d ", cells[i]);
  fputc ('\n', qerr);

  // Coord reduction
  coord p = {0};
  mat3 p2 = {0};
  foreach (reduction(+:p) reduction(+:p2))
    foreach_dimension() {
      p.x++;
      p2.x.x++;
    }
  fprintf (qerr, "%g %g %g %g\n", p.x, p.y, p2.x.x, p2.y.y);

  // test array of coord and mat3 reduction 
  #define arr_size 10
  coord P[arr_size] = {{0}};
  mat3 P2[arr_size] = {{{0}}};
  foreach (reduction(+:P[:arr_size]) reduction(+:P2[:arr_size]))
    foreach_dimension() {
      P[(int)(10*fabs(x))].x++;
      P2[(int)(10*fabs(x))].x.x++;
    }
  
  char * _x="x",* _y="y";
  foreach_dimension(){
    fprintf (qerr,"P.%s : ", _x);
    for (int i = 0; i < arr_size; i++) 
        fprintf (qerr, "%g ", P[i].x);
    fputc ('\n', qerr);
    fprintf (qerr,"P.%s.%s : ", _x,_x);
    for (int i = 0; i < arr_size; i++) 
        fprintf (qerr, "%g ", P2[i].x.x);
    fputc ('\n', qerr);
  }

  // test that MPI sum reductions work also when the initial value is
  // not zero
  double sum = 1.;
  foreach(reduction(+:sum));
  assert (sum == 1.);

  {
    coord p = {1,2,3};
    foreach(reduction(+:p));
    assert (p.x == 1 && p.y == 2 && p.z == 3);
  }
}
