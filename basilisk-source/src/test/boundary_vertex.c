/**
Checks boundary conditions for vertex fields in parallel. */

#define shape(x,y) (sq(x) + sq(y))

int main()
{
  size (16);
  init_grid (4);
  origin (-L0/2., -L0/2.);
  refine (sq(x) + sq(y) < sq(L0/4.) && level < 3);

  scalar s[];
  foreach_dimension() {
    s[left] = dirichlet(shape(x,y));
    s[right] = dirichlet(shape(x,y));
  }
  foreach()
    s[] = shape(x,y);
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (s[] + s[-1] + s[0,-1] + s[-1,-1])/4.;
  
  foreach() {
    fprintf (qerr, "%g %g %g\n", x - Delta/3., y - Delta/3., phi[]);
    fprintf (qerr, "%g %g %g\n", x - Delta/3., y + Delta/3., phi[0,1]);
    fprintf (qerr, "%g %g %g\n", x + Delta/3., y - Delta/3., phi[1,0]);
    fprintf (qerr, "%g %g %g\n", x + Delta/3., y + Delta/3., phi[1,1]);
  }

#if _MPI
  debug_mpi (NULL);
#else
  output_cells();
#endif
  
}
