/**
# Vertex boundary conditions on multigrids

it [works](vertex_boundaries/out)
*/
vertex scalar v[]; 

int main() {
  init_grid (4);
  face vector vf;
  foreach_dimension()
    vf.x = v;
  v.v = vf;
  v.face = true;
  v[left] = 1;
  v[right] = 1;
  // boundary ({v}); // automatic!
  foreach_vertex()
    printf ("%g\t%g\n", x, v[]);
}

