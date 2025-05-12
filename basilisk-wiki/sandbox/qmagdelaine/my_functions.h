/**
Few functions that I often use.

A function to rescale normals so that they are unit vectors w.r.t. the
2-norm (by default, the 1-norm is adopted for efficiency purposes). */

coord normal (Point point, scalar c) {
  coord n = mycs (point, c);
  double nn = 0.;
  foreach_dimension()
    nn += sq(n.x);
  nn = sqrt(nn);
  foreach_dimension()
    n.x /= nn;
  return n;
}

/**
A function to compute 2-norm normal in every cell. */

void compute_normal (scalar f, vector normal_vector) {
  foreach() {
    coord n = normal (point, f);
    foreach_dimension() 
      normal_vector.x[] = n.x;
  }
  boundary((scalar*){normal_vector});
}

/**
A function to suppress glitches after an advection. */

void magnet (scalar f, double error) {
  foreach() {
    f[] = clamp(f[], 0., 1.);
    f[] = (f[] < error ? 0. : (f[] > 1. - error ? 1. : f[]));
  }
  boundary ({f});
}

/**
A function to compute in each point the divergence of a gradient based flux. */

void my_laplacian (scalar f, scalar l, face vector D) {
  boundary({f, D});
  foreach() {
    l[] = 0.;
    foreach_dimension()
      l[] += (f[1] - f[0])*D.x[1] - (f[] - f[-1])*D.x[];
    l[] /= sq(Delta);
  }
  boundary({l});
}

/**
A function to mesure the length of the interface in the cell. Warning: the
length is normalised by the size of the cell. To get the real interface length
you have to multiplie it by the cell size $\Delta$. */

double interface_length (Point point, scalar c)
{
  coord n = mycs (point, c);
  double alpha = line_alpha (c[], n);
  coord coord_centroid = {0, 0};
  return line_length_center(n, alpha, &coord_centroid);
}
