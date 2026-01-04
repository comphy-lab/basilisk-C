/**
# Automatic stencils / boundary conditions */

int main()
{
  init_grid (8);

  scalar s[];

  /**
  This should warn (once) that we are assigning a scalar multiple
  times. */

  for (int i = 1; i < 2; i++)
    foreach_face()
      s[] = 1;

  /**
  This should set s as the y-component of a face vector. */

  foreach_face(y)
    s[] = 1;

  /**
  This should set v to be a face vector. */
  
  vector v[];
  foreach_face()
    v.x[] = 1;

  /**
  Same thing for the components of a tensor. */

  tensor Fq[];
  foreach_face()
    Fq.x.x[] = Fq.y.x[] = 0.;

  /**
  This should set Fv.x.y as a vertex scalar. */

  tensor Fv[];
  foreach_vertex()
    Fv.x.y[] = 1;
  
  /**
  Here flux boundary conditions should be imposed on s[]. */
  
  scalar a[];
  foreach()
    a[] = (s[0,1] - s[])/Delta;

  /**
  Here we use 'tangential' values of the face vector y-component 's',
  so central boundary conditions need to be imposed. */
  
  foreach()
    a[] = (s[1] - s[])/Delta;

  /**
  If normal boundary conditions are set, central boundary conditions
  should be imposed (not just fluxes). */

  foreach_face(y)
    s[] = 1;

  s[top] = 0.;
  foreach()
    a[] = (s[0,1] - s[])/Delta;  

  /**
  The same but for face vectors. */

  face vector f[];
  foreach_face()
    f.x[] = 1;

  /**
  Flux BCs here. */
  
  foreach()
    a[] = (f.y[0,1] - f.y[])/Delta;

  /**
  Central BCs here. */
  
  f.n[top] = 0.;
  foreach()
    a[] = (f.y[0,1] - f.y[])/Delta;  
  
  /**
  This should warn that we are assigning a face vector using a
  foreach() loop. */

  foreach()
    f.x[] = 2;

  /**
  This should warn that we are assigning a vertex scalar using a
  foreach() loop. */

  vertex scalar b[];
  foreach()
    b[] = 2;

  /**
  This should warn that we are using a global variable without
  reduction. */

  double t = 0;
  foreach()
    t++;
}
