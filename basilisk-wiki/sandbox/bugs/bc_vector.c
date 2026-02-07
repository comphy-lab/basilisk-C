/**
# Boundary conditions on a single component of a vector field are
  incorrect on trees */

int main(){
  init_grid(1);

  vector u[];
  u.n[bottom] = dirichlet(0);
  
  foreach()
    u.y[] = 1;

  /**
  Note that this is not really due to automatic BCs, which are
  correctly applied (on a single component), but to the implementation
  of boundary conditions on trees (this works on multigrid). This is
  most probably due to the deprecated `mask()` feature. */
  
  // boundary ({u.y}); // this still does not work
  // boundary ((scalar *){u}); // this works
  
  foreach() {
    fprintf (stderr, "%g\n", u.y[0,-1]);
    assert ((u.y[0,-1] == -1.));
  }
}
