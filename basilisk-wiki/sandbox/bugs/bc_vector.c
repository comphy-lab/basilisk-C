vector u[];

u.n[bottom] = dirichlet(0);

int main(){
  init_grid(4);

  foreach() {
    u.x[] = x;
    u.y[] = x;
  }

  /* This loop fails to trigger the automatic BC, since the value of u.y in 
    the ghost cell is incorrect. */ 
  double dd=L0/N;
  foreach(){
    if (fabs(x - 0.5) < dd && fabs(y) < dd) {
      printf("loop1, pid:%d x:%g y: %g uy:%g %g\n", pid(), x, y, 
      u.y[], u.y[0, -1]);
    }
  }

  printf("\n");
  /* However, this loop works. In short, only the x component of u triggers
    the automatic BC. */ 
  foreach() {
    if (fabs(x - 0.5) < dd && fabs(y) < dd) {
      printf("loop2, pid:%d x:%g y: %g uy:%g %g\n", pid(), x, y, 
      u.y[], u.y[0, -1]);
      printf("loop2, pid:%d x:%g y: %g ux:%g %g\n", pid(), x, y, 
      u.x[], u.x[0, -1]);
    }
  }
  
  return 0;
}