/**
More tests of undefined conditions and evaluations of branches. */

int main()
{
  init_grid (1);
  scalar s[];
  foreach() {
    interpreter_verbosity (4);
    coord a = {33};
    if (max(0,s[]) > 0) {
      double q = 1;
      coord u = {0};
      if (q)
	a = u;
      q = a.x;
      if (q);	  
    }
    else {
      a.x; // this must be 33
      a.x = 2;	
    }
    if (a.x > 0) // this must be (unset)
      ;
  }
}
