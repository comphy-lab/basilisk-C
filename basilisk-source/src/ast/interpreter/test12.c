/**
Checks that minmod2() is always undefined. */

#include "utils.h"

int main()
{
  init_grid (1);
  scalar s[];
  {
    interpreter_verbosity (2);
    foreach() {
      if (s[-1]) {
	double a = minmod2 (s[-1], s[], s[1]);
	if (a);
      }
    }
  }
}
