#include "grid/multigrid3D.h"
#include "utils.h"

#define LEVEL 7

scalar f[];

int main()
{
  init_grid(1 << LEVEL);
  foreach()
    f[] = 1.0;
  dump("test",{f});
  fprintf(ferr,"Dump complete.\n");
  return 1;
}
