#include "grid/multigrid3D.h"
#include "utils.h"

#define LEVEL 7

scalar f[];

int main()
{
  init_grid(1 << LEVEL);
  restore("test");
  return 1;
}
