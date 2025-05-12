#include "grid/multigrid.h"
#include "utils.h"
#include "fractions.h"

scalar f[];

int main()
{
  init_grid(1);
  restore ("dump");
  output_facets (f);
}
