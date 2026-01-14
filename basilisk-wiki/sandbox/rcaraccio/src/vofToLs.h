/**
Simple conversion of a VOF field to a level-set field.
*/

#include "redistance.h"

void vof_to_ls (scalar f, scalar levelset, int imax = 3)
{
  double deltamin = L0/(1 << grid->maxdepth);
  foreach()
    levelset[] = -(2.*f[] - 1.)*deltamin*0.75;
#if TREE
  restriction({levelset});
#endif
  redistance (levelset, imax);
}
