#define dimension 2
#define GRIDNAME "MyQT"
typedef double real;
#include "mytree.h"
#include "prolongate_halo_flux.h"
#include "mytree-utils.h"
void myQT_methods() {
  tree_methods();
  boundary_face = dummy;
  //boundary_flux = halo_flux_prolongate;
  boundary_level = tree_boundary_level;
}
