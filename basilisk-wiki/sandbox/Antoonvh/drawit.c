#include "view.h"

#define OFFSET(l)  (((1 << ((2*l)+2)) - 1)/3)
#define _O (-GHOSTS)
#define C_IND(i,j,l) ((i + _O) + (1 << l)*(j + _O))
#define INDEX (OFFSET(level - 1) + C_IND(point.i, point.j, level)) 

scalar b1[], b3[], b5[], bt[];

int main() {
  L0 = 30;
  X0 = Y0 = -L0/2;
  init_grid (1 << 9);
  FILE * fp1 = fopen ("1-risedata", "rb");
  FILE * fp3 = fopen ("3-risedata", "rb");
  FILE * fp5 = fopen ("5-risedata", "rb");

  foreach() {
    fseek (fp1, INDEX*sizeof(double), SEEK_SET);
    fread(&b1[], 1, sizeof(double), fp1);
    fseek (fp3, INDEX*sizeof(double), SEEK_SET);
    fread(&b3[], 1, sizeof(double), fp3);
    fseek (fp5, INDEX*sizeof(double), SEEK_SET);
    fread(&b5[], 1, sizeof(double), fp5);
  }
  view (fov = 12, ty = -1/30., samples = 4, width = 600, height = 750);
  isoline ("b1", .5, lw = 10, lc = {31./255., 119./255., 180./255.});
  isoline ("b3", .5, lw = 10, lc = {255/255., 127./255., 14/255.});
  isoline ("b5", .5, lw = 10, lc = {44./255., 160./255., 44/255.});
  draw_string ("(a)", lw = 2, pos = 1, size = 15);
  foreach()
    bt[] = b5[] - b1[];
  squares("bt", map = cool_warm, min = -.8, max = .8);
  save ("Drawing.png");
  dump ("bs");
}
