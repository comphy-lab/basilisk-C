#define GHOSTS 2
#include "utils.h"

#define OFFSET(l)  (((1 << ((2*l)+2)) - 1)/3)
#define _O (-GHOSTS)
#define C_IND(i,j,l) ((i + _O) + (1 << l)*(j + _O))
#define INDEX (OFFSET(level - 1) + C_IND(point.i, point.j, level)) 

scalar b[];

int main() {
  foreach_dimension()
    periodic(left);
  for (int t = 1; t<= 5; t += 2) {
    char fname[99];
    sprintf(fname, "%d-dump%d", t, 9);
    restore (fname);
    b.restriction = restriction_vertex;
    restriction ({b});
    printf ("%d\n", depth());
    sprintf(fname, "%d-risedata", t);
    FILE * fp = fopen (fname, "wb");
    foreach_cell() {
      fseek (fp, INDEX*sizeof(double), SEEK_SET);
      fwrite (&b[], 1, sizeof(double), fp);
    }
    fclose (fp);
  }
}
