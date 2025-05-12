#include "utils.h"

#define OFFSET(l)  (((1 << ((2*l)+2)) - 1)/3)
#define _O (-GHOSTS)
#define C_IND(i,j,l) ((i + _O) + (1 << l)*(j + _O))
#define INDEX (OFFSET(level - 1) + C_IND(point.i, point.j, level))

scalar b[];

int main() {
  foreach_dimension()
    periodic(left);
  int lm = 9;
  for (int l = 5; l < lm; l++) {
    for (int t = 1; t<= 5; t += 2) {
      char dname[99], fname[99];
      sprintf (dname, "%d-dump%d", t, l);
      sprintf (fname, "%d-risedata", t);
      restore (dname);
      FILE * fp = fopen (fname, "rb");
      double e = 0, em = -1;
      if (t == 1)
	printf ("%d ", 1 << l);
      foreach(reduction(+:e) reduction(max:em)) {
	fseek (fp, INDEX*sizeof(double), SEEK_SET);
	double bi;
	fread (&bi, 1, sizeof(double), fp);
	double el = fabs(b[] - bi);
	e += dv()*el;
	if (el > em)
	  em = el;
	b[] = bi;
      }
      printf ("%g %g ", e, em);
      fclose (fp);
    }
    printf ("\n");
  }
}
