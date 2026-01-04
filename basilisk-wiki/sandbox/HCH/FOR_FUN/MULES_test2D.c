/**
This is a bubble moving test case for MULES. 
*/
#include "grid/cartesian.h"
#include "run.h"
#include "fractions.h"
#include "./TVD_nosplit.h"
#include "./TVD_split.h"
#include "./MULES.h"

scalar f1[], f2[], f3[], f4[];
			     
face vector uf[];

int main() {
  origin (-0.5, -0.5);
  L0 = 1.;
  N = 256;

  FILE *output_file = fopen("moving_output", "w");
  if (!output_file) {
    fprintf(stderr, "Error opening output file\n");
    return 1;
  }

  run();

}

#define circle(x,y) (sq(0.1) - sq(x + 0.3) - sq(y + 0.3))

event init(i = 0)
{
  foreach_face(x) uf.x[] = 0.1;
  foreach_face(y) uf.y[] = 0.0;

/**
1:MULES, 2:no split TVD advection, 3: dimensional split TVD advection, 4: pure upwind advection
*/
  fraction (f1, circle(x,y));
  fraction (f2, circle(x,y));
  fraction (f3, circle(x,y));
  fraction (f4, circle(x,y));
}

event velocity (i++, last) {
  dt = dtnext (0.001);
}

event advection (i++, last) {
  scalar * tracer1 = {f1};
  scalar * tracer2 = {f2};
  scalar * tracer3 = {f3};
  scalar * tracer4 = {f4};
  MULES_advection (tracer1, uf, dt);
  TVD_nosplit_advection (tracer2, uf, dt);
  TVD_advection (tracer3, uf, dt, i);
  upwind_advection (tracer4, uf, dt);
}

event dtprint (i += 100, last) {
  fprintf(stderr, "i = %d \t dt = %f\n", i, dt);
}

event output_event (i+=100) {
  double allerr1 = 0., allerr2 = 0., allerr3 = 0.;
  foreach() {
    allerr1 += fabs(f3[] - f1[]);
    allerr2 += fabs(f3[] - f2[]);
    allerr3 += fabs(f1[] - f4[]);
  }

  FILE *output_file = fopen("moving_output", "a");  // Append mode
  if (output_file) {
    fprintf(output_file, "t=%g MULES error is %g, nosplit TVD error is %g, difference between upwind is %g\n", t, allerr1, allerr2, allerr3);
    fclose(output_file);
  } 
}

event display(i+=100)
{
  /**
  We recompute $\lambda$ here since it will be freed after the call of the function
  */
  face vector flux_corr_t[];
  face vector flux_low_t[];
  face vector lambda_t[];

  upwind (f1, uf, flux_low_t);
  corr_flux (f1, uf, flux_corr_t);
  limiter(f1, flux_low_t, flux_corr_t, lambda_t, dt);

  char name0[100];
  sprintf(name0, "lambda-%d", i);
  FILE *lambdaout = fopen(name0, "w");
  foreach_face()
  {
    fprintf(lambdaout, "%f \t %f \t %f\n", x, y, lambda_t.x[]);
  }
  fclose(lambdaout);

  char name1[100];
  sprintf(name1, "MULES-%d", i);
  FILE *face1 = fopen(name1, "w");
  output_facets(f1, face1);
  fclose(face1);

  char name2[100];
  sprintf(name2, "no_split_TVD-%d", i);
  FILE *face2 = fopen(name2, "w");
  output_facets(f2, face2);
  fclose(face2);

  char name3[100];
  sprintf(name3, "TVD-%d", i);
  FILE *face3 = fopen(name3, "w");
  output_facets(f3, face3);
  fclose(face3);

  char name4[100];
  sprintf(name4, "upwind-%d", i);
  FILE *face4 = fopen(name4, "w");
  output_facets(f4, face4);
  fclose(face4);
}

event stop (t = 4);