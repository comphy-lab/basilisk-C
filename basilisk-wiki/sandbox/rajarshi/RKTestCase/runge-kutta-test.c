/**
#RUNGE KUTTA TEST CASE DEMONSTRATION.
*/

#include "grid/multigrid.h"
#define dimension 2

#include "utils.h"
#include "../Header_File/runge-kutta-new.h"

static void du (scalar * ul, double t, scalar * kl){
  vector * ulf = NULL;
  vector * klf = NULL;
  scalar * ulc = NULL;
  scalar * klc = NULL;

  scalar f,g;
  for (f,g in ul,kl){
     if(f.face){
        ulf = vectors_add(ulf,f.v);
        klf = vectors_add(klf,g.v);
     }
     else{
        ulc = list_add(ulc,f);
        klc = list_add(klc,g);
     }
  }

  if(ulc!=NULL){
    foreach(){
      for (f,g in ulc,klc) 
         g[] = t*f[];  
    }
  }
  if(ulf!=NULL){
    vector h,i;
    foreach_face(){
      for(h,i in ulf,klf)
         i.x[] = t*h.x[]; 
     }      
  }
 boundary(kl);
 free(ulf);
 free(ulc);
 free(klf);
 free(klc);
}


int main()
{
  init_grid (1);

  for (int order = 1; order <= 4; order *= 2)
    for (double dt = 1e-2; dt <= 8e-2; dt *= 2) {
      scalar a[];
      vector b[];
      face vector c[];
      foreach(){
	a[] = 1.;
        foreach_dimension()
          b.x[] = 2.;
      }
      foreach_face()
          c.x[] = 3.;

      double ea_max = 0.;
      double eb_max = 0.;
      double ec_max = 0.;

      for (t = 0; t <= 2.; t += dt) {
	foreach() {
            double e = fabs (a[] - exp(t*t/2.));
	    if (e > ea_max)
	        ea_max = e;
	  
            foreach_dimension(){
                double f = fabs(b.x[] - 2.*exp(t*t/2.));
                if (f > eb_max)
                   eb_max = f; 
            }
        }

        foreach_face() {
            double g = fabs(c.x[] - 3.*exp(t*t/2.));          
            if (g > ec_max)
               ec_max = g;
        }  

	runge_kutta ({a,b,c}, t, dt, du, order);
      }
      fprintf (stderr, "%g \t\t %g \t\t %g \t\t %g \t\t %d\n", dt, ea_max, eb_max, ec_max, order);
    }
}

/**
~~~gnuplot Error convergence for different orders
set xlabel 'dt'
set ylabel 'error'
set logscale
set key outside
plot "< grep '1$' log" u 1:2 t 'scalar-O1', 15.*x t '15 dt',  \
     "< grep '1$' log" u 1:3 t 'vector-O1', "< grep '1$' log" u 1:4 t 'face vector-O1',   \
     "< grep '2$' log" u 1:2 t 'scalar-O2', 4.*x*x t '4 dt^2', \
     "< grep '2$' log" u 1:3 t 'vector-O2', "< grep '2$' log" u 1:4 t 'face vector-O2',   \
     "< grep '4$' log" u 1:2 t 'scalar-O4', x**4/2. t 'dt^4/2', \
     "< grep '4$' log" u 1:3 t 'vector-O4', "< grep '4$' log" u 1:4 t 'face vector-O4'

~~~
*/
