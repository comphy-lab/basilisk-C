/**
# Height function
*/
#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"

scalar f[], * interfaces = {f};

vector h[];

int main()
{
  //  size (5);

    f.height = h;

   run();
}

event init (t = 0)
{
  fraction (f, - (sq(x) + sq(y) - sq(0.5)));
}

event print (t=end)
{
  static FILE * fp = fopen ("facets", "w");
  output_facets (f, fp);
  output_cells(stdout);
  
   foreach()
     {
       if (f[] > 1e-6 && f[] < 1. - 1e-6) {
	 fprintf(stderr,"%lf %lf %lf %lf %f\n",x,y,f[],
                 y + height(h.y[]),height(h.y[]));
       }
     }
}

/**
~~~gnuplot 
plot [0:0.6][0:0.6] "out" u 1:2 w l t "grid","log" u 1:2 w p t "center","facets" u 1:2 w l t "facets"
~~~

 */