/**
#Refine Function Bug
*/

#include "grid/bitree.h"
#include "utils.h"

#if dimension == 1
double SolutionCompact (double x, double Delta){
  double xq1, xq2;
  xq1 = x - Delta/(2.*sqrt(3));
  xq2 = x + Delta/(2.*sqrt(3));
  if(sq(x) <=1)
    return ( 0.5 * ( pow((1-sq(xq1)),5) + pow((1-sq(xq2)),5) )  );
  else
    return (0);
}
#endif

int main (){

  FILE *fp;
  L0 = 4;
  origin (-2);
  init_grid(1<<4);
  scalar s[];
 
  fp = fopen("BeforeRefine.dat","w");
  foreach(){
     s[] = SolutionCompact (x,Delta);
     fprintf(fp,"%g %g\n",x,s[]);
  }
  fclose(fp);

  foreach_dimension(){
     s[left] = neumann(0);
     s[right] = neumann(0);
   }
  boundary({s});

#if 0  
  refine((level < 5 ) && (sq(x) <= 0.25) );
#else
  do {									
    int refined;								
    do {									
      refined = 0;							
      tree->refined.n = 0;						
      foreach_leaf()							
	if ((level < 5 ) && (sq(x) <= 0.25)) {
	  refine_cell (point, all, 0, &tree->refined);			
	  refined++;							
	  continue;							
	}									
      mpi_all_reduce (refined, MPI_INT, MPI_SUM);				
      if (refined) {						
	mpi_boundary_refine (all);					
	// mpi_boundary_update (all);					
      }
    } while (refined);
  } while(0);
#endif
  
  fp = fopen("AfterRefine.dat","w");
  foreach_level(4)
    fprintf(fp,"%g %g\n",x,s[]);
  fclose(fp);
}

/**
## Results

~~~gnuplot Coarse cell values - Before and After refine
set key outside
set output 'Before_After_Refine.png'
set xlabel 'X'
set ylabel 'Solution(X)'
set grid
plot 'BeforeRefine.dat' u 1:2 w p t 'Before-Refine', 'AfterRefine.dat' u 1:2 w p t 'After-Refine'
~~~

*/
