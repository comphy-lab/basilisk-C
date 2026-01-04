/** 
This simple code shows a problem with the usage of foreach_boundary when periodic boundary condition is applied in only one direction.
*/
#define BGHOSTS 1
#include "grid/quadtree.h"

int main()
{
  L0 = 8;
  periodic(top);
  init_grid(1 << 2);

  FILE * fp = fopen ("out.txt", "w");
  fprintf(fp, "%s \n","----------boundary iterator for left boundary-------------");
  foreach_boundary(left)
  {
  	fprintf(fp, "x:%g\ty:%g \n",x,y);
  }
  fprintf(fp, "%s\n","-----------actual left boundary------------" );
  foreach_face(x)
  {
  	if(0 == x)
  		fprintf(fp, "x:%g\ty:%g \n",x,y);
  }
  fclose (fp);
}
/**
This is the output:

![](foreach_bnd/out.txt){width="100%"}

If the foreach_boundary() iterator works correctly, the position information between splitlines should be identical. However, the foreach_boundary() iterator trys to access the ghost cell instead of the actual cell in the domain. Clearly an offset exists in the foreach_boundary() iterator and its value is related to a macro named BGHOSTS. */
