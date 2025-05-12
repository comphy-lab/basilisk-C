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

  printf("%s \n","----------boundary iterator for left boundary-------------");
  foreach_boundary(left)
  {
  	printf("x:%g\ty:%g \n",x,y);
  }
  printf("%s\n","-----------actual left boundary------------" );
  foreach_face(x)
  {
  	if(0 == x)
  		printf("x:%g\ty:%g \n",x,y);
  }
}
/**
Running the above piece of code gives
*/
----------boundary iterator for left boundary-------------
x:0    y:-1
x:0    y:1
x:0    y:3
x:0    y:5
-----------actual left boundary------------
x:0    y:1
x:0    y:3
x:0    y:5
x:0    y:7 
/**
If the foreach_boundary() iterator works correctly, the position information between splitlines should be identical. However, the foreach_boundary() iterator trys to access the ghost cell instead of the actual cell in the domain. Clearly an offset exists in the foreach_boundary() iterator and its value is related to a macro named BGHOSTS. */
