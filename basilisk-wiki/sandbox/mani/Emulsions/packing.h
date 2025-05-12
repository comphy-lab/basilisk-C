
/**
# Packing geometry for droplets in 2D 

   This works with periodic boundaries on left and right and origin
   centred at the default left bottom corner. For other boundaries,
   one should change the values in array **ar** by adding the shift
   amount to all of the elements. For example ar[] = {1, 2, 1.5,
   3.5}. This can also be used to change the arrangement of droplets
   in the x-direction.

   Define **D**: Diameter of the droplet before adding this file.

   Change other parameters **(i0, num, etc.)** in **main()** function

   Total number of droplets will be **num*numy**. Domain size must be
   twice the num: size(2*num)
*/

static double cell_gap(int gap, double length) 
{
return gap*length/(1 << LEVEL);
}

double i0 = 1.*D; //initial height
int num = 4; // number of droplets in x direction
int numy = 3; //number of droplet layers in y direction
/** Distance between one layer of droplets to another can be
    controlled by iy0 and ix0. This can be changed according to once
    preference. */

double iy0 = 2.*D;
double ix0 = 2.*D;
int pk = 4;

/** Arrangement of droplets in the x-direction for different
    layers. The third and fourth element of **ar** can be slightly
    modified based on the level of refinement so that the interface is
    properly initiated close to the boundary. */
double ar[] = {0, 1, 0.5, 1.5};
double packing_geometry(double x, double y)
{
  int Num = num;
  double maxi, max;
  for(int i = 0 ; i < numy ; i++)
    {
      if((i%pk == 0) ? (Num = num) : (Num = num - 1));
      for(int j = 0 ; j <= Num ; j++)
	{
	  maxi =  -(sq(x - (ar[i%pk]*D + j*ix0)) +
		    sq(y - (i0 + i*iy0)) - sq(0.5*D));
	  if((!i && !j) || (maxi >= max))
	    max = maxi;
	}
    }
  return max;
}




