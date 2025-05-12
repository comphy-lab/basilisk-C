/**
# Packing geometry for droplets in 2D 

   This function arranges droplets/bubbles in either hexagonal packing
   or square packing based on initially assigned pk value (2 and 1
   respectively).

   This works with periodic boundaries on left and right. For other
   boundaries, one should change the values in array **ar** by adding
   the shift amount to all of the elements.

   One can also use rotated arrangement of a hexagonal packing. This
   code also facilitates one to introduce a pertubation to teh regular
   arrangement.

## User interface
   Copy and paste this in your main code before adding this file
   and choose the appropriate values for these variables:
*/
/*

#define D 1. //Diameter of the each droplet
#define LEVEL 8
#define NUM 10 // Number of droplets in x direction 
#define NUMY 10 // Number of droplets in y direction
#define DIST (0.2 + D)
#define PK 2 // Hexagonal packing is 2 and Square packing is 1
#define I0 1.25*D // Distance between bottom to the center of bottom
		// most droplets
#define theta_pac 0. //0 for 0 degrees rotaion, 0.523599 for 30
			   //degrees rotation
#define rand_magn 5 //5 represents max(0.05*D) amplitude random
		    //perturbation, for no perturbation set it to 0

*/

/**
  Size of the domain must be declared as equal to **LENGTH**
  (size(LENGTH)). Total number of droplets will be **NUM\*NUMY**.
Calls in the main function should be as follows:
/* 
size(LENGTH);
  init_grid (1 << LEVEL);
  rval(rval1, rval2, NUM, NUMY, rand_magn);
  fraction (f, packing_geometry(x,y));
*/
*/

#if !(PK - 2)
#define XDIST DIST*2.*cos(1.0472 - theta_pac)
#define YDIST DIST*sin(1.0472 - theta_pac)
#elif !(PK - 1)
#define XDIST DIST
#define YDIST DIST
#endif
#define LENGTH NUM*XDIST


//double ar[] = {0, 0.5*(1 + CELL_GAP(GAP)/D)};
//double ar[] = {0, 0.5*(1 + CELL_GAP(GAP)/D) + (CELL_GAP(0.5*GAP)/D)};
double ar[] = {0, 0.5*XDIST};
double rval1[(NUM+1)*NUMY];
double rval2[(NUM+1)*NUMY];
void rval(double * rval1, double * rval2, int numx, int numy, int ravan)
{
  for(int i = 0 ; i < numy ; i++)
    {
      for(int j = 0 ; j <= numx ; j++)
	{
	  if(!ravan)
	    {
	      rval1[i*(numx+1) + j] = 0;
	      rval2[i*(numx+1) + j] = 0;
	    }
	  else
	    {
	      rval2[i*(numx+1) + j] = 0.1*((rand() % (2*ravan)) - ravan);
	      if(!i)
		{
		  rval1[i*(numx+1) + j] = 0;
		  rval2[i*(numx+1) + j] = 0;
		}
	      if(i%PK == 0 && (j == 0 || j == numx))
		{ rval1[i*(numx+1) + j] = 0;
		  rval2[i*(numx+1)] = rval2[i*(numx+1) + numx];
		}
	      else
		{
		  rval1[i*(numx+1) + j] = 0.1*((rand() % (2*ravan)) - ravan);
		}
	    }
	}
    }
}

double packing_geometry(double x, double y)
{
  int num = NUM;
  double maxi = 0, max = 0;
  for(int i = 0 ; i < NUMY ; i++)
    {
      if((i%PK == 0) ? (num = NUM) : (num = NUM - 1));
      for(int j = 0 ; j <= num ; j++)
	{
	  double rval11 = rval1[i*(NUM+1) + j];
	  double rval22 = rval2[i*(NUM+1) + j];
	  maxi =  -(sq(x - (ar[i%PK]*D + j*XDIST + rval11*0.1*D)) +
		    sq(y - (I0 + i*YDIST + rval22*0.1*D)) - sq(0.5*D));
	  if((!i && !j) || (maxi >= max))
	    max = maxi;
	}
    }
  return max;
}
