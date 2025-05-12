/** 
#Vertical (y-direction) profile function for scalar fields on quad and octree

This function can provide profiles of scalarfields in the y-direction. It probes every possible y-level, and calculates the volume(!) averaged value at that height, irrespective of the amount of gridcells (can be an average over one grid cell only). Therfore, the amount of gridcells and volume at that height is also noted. Note this also works on quadtrees, then read 'surface' instead of 'volume'.  

the function requires a scalarfield, the and the upperbound for the height of your desired profile. 

The function prints the result to the terminal and is MPI domain decomposition friendly.  

I do not suggest to use this function as it only works well under special conditions. Curretly [this](profil5.h) is my best effort. 
*/ 

void profile(scalar A,double hmax)
{
  int max = depth();
  printf("y\tVAR\tvolume\t#Cells\n");
  for (double yp = Y0; yp<h; yp+=(L0/pow(2,max+1)))
  {
    double vol = 0;
    double s = 0;
    int n = 0; 
    foreach(reduction(+:vol) reduction(+:s) reduction(+:n))
    {	
      if (yp == y)
      {
        vol +=dv();
        s+=(dv()*A[]);
        n+=1;
      }
    }
    if (n>0)
      printf("%g\t%g\t%g\t%d\n",yp,s/vol,vol,n);
  }
}
/**
Remember not to use it, the main issue is that at each height there exist only cells at a single level. 
*/