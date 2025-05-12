/** 
Vertical (y-direction) profile function for scalar fields on quad and octree
This function is an inefficient way of extracting profiles of scalarfields in the y-direction (for me, x and z are homogeneous directions). It probes every possible y-level, and calculates the volume(!) averaged value at that height, irrespective of the amount of gridcells (can be an average over only one grid cell). Therfore, the amount of gridcells and volume at that height is also noted. Note this also works on quadtrees, then read 'surface' instead of 'volume'  

the function requires a scalarfield, the lower and upperbounds of the height of your desired profile, and the maximum resolution (in levels) 

It now prints to the terminal 
*/ 

void profile(scalar A,double ym,double h, int max)
{
  printf("y\tVAR\tvolume\t#Cells\n");
  
  for (double yp = ym; yp<h; yp+=(L0/pow(2,max+1)))
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
note: This method is generally incorrect. The best practice would be to interpolate all field values on a equidistant grid with the maximum resolution and then simply extraxt profiles from this regular grid. The horizontal avereging is now done with volume (in 3D) weights and this is rather arbitrarry. one might prefer to take surface weights (That would represent a slice of volume). these three methods will proreduce three different results. The extend of the differences depends on the octree structure and the distribution of the scalar field. Furthermore, some levels may only contain a few gridcells, according to statistics, this may cause the profiles to show some 'noise'.  
*/