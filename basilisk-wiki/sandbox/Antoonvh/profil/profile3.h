/**
# Profiles in the y-direction on tree-grids.
This function calculates vertical (y-direction) profiles of a list of scalars. The output it written in a ".dat" file. This function is a sligth improvement over [an earlier version](profile2.h). 

As an unput the function requires a list of scalars, the minimum and maximum height between which it will calculate the profiles, the level of the resolution $\Delta_i = \mathrm{L0}/2^{\text{max}}$ at which the profiles should be calculated and a double specifier, e.g. the physical time.

## Method
Starting from $y = ym + \Delta_{i}/2 $ the solution is interpolated and summed-up over a horizontal slice on a regular grid with $\Delta_i$ resolution. This average is printed for each scalar in the list. After each interation the height is increased by $\Delta_i$ until the maximum height is reached.

Since typically most of the time is spend in the interpolation routine I consider this function to be parallelized. It works best on a $\mathrm{L0} \times \mathrm{L0}$ domain. By default it writes the output in a "data" folder that is not automatically created. 
*/

void profile(scalar * list,double ym,double h, int max,double t)
{
  char names[100];
  sprintf(names,"./data/profilest=%gl=%d.dat",t,max);
  FILE * fpprof = fopen(names,"w");
  int nn = (pow(2,max));
  int len = list_len(list);
  double aver[len];
  double stp = L0/nn;
  // Header
  fprintf(fpprof,"y\t");
  for(scalar s in list)
    fprintf(fpprof,"%s\t",s.name);
  fprintf(fpprof,"\n");
  
  for (int l = 0;l < nn/(L0/(h-ym)); l++)
    {
      double yp = stp * l +ym +stp/2;
      memset(aver, 0., len * sizeof(double));
      for (int i = 0; i < nn; i++)
	{
	  double xp = stp*i + X0 + stp/2.;
	  for (int j = 0; j < nn; j++) 
	    {
	      double zp = stp*j + Z0 + stp/2.;
	      Point point = locate (xp,yp,zp);
	      int k = 0;
	      for (scalar s in list)
		{
		  aver[ k++ ] += point.level >= 0 ? interpolate (s, xp, yp,zp)  : 0;
		}
	    }
	}
      
      // MPI reduction
      if (pid() == 0)
	{ // master
	  @if _MPI
	    MPI_Reduce (MPI_IN_PLACE, aver[0], len, MPI_DOUBLE, MPI_SUM, 0,
			MPI_COMM_WORLD);
	  @endif
	    int k = 0;
	  for (scalar s in list)
	    {
	      if (k ==0)
		fprintf(fpprof,"%g\t%g",yp,aver[k]/(nn*nn));
	      if (k > 0)
		fprintf(fpprof,"\t%g",aver[k]/(nn*nn));
	      k++;
	      if (k == len)
		{
		  fprintf(fpprof,"\n");
		}
	    }
	}
      @if _MPI
      else // slave
	MPI_Reduce (aver[0], NULL, len, MPI_DOUBLE, MPI_SUM, 0,
		    MPI_COMM_WORLD);
      @endif
	}
  fclose(fpprof);
}
/**
This method has been updated to be "grid-adaptive". [See here](profile4.h) 
*/

