/**
# Profiles in the y-direction on tree-grids.
This function calculates vertical (y-direction) profiles of a list of scalars. The output it written in a ".dat" file. 

As an unput the function requires a list of scalars, the minimum and maximum height between which it will calculate the profiles, the level of the resolution $\Delta_i = \mathrm{L0}/2^{\text{max}}$ at which the profiles should be calculated and a double specifier, e.g. the physical time.

## Method
Starting from $y = ym + \Delta_{i}/2 $ the solution is evaluated on a horizontal slice by interpolation on a regular grid with $\Delta_i$ resolution. The average of this height is calculated and printed for each scalar in the list. After each interation the height is increased by $\Delta_i$ until the maximum height i reached.

This function is compatible with MPI domain decomposition, however it is rather inefficient. It works best on a $\mathrm{L0} \times \mathrm{L0}$ domain. By default it writes the output in a "data" folder that is not automatically created. 
*/
void profile(scalar * list,double ym,double h, int max,double t)
{
  char names[100];
  sprintf(names,"./data/profilest=%gl=%d.dat",t,max);
  FILE * fpprof = fopen(names,"w");
  int mm = 0;
  int nn = (pow(2,max));
  int len = list_len(list);
  double stp = L0/nn;
  fprintf(fpprof,"y\t");
  for(scalar s in list)
    {      
      fprintf(fpprof,"%s\t",s.name);
      mm++;
    }
  fprintf(fpprof,"\n");
  for (int l = 0;l < nn/(L0/(h-ym)); l++)
    {
    double yp = stp * l +ym +stp/2;
    double ** field = matrix_new (nn, nn, len*sizeof(double));
    for (int i = 0; i < nn; i++)
    {
      double xp = stp*i + X0 + stp/2.;
      for (int j = 0; j < nn; j++) 
      {
        double zp = stp*j + Y0 + stp/2.;
        Point point = locate (xp,yp,zp);
	      int k = 0;
	      for (scalar s in list)
                field[i][len*j + k++] = point.level >= 0 ? interpolate (s, xp, yp,zp)  : nodata;
      }
    }
    if (pid() == 0)
    { // master
      @if _MPI
        MPI_Reduce (MPI_IN_PLACE, field[0], len*nn*nn, MPI_DOUBLE, MPI_MIN, 0,
                    MPI_COMM_WORLD);
	  @endif
	    int k = 0;
      for (scalar s in list)
      {
        double sum = 0;
        for(int i = 0 ; i <nn ; i++)
        {
          for(int j=0;j<nn ; j++)
          {
		      sum+=field[i][j*len+k];   
          }
        }
        if (k ==0)
          fprintf(fpprof,"%g\t%g",yp,sum/(nn*nn));
        if (k > 0)
          fprintf(fpprof,"\t%g",sum/(nn*nn));
        k++;
        if (k == len)
        {
          fprintf(fpprof,"\n");
        }
      }
    }
    @if _MPI
      else // slave
      MPI_Reduce (field[0], NULL, len*nn*nn, MPI_DOUBLE, MPI_MIN, 0,
                  MPI_COMM_WORLD);
      @endif
      matrix_free (field);
  }
  fclose(fpprof);
}