/**
# A Grid Adaptive Profiling Function for Oct-tree Grids.
This function prints vertical (y-direction) profiles of a list of scalars. It incorporates a grid adaptive method of an [earlier version](profile3.h). 

The adaptivity works as follows: For a given height and a given x-coordinate the solution will be probed in the z-dicection (i.e. a z-track). The distance between probe locations is determined by the local grid box size. The average grid box size probed on this z-track at that will be the step taken in the x-direction for the next z-track. Once the full x-z plane has been probed the average (over the probe locations) is calculated and printed. The next y-height is determined by an (estimated) average grid box size of the current x-z plane. 

This function can be much faster than earlier profiling functions. Most notably, the effort required for this function is now approximately proportinal to the number of grid cells. Furthermore, an additional argument ($rf$), the "reduction factor" can be given. This value is a multiplier for the step-size taken in all directions. E.g. $rf=1.71$ will result in $1.71^3\approx 5$ times fewer probe locations, the default value is 1.    

Note that the asymetry between x and z direction could potentially cause unwanted effects under some conditions. Also the asymetry between locations close the the left and front boundary could be a topic of concern.
*/
struct prof {
  scalar * list; // List of scalars
  double ym;     // Start height for profile
  double h;      // Max height of profile (h > ym )
  double t;      // specefier for output file name (e.g. time) 
  int max;       // Maximum level of refinement
  double rf;     // Reduction factor
};

void profile(struct prof p){
  char names[100];
  sprintf(names,"./profilest=%g.dat",p.t);
  FILE * fpprof = fopen(names,"w");
  if (!p.rf)
    p.rf=1;
  if (!p.max)
    p.max=depth();
  int ny,nz;
  int len = list_len(p.list);
  double aver[len];
  double stp = L0/pow(2,p.max);
  // Header
  fprintf(fpprof,"y\t");
  for(scalar s in p.list)
    fprintf(fpprof,"%s\t",s.name);
  fprintf(fpprof,"\n");
  double yp = p.ym + stp/2;
  
  while(yp<=p.h){
    ny=0;
    double xp = X0 + stp/2;
    memset(aver, 0., len * sizeof(double));
    while(xp<(X0+L0)){
      nz=0;
      double zp = Z0 + stp/2.;
      while(zp<(Z0+L0)){
	Point point = locate (xp,yp,zp);
	int k = 0;
	for (scalar s in p.list)
	  aver[ k++ ] += point.level >= 0 ? interpolate (s, xp, yp,zp)  : 0;
	nz++;
	ny++;
	zp+=p.rf*Delta;
      }
      xp+=L0/nz;
    }
    // MPI reduction
    if (pid() == 0){ // master
      @if _MPI
	MPI_Reduce (MPI_IN_PLACE, aver[0], len, MPI_DOUBLE, MPI_SUM, 0,
		    MPI_COMM_WORLD);
      @endif
	int k = 0;
      for (scalar s in p.list){
	if (k ==0)
	  fprintf(fpprof,"%g\t%g",yp,aver[k]/(ny));
	if (k > 0)
	  fprintf(fpprof,"\t%g",aver[k]/(ny));
	k++;
	if (k == len)
	  fprintf(fpprof,"\n");
      }
    }
    @if _MPI
    else // slave
      MPI_Reduce (aver[0], NULL, len, MPI_DOUBLE, MPI_SUM, 0,
		  MPI_COMM_WORLD);
    @endif
      yp+=L0/sqrt(ny);
  }
  fclose(fpprof);
}
/**
## Do not use
The horizontal averages are biased towards the values that are at a higher level of refinement. This in not a nice feature...
You could follow the link for a more consistent profiling function [i.e. here](profil5.h).