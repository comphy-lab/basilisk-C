/**
# My Fifth profiling function 
[After](profile1.h) [some](profile2.h) [problematic](profile3.h) [earlier](profile4.h) versions, i am happy to present a new profiling function.

This function slices the octree at various heights (y-coordinates) and interpolates the values of the scalar fields in the list. The function also works for non cubic domains. 
A drawback is that all cells are checked for each height over and over again. Atleast this drawback scales with the total number of gridcells. 

Since the distance between the probe heights is a certain input factor times the average grid cell height at a given height, this function can be called grid-adaptive. Furthermore, since most of the calculations are done in the foreach() loop, its also parallelized.  

syntax,

- list is a list of scalars that you want to have profiles of. 
- ym is the starting height of the profile
- h is the maximum height of the profile (h>ym)
- t is a double that is printed in the file name (e.g. time)
- rf is the reduction factor of the probed locations (default value is 1)

e.g.
profile({A,B,C}, 0.1, 1000., t, 1.5);
*/

struct prof {
  scalar * list;
  double ym;
  double h;
  double t;
  double rf;
};

void profile(struct prof p){
  char names[100];
  sprintf(names,"profilest=%g.dat",p.t);
  FILE * fpprof = fopen(names,"w");
  if (!p.rf)
    p.rf=1;
  int len = list_len(p.list);
  
  // Header
  fprintf(fpprof,"y\t");
  for(scalar s in p.list)
    fprintf(fpprof,"%s\t",s.name);
  fprintf(fpprof,"\n");
  double Dolt =1;
  for (double yp=p.ym;yp<=p.h;yp+=Dolt){
    double aver[len];
    for (int i=0;i<len;i++)
      aver[i]=0.;
    int m=0;
    double a=0;
    foreach(reduction(+:a) reduction(+:m)){
      if ((fabs(y-yp)<=(Delta/2))){
        m++;
        a+=sq(Delta);
        int k = 0;
        for (scalar s in p.list)
          aver[k++] += point.level >= 0 ? interpolate (s, x, yp,z)*sq(Delta)  : 0; // I think the foreach loop may automatically ensure that point.level>=0 
      }
    }
    // MPI reduction
    if (pid() == 0){ // master
      @if _MPI
        MPI_Reduce (MPI_IN_PLACE, &aver[0], len, MPI_DOUBLE, MPI_SUM, 0,
                    MPI_COMM_WORLD);
      
      @endif
        int k = 0;
      for (scalar s in p.list){
	if (k ==0)
          fprintf(fpprof,"%g\t%g",yp,aver[k]/a);
	if (k > 0)
          fprintf(fpprof,"\t%g",aver[k]/a);
        k++;
        if (k == len)
          fprintf(fpprof,"\n");
      }
    }
    @if _MPI
    else // slave
      MPI_Reduce (&aver[0], NULL, len, MPI_DOUBLE, MPI_SUM, 0,
                  MPI_COMM_WORLD);
    @endif
      // Increment of height
      Dolt=p.rf*(sqrt(a))/sqrt(m);
  }
  fclose(fpprof);
}
