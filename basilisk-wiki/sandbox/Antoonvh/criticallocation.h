/**
Find a location and grid size that is critical with respect to the CFL criterion.
*/

void find_a_critical_location(double xcrit,double ycrit,double zcrit,double lev)
{
  double xyzdt[5];
  double dtmax = 100.0*DT;
  int n=0;

  @if _MPI
    double xyzdtglob[5*npe()];
  @endif
    
  foreach_face(reduction(min:dtmax)){
    if (u.x[] != 0.) {
      double dt = Delta*cm[]/fabs(u.x[]);
      if (dt < dtmax){
        xyzdt[0]=x;
        xyzdt[1]=y;
        xyzdt[2]=z;
        xyzdt[3]=dt;
        xyzdt[4]=(double)level;
        dtmax = dt;
      }
    }
  }
  
  @if _MPI    
    MPI_Gather(&xyzdt,5,MPI_DOUBLE,xyzdtglob,5,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for (int m = 0;m<npe();m++){
      if (xyzdtglob[5*m+3]==dtmax){
        n=m; 
        // fprintf(ferr,"Critical cell proc. ID = %d\n",n)
        xcrit=xyzdtglob[5*n];
        ycrit=xyzdtglob[(5*n)+1];
        zcrit=xyzdtglob[(5*n)+2];
        lev=xyzdtglob[(5*n)+4];
        m+=npe();
      }
    }
  @endif
    
  @if !_MPI  
    xcrit=xyzdt[0];
    ycrit=xyzdt[1];
    zcrit=xyzdt[2];
    lev=xyzdt[4];
  @endif
}