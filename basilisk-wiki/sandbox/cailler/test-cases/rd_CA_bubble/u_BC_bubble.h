//                         /* VELOCITY VECTOR BC */

// Projection onto e_z:
double proj_ez (double fR, double fth, double th){
  return (
    fR*cos(th) - fth*sin(th)
  );
}

// Projection onto e_r:
double proj_er (double fR, double fth, double th){
  return (
    fR*sin(th) + fth*cos(th)
  );
}


#if AXI
  // AXI SOURCE:
  // ===========

  double ur_A (double x, double y){
    // normal component
    double R = sqrt(sq(x) + sq(y));
    double th = atan2(y,x); 
    double uR = sq(0.5*SIZE)/sq(R) ;
    double uth = 0. ;
    double ur = proj_er(uR, uth, th);
    
    return (ur);
  }

  double uz_A (double x, double y){
    // tangential component
    double R = sqrt(sq(x) + sq(y));
    double th = atan2(y,x); 
    double uR = sq(0.5*SIZE)/sq(R) ;
    double uth = 0. ;
    double uz = proj_ez(uR, uth, th);
    
    return (uz);
  }

  // analytical try for pressure: 28/02/2024
  double p_axi_ext (double x, double y){
    double R = sqrt(sq(x) + sq(y));
    return ( - pow(0.5*SIZE/R, 4.) );
  }

  double p_axi_bubble (double x, double y){
    double R = sqrt(sq(x) + sq(y));
    double p1 = 2./R;
    double p2 = - pow(0.5*SIZE/R, 4.);
    return (p1 + p2);
  }

#else // !AXI --> 2D
// 2D SOURCE:
// ==========

  double ur_A (double x, double y){
    // normal component
    double R = sqrt(sq(x) + sq(y));
    double th = atan2(y,x); 
    double uR = (0.5*SIZE)/R ;
    double uth = 0. ;
    double ur = proj_er(uR, uth, th);
    
    return (ur);
  }

  double uz_A (double x, double y){
    // tangential component
    double R = sqrt(sq(x) + sq(y));
    double th = atan2(y,x); 
    double uR = (0.5*SIZE)/R ;
    double uth = 0. ;
    double uz = proj_ez(uR, uth, th);
    
    return (uz);
  }
#endif // !AXI



