/**
# Edge Based Interface Tracking solver

We wish to advect an interface identified by marker points that are constrained
to move only along the grid lines. We refer to it as EBIT method and
the following assumptions make it a Semushin's method:
The underlying grid is a 2D square grid, the intersections are at most two
per square edge of the grid and the interpolation
between the marker points is linear.
*/

#include "vof.h"

@define exist(val) (val > 0. && val < 1.)
@define existNeg(val) (val > -1. && val < 0.)
#define myface_value(a,i) ((a[0,i] + a[0,i-1])/2.)

face vector s[], snew[], s_tmp[], s_flag[];
double tTime = 0.;

/**
  Given the 4 coordinates of 2 marker points, it returns the y coordinate
  of the point intersecting the vertical grid line at $x = 1$.
  */
static double intersect(double x1, double y1, double x2, double y2){
  double m = (y2 - y1)/(x2 - x1);
  return (m*1. + (y1 - m*x1));
}

/** Injection and restriction operators more or less work when using adaptivity.
Note that there are some workarounds that will be removed in future versions. */

static void myface_injection (Point point, vector v){ 
  foreach_dimension(){
    
    for (int i=-2; i<=2; i++)  // At the moment, this is needed
      for (int j=-2; j<=2; j++)
        if (fine(v.x,i,j) > 0.)   // Needed?
          v.x[] = 0.;
        
        if (fine(v.x,0,0) > 0.)
          v.x[] = fine(v.x,0,0)/2.;
        
        if (fine(v.x,0,1) > 0.)
          v.x[] = 0.5 + fine(v.x,0,1)/2.;
  }    
}


static inline void myrestriction_face (Point point, scalar s)
{
  myface_injection (point, s.v);
}

foreach_dimension()
  static void myrefine_face_x (Point point, scalar s)
  {
    vector v = s.v;

    if (!is_refined(neighbor(-1)) &&
      (is_local(cell) || is_local(neighbor(-1)))) {
      if (v.x[0] >= 0.5){
        fine(v.x,0,1) = (v.x[0] - 0.5)*2.;
//         printf("v: %g   fine: %g \n",v.x[],fine(v.x,0,1));
      }
      else{
        fine(v.x,0,0) = v.x[0]*2.;
//         printf("v: %g   fine: %g \n",v.x[],fine(v.x,0,0));
      }
      }  
      
      if (!is_refined(neighbor(1)) && neighbor(1).neighbors &&
        (is_local(cell) || is_local(neighbor(1)))) {
        if (v.x[1] >= 0.5){
          fine(v.x,2,1) = (v.x[1] - 0.5)*2.;
//           printf("v: %g   fine: %g \n",v.x[1],fine(v.x,2,1));
        }
        else {
          fine(v.x,2,0) = v.x[1]*2.;
//           printf("v: %g   fine: %g \n",v.x[1],fine(v.x,2,0));
        }
        }
        
        if (is_local(cell)){     //Mid line
          if (v.x[] > 0. && v.x[1] > 0.){  //Both left and right exist and mid is easy
            if (v.x[] + v.x[1] > 1 ){   //Average is > 0.5
              fine(v.x,1,1) = (v.x[] + v.x[1]) - 1.;
              if ((v.x[] - 0.5)*(v.x[1] - 0.5) < 0.){  // One above and one below
                if ((v.x[] - 0.5)>0.){
                  fine(v.y,1,1) = ((0.5 - v.x[])/(v.x[1] - v.x[]) - 0.5)*2.;
//                   printf("a1 %g\n", fine(v.y,1,1));
                }
                else {
                  fine(v.y,0,1) = (0.5 - v.x[])/(v.x[1] - v.x[])*2.;
//                   printf("a2 %g\n", fine(v.y,0,1));
                }
              }
            }
            else {    //Average is <= 0.5
              fine(v.x,1,0) = (v.x[] + v.x[1]);
              if ((v.x[] - 0.5)*(v.x[1] - 0.5) < 0.){  // One above and one below
//                 printf(" %g %g \n",v.x[],v.x[1]+v.x[]);
                if ((v.x[] - 0.5)<0.){
                  fine(v.y,1,1) = ((0.5 - v.x[])/(v.x[1] - v.x[]) - 0.5)*2.;
//                               printf("a1 %g\n", fine(v.y,1,1));
                }
                else {
                  fine(v.y,0,1) = (0.5 - v.x[])/(v.x[1] - v.x[])*2.;
//                               printf("a2 %g\n",fine(v.y,0,1));
                }
              }
            }
          }
          
          else if (exist(v.x[]) && exist(v.y[0,1])){  //Both left and right exist and mid is easy
            double m = (1. - v.x[])/(v.y[0,1] - 0.);
            double val = m*0.5 + v.x[];
//                     printf("val: %g, v.x: %g, v.y: %g m %g %s\n",val,v.x[],v.y[0,1],m,v.x.name);
            if (val > 0.5 && val <= 1.){
              fine(v.x,1,1) = (val - 0.5)*2.;
//                         printf("val1: %g, v.x: %g, v.y: %g fine %g\n",val,v.x[],v.y[0,1],fine(v.x,1,1));
            }
            else if (val >= 0. && val <= 0.5){
              fine(v.x,1,0) = val*2.;
//                         printf("val2: %g, v.x: %g, v.y: %g fine %g %s\n",val,v.x[],v.y[0,1],fine(v.x,1,0),v.x.name);
            }
            double x3 = (0.5 - v.x[])/(m);
            if (x3 >=0 && x3 < 0.5){
              fine(v.y,0,1) = x3*2;
//                         printf("val32: v.x: %g, v.y: %g fine %g %s\n",v.x[],v.y[0,1],fine(v.y,0,1),v.x.name);
            }
          }
          
          else if (exist(v.x[1,0]) && exist(v.y[0,1])){  //Both left and right exist and mid is easy
            double m =  (1. - v.x[1,0])/(1. - v.y[0,1]);
            double val = 0.5*m + v.x[1,0];
//                     printf("val: %g, v.x: %g, v.y: %g m %g %s\n",val,v.x[1,0],v.y[0,1],m,v.x.name);
            if (val > 0.5 && val <= 1.){
              fine(v.x,1,1) = (val - 0.5)*2.;
//                         printf("val1: %g, v.x: %g, v.y: %g fine %g\n",val,v.x[1,0],v.y[0,1],fine(v.x,1,1));
            }
            else if (val >= 0. && val <= 0.5){
              fine(v.x,1,0) = val*2.;
//                         printf("val2: %g, v.x: %g, v.y: %g fine %g %s\n",val,v.x[1,0],v.y[0,1],fine(v.x,1,0),v.x.name);
            }
            double x3 = (0.5 - v.x[1,0])/(m);
            if (x3 >=0.5 && x3 < 1.){
              fine(v.y,1,1) = (x3 - 0.5)*2;
//                         printf("val33: v.x: %g, v.y: %g fine %g %s x3 %g val %g m %g\n", v.x[1,0], v.y[0,1], fine(v.y,1,1), v.y.name, x3, val, m);
            }
          }
          
          else if (exist(v.x[]) && exist(v.y[0,0])){  //Both left and right exist and mid is easy
            double m = (0. - v.x[])/(v.y[0,0] - 0.);
            double val = m*0.5 + v.x[];
//                     printf("val1: %g, v.x: %g, v.y: %g m %g %s\n",val,v.x[],v.y[0,1],m,v.x.name);        
            if (val > 0.5 && val <= 1.){
              fine(v.x,1,1) = (val - 0.5)*2.;
//                         printf("val1: %g, v.x: %g, v.y: %g fine %g\n",val,v.x[],v.y[0,1],fine(v.x,1,1));
            }
            else if (val >= 0. && val <= 0.5){
              fine(v.x,1,0) = val*2;
//                         printf("val2: %g, v.x: %g, v.y: %g fine %g\n",val,v.x[],v.y[0,1],fine(v.x,1,0));
            }
          }  
        }
  }

/**
  This event perform the 1-D advection of the Semushin marker points. */

event advect_x (i++) {
  foreach_dimension(){
    double tmp = 0;
    double y2, y1, x2, x1;
    y2 = y1 = x2 = x1 = 0.;
    
    foreach_face(){
      s_tmp.x[] = 0.;
      s_flag.x[] = 0.;
    }
    boundary ((scalar *){s_tmp, s_flag});
    
    foreach_face(y){ // This performs the horizontal advection
      if (exist(s.y[]))
        snew.y[] = s.y[] + myface_value(u.x,0)*dt/Delta;
    }
    boundary ((scalar *){snew});
    
    foreach_face(y){
      if (snew.y[] > 1.){  // The point is crossing a vertical line with $u>0$.
        tmp = snew.y[] - 1.;
        snew.y[] = 0.;
        s_tmp.y[] = tmp;
      }
      else if (snew.y[] < 0.){  // The point is crossing a vertical line with $u<0$.
        tmp = 1. + snew.y[];
        snew.y[] = 0.;
        s_tmp.y[] = -tmp;
      }
    }
    boundary ((scalar *){snew, s_tmp});
    
    foreach_face(y){
      if(exist(s_tmp.y[-1,0]))
        snew.y[] = s_tmp.y[-1,0];
      else if(existNeg(s_tmp.y[1,0])){
        snew.y[] = -s_tmp.y[1,0];
      }
    }
    boundary ((scalar *){snew});
    
    foreach_face(x){
      if(exist(s_tmp.y[-1,0])){    
        if (exist(s.x[-1,0])){
          y2 = s.x[-1,0]; y1 = 0.; x2 = myface_value(u.x,0)*dt/Delta; x1 = 1 + snew.y[];
          snew.x[] = intersect(x1, y1, x2, y2);
          //printf("i: %d y3 top is %g x1 is %g x2 is %g y1 is %g y2 is %g pid %d xy %g %g\n", i, snew.x[], x1, x2, y1, y2, pid(), x, y); 
        }
        
        if (exist(s.x[-1,-1])){
          y2 = 1 - s.x[-1,-1]; y1 = 0.; x2 = myface_value(u.x,0)*dt/Delta; x1 = 1 + snew.y[];
          s_tmp.x[] = 1. - intersect(x1, y1, x2, y2);
          //printf("i: %d y3 bot is %g x1 is %g x2 is %g y1 is %g y2 is %g pid %d xy %g %g\n", i, s_tmp.x[], x1, x2, y1, y2, pid(), x, y); 
        }
        
        if (exist(snew.y[-1,-1])){
          y2 = 1. ; y1 = 0.; x2 = snew.y[-1,-1]; x1 = 1 + snew.y[];
          s_tmp.x[] = 1. - intersect(x1, y1, x2, y2);
          //printf("i: %d y3 HBOT is %g x1 is %g x2 is %g y1 is %g y2 is %g pid %d xy %g %g \n", i, s_tmp.x[], x1, x2, y1, y2, pid(), x, y);       
        }
        
        if (exist(snew.y[-1,1])){
          y2 = 1. ; y1 = 0.; x2 = snew.y[-1,1]; x1 = 1 + snew.y[];
          snew.x[] =  intersect(x1, y1, x2, y2); 
          //printf("i: %d y3 HTOP is %g x1 is %g x2 is %g y1 is %g y2 is %g pid %d xy %g %g \n", i, snew.x[], x1, x2, y1, y2, pid(), x, y);   
        }
      }
      
      else if (existNeg(s_tmp.y[])){ //The point is crossing a vertical line with $u>0$.
        if (exist(s.x[1,0])){
          y2 = s.x[1,0]; y1 = 0.; x2 = 2. + myface_value(u.x,0)*dt/Delta; x1 =  snew.y[-1,0] ;
          snew.x[] = intersect(x1, y1, x2, y2);
          //printf("i: %d y3 top is %g x1 is %g x2 is %g y1 is %g y2 is %g\n", i, snew.x[0,0], x1, x2, y1, y2); 
        }
        
        if (exist(s.x[1,-1])){
          y2 = 1 - s.x[1,-1]; y1 = 0.; x2 = 2. + myface_value(u.x,0)*dt/Delta; x1 = snew.y[-1,0];
          s_tmp.x[] = 1. - intersect(x1, y1, x2, y2);
          //printf("i: %d y3 bot is %g x1 is %g x2 is %g y1 is %g y2 is %g\n", i, s_tmp.x[], x1, x2, y1, y2); 
        }
        
        if (exist(snew.y[0,-1])){
          y2 = 1. ; y1 = 0.; x2 = 1. + snew.y[0,-1]; x1 = snew.y[-1,0];
          s_tmp.x[] = 1. - intersect(x1, y1, x2, y2);
          //printf("i: %d y3 HBOT is %g x1 is %g x2 is %g y1 is %g y2 is %g\n", i, s_tmp.x[], x1, x2, y1, y2);       
        }
        
        if (exist(snew.y[0,1])){
          y2 = 1. ; y1 = 0.; x2 = 1. + snew.y[0,1]; x1 = snew.y[-1,0];
          snew.x[] = intersect(x1, y1, x2, y2);
          //printf("i: %d y3 HTOP is %g x1 is %g x2 is %g y1 is %g y2 is %g\n", i, snew.x[0,0], x1, x2, y1, y2);   
        }      
      }  
    }
    boundary ((scalar *){snew, s_tmp});
    
    foreach_face(x){
      if(exist(s_tmp.y[-1,0])){    
        //  Remove already existing points on the "crossed" line if they exist  
        if (exist(s.x[]))  
          s.x[] = 0.;
        if (exist(s.x[0,-1])) 
          s_flag.x[] = 0.5;  // Temporary(?) workaround
      }
      else if(existNeg(s_tmp.y[])) {    
        //  Remove already existing points on the "crossed" line if they exist  
        if (exist(s.x[]))  
          s.x[] = 0.;
        if (exist(s.x[0,-1])) 
          s_flag.x[] = 0.5;  // Temporary(?) workaround
      }
    }
    boundary ((scalar *){s, s_flag}); 
    
    foreach_face(x){
      if(exist(s_flag.x[0,1]))
        s.x[] = 0;
    }
    boundary ((scalar *){s});   
    
    foreach_face(x)
      if(exist(s_tmp.x[0,1]))
        snew.x[] = s_tmp.x[0,1];
      boundary ((scalar *){snew});
    
    /** We now have to compute the remaining intersections.*/
    foreach_face(x){
      if(myface_value(u.x,0)>=0. ){
        if (exist(s.x[0,0]) && exist(snew.y[-1,0])){ //TL
          y2 = s.x[0,0]; y1 = 0.; x2 = 1. + myface_value(u.x,0)*dt/Delta; x1 = snew.y[-1,0] ;
          snew.x[] = intersect(x1, y1, x2, y2);
          // printf("i: %d y3 TL is %g x1 is %g x2 is %g y1 is %g y2 is %g %d\n", i, snew.x[], x1, x2, y1, y2, pid());   
        }
        if (exist(s.x[0,0]) && exist(snew.y[-1,1])){ //BL
          y2 = s.x[0,0]; y1 = 1.; x2 = 1. + myface_value(u.x,0)*dt/Delta; x1 = snew.y[-1,1];
          snew.x[] =  intersect(x1, y1, x2, y2);
          // printf("i: %d y3 BL is %g x1 is %g x2 is %g y1 is %g y2 is %g %d \n", i, snew.x[], x1, x2, y1, y2, pid());   
        }
        if (exist(s.x[0,0]) && exist(s.x[-1,0])){ //VERT
          y2 = s.x[0,0]; y1 = s.x[-1,0]; x2 = 1. + myface_value(u.x,0)*dt/Delta; x1 = myface_value(u.x,-1)*dt/Delta;
          snew.x[] = intersect(x1, y1, x2, y2);
          // printf("i: %d y3 VERT is %g x1 is %g x2 is %g y1 is %g y2 is %g %d\n", i, snew.x[], x1, x2, y1, y2, pid());    
        }
        if (exist(s.x[0,0]) && exist(snew.x[-1,0])){ //VERT
          y2 = s.x[0,0]; y1 = snew.x[-1,0]; x2 = 1. + myface_value(u.x,0)*dt/Delta; x1 = 0.;
          snew.x[] = intersect(x1, y1, x2, y2);
          // printf("i: %d y3 VERT is %g x1 is %g x2 is %g y1 is %g y2 is %g %d\n", i, snew.x[], x1, x2, y1, y2, pid());    
        }
      }
      else{
        if (exist(s.x[0,0]) && exist(s.x[1,0])){ //VERT                               
          y2 = s.x[0,0]; y1 = s.x[1,0]; x2 =  1. + myface_value(u.x,0)*dt/Delta; x1 = 2. + myface_value(u.x,1)*dt/Delta;                             
          snew.x[] = intersect(x1, y1, x2, y2);                                                                      
          // printf("i: %d y3 VERT is %g x1 is %g x2 is %g y1 is %g y2 is %g\n", i, snew.x[], x1, x2, y1, y2);    
        }
        if (exist(s.x[0,0]) && exist(snew.x[1,0])){ //VERT                               
          y2 = s.x[0,0]; y1 = snew.x[1,0]; x2 =  1. + myface_value(u.x,0)*dt/Delta; x1 = 2.;                             
          snew.x[] = intersect(x1, y1, x2, y2);                                                                      
          // printf("i: %d y3 VERT is %g x1 is %g x2 is %g y1 is %g y2 is %g\n", i, snew.x[], x1, x2, y1, y2);    
        }          
        if (exist(s.x[0,0]) && exist(snew.y[0,0])){ //TR                              
          y2 = s.x[0,0]; y1 = 0.; x2 = 1. + myface_value(u.x,0)*dt/Delta; x1 = 1. + snew.y[0,0];                           
          snew.x[] = intersect(x1, y1, x2, y2);                                                                  
          // printf("i: %d y3 TR is %g x1 is %g x2 is %g y1 is %g y2 is %g\n", i, snew.x[], x1, x2, y1, y2);                  
        }                                                                                                               
        if (exist(s.x[0,0]) && exist(snew.y[0,1])){ //BR                              
          y2 = s.x[0,0] ; y1 = 1.; x2 = 1.+ myface_value(u.x,0)*dt/Delta; x1 = 1. + snew.y[0,1];                   
          snew.x[] = intersect(x1, y1, x2, y2);                                                                  
          // printf("i: %d y3 BR is %g x1 is %g x2 is %g y1 is %g y2 is %g \n", i, snew.x[], x1, x2, y1, y2);    
        }
      }
    }
    boundary ((scalar *){snew});
    
    foreach_face(x)
      s.x[] = snew.x[];
    foreach_face(y)
      s.y[] = snew.y[];
    boundary ((scalar *){s});
    
    foreach_face()
      snew.x[] = 0.;
    boundary ((scalar *){snew});
  }
}

/** The interface can be shown with gnuplot */
#ifdef OUT_GNUPLOT

event output_interface(i++;i <= IT) {  

  char out[100];
  const char testName[20] = TEST;
  snprintf(out, sizeof(out), "%s_%d_%d.dat", testName, Nl, i);

  FILE * fp1;
  fclose(fopen(out, "w"));
  fp1 =fopen(out, "a");
  
  foreach(serial){
    double xx[2], yy[2];
    xx[0] = xx[1] = yy[0] = yy[1] = 0.;
    int ii = 0;
    if (exist(s.x[0,0]) || exist(s.x[1,0]) || exist(s.y[0,0]) || exist(s.y[0,1])){
      if (exist(s.x[0,0])){
        yy[ii] = y - Delta/2. + s.x[0,0]*Delta;
        xx[ii] = x - Delta/2.;
        ii++;
      }
      if (exist(s.x[1,0])){
        yy[ii] = y - Delta/2. + s.x[1,0]*Delta;
        xx[ii] = x + Delta/2.;
        ii++;
      }
      if (exist(s.y[0,0])){
        xx[ii] = x - Delta/2. + s.y[0,0]*Delta;
        yy[ii] = y - Delta/2.;
        ii++;
      }
      if (exist(s.y[0,1])){
        xx[ii] = x - Delta/2. + s.y[0,1]*Delta;
        yy[ii] = y + Delta/2.;
        ii++;
      }      
      fprintf (fp1, "%g %g\n%g %g\n\n", xx[0], yy[0], xx[1], yy[1]);
      fflush(fp1);
    }
  }
  fclose(fp1);  
}
#endif


/**
The representation of the interface with Semushin markers can be used to create
the associated VOF field. We use the tag function to identify the 3 regions:
phase 1, phase 0 and the cells cut by the interface.
The VOF field is here used to easily compute the area of the reference phase. */

#ifdef SEMU2VOF 
#include "tag.h"
  scalar ff[], f1[];
event semu2vof(i++) {

  foreach(){
    if (!exist(s.x[0,0]) && !exist(s.x[1,0]) && !exist(s.y[0,0]) && !exist(s.y[0,1])){  //Full or empty cells
      f1[] = 2.;
    }
    else 
      f1[] = 0.;
  }

  foreach()
    ff[] = f1[] > 1.;
  int n = tag(ff);
  foreach()
    ff[] = ff[] - 1.;
  assert(n == 2); 
  
  foreach_dimension(){
    foreach(){
      if (exist(s.x[0,0]) && exist(s.x[1,0])){  // Opposite sides (L-R + T-B)    
        if (ff[0,-1] > 0.5)
          ff[] = (s.x[0,0] + s.x[1,0])/2.;
        else
          ff[] = 1. - (s.x[0,0] + s.x[1,0])/2.; 
      }
      
      if (exist(s.x[0,0]) && exist(s.y[0,0])){  // Left and bottom   
        if (ff[-1,-1] > 0.5)
          ff[] = (s.x[0,0]*s.y[0,0])/2.;
        else if (ff[-1,-1] > -1. + 1.e-5)
          ff[] = 1. - (s.x[0,0]*s.y[0,0])/2.;
        else if (ff[1,1] > 0.5)
          ff[] = 1. - ((1. - s.x[0,0])*s.y[0,0])/2.;
        else 
          ff[] = ((1. - s.x[0,0])*s.y[0,0])/2.;        
      }
      

      if (exist(s.x[0,0]) && exist(s.y[0,1])){  // Left and top + right and bot   
        if (ff[-1,1] > 0.5)               
          ff[] = ((1. - s.x[0,0])*s.y[0,1])/2.;
        else if (ff[-1,1] > -1. + 1.e-5)
          ff[] = 1. - ((1. - s.x[0,0])*s.y[0,1])/2.;
        else if (ff[1,-1] > 0.5)
          ff[] = 1. - ((1. - s.x[0,0])*s.y[0,1])/2.;
        else 
          ff[] = ((1. - s.x[0,0])*s.y[0,1])/2.;
      } 
      
      if (exist(s.x[1,0]) && exist(s.y[0,1])){  // Right and top   
        if (ff[1,1] > 0.5)
          ff[] = ((1. - s.x[1,0])*(1. - s.y[0,1]))/2.;
        else if (ff[1,1] > -1. + 1.e-5)
          ff[] = 1. - ((1. - s.x[1,0])*(1. - s.y[0,1]))/2.;
        else if (ff[-1,-1] > 0.5)
          ff[] = 1. - ((1. - s.x[1,0])*s.y[0,1])/2.;
        else 
          ff[] = ((1. - s.x[1,0])*s.y[0,1])/2.;        
      }        
    }
    foreach()
      f[] = ff[];
  }
  double area = 0.;
  foreach(reduction(+:area)){
    area += f[]*sq(Delta);
  }
  printf("Area is %f \n", area);
  
}

#endif