/**
# Header File - Direct Laplacian Computation

We compute the Direct Laplacian of a function using a 4th order scheme & a 2nd order scheme.
The computation is done both for a Multigrid as well as a Quadtree */

static inline double bicubic(Point point, scalar A){

 #if dimension == 1
  return((-7*coarse(A,-child.x)+105*coarse(A,0)+35*coarse(A,child.x)-5*coarse(A,2*child.x))/128.);

 #elif dimension == 2
  return(( -7.*(-7*coarse(A,-child.x,-child.y,0)+105*coarse(A,0,-child.y,0)+35*coarse(A,child.x,-child.y,0)-5*coarse(A,2*child.x,-child.y,0)) + 105.*(-7*coarse(A,-child.x,0,0)+105*coarse(A,0,0,0)+35*coarse(A,child.x,0,0)-5*coarse(A,2*child.x,0,0)) + 35.*(-7*coarse(A,-child.x,child.y,0)+105*coarse(A,0,child.y,0)+35*coarse(A,child.x,child.y,0)-5*coarse(A,2*child.x,child.y,0)) - 5.*(-7*coarse(A,-child.x,2*child.y,0)+105*coarse(A,0,2*child.y,0)+35*coarse(A,child.x,2*child.y,0)-5*coarse(A,2*child.x,2*child.y,0)) )/(128.*128.));
 #endif

}

static inline void refine_bicubic (Point point, scalar s)
{
  foreach_child()
    s[] = bicubic(point, s);
}


static inline double biquintic(Point point, scalar A){

#if dimension == 1
  printf("\n %g %g",x,y);
  return((35.*coarse(A,-2*child.x) - 252.*coarse(A,-child.x) + 1890.*coarse(A,0) + 420.*coarse(A,child.x) - 45.*coarse(A,2*child.x))/2048.);
   
#elif dimension == 2
  return 0;
#endif
}

static inline void refine_biquintic(Point point, scalar s){
  foreach_child()
    s[] = biquintic(point,s);
}


static inline double biquartic(Point point, scalar A){

 #if dimension == 1
   return ((105.*coarse(A,-2*child.x,0,0)-756*coarse(A,-child.x,0,0)+5670*coarse(A,0,0,0)+1260*coarse(A,child.x,0,0)-135*coarse(A,2*child.x,0,0))/(24*256));

 #elif dimension == 2
  return((105.*(105.*coarse(A,-2*child.x,-2*child.y,0)-756*coarse(A,-child.x,-2*child.y,0)+5670*coarse(A,0,-2*child.y,0)+1260*coarse(A,child.x,-2*child.y,0)-135*coarse(A,2*child.x,-2*child.y,0)) -756.*(105.*coarse(A,-2*child.x,-child.y,0)-756*coarse(A,-child.x,-child.y,0)+5670*coarse(A,0,-child.y,0)+1260*coarse(A,child.x,-child.y,0)-135*coarse(A,2*child.x,-child.y,0)) + 5670.*(105.*coarse(A,-2*child.x,0,0)-756*coarse(A,-child.x,0,0)+5670*coarse(A,0,0,0)+1260*coarse(A,child.x,0,0)-135*coarse(A,2*child.x,0,0)) + 1260.*(105.*coarse(A,-2*child.x,child.y,0)-756*coarse(A,-child.x,child.y,0)+5670*coarse(A,0,child.y,0)+1260*coarse(A,child.x,child.y,0)-135*coarse(A,2*child.x,child.y,0)) - 135.*(105.*coarse(A,-2*child.x,2*child.y,0)-756*coarse(A,-child.x,2*child.y,0)+5670*coarse(A,0,2*child.y,0)+1260*coarse(A,child.x,2*child.y,0)-135*coarse(A,2*child.x,2*child.y,0)))/(256.*256.*24.*24.));
 #endif

}

static inline void refine_biquartic (Point point, scalar s)
{
  foreach_child()
    s[] = biquartic(point, s);
}

void Laplacian2nd (scalar *al, scalar *lapl){

   scalar a = al[0], Lap = lapl[0];
  
#if TREE
  
  face vector g[];
  foreach_face()
    g.x[] = (a[] - a[-1])/Delta;
  boundary_flux ({g});
  foreach () {
    Lap[] = 0;
    foreach_dimension()
      Lap[] += (g.x[1] - g.x[])/Delta;
     }
  
#else

  foreach () {
    Lap[] = 0;
    foreach_dimension()
      Lap[] += (a[1] - 2*a[] + a[-1])/sq(Delta);
  }

#endif

boundary({Lap});
}

void Laplacian4th (scalar * al, scalar * lapl){

   scalar a = al[0], Lap = lapl[0];

#if TREE
 
   face vector g[];
   foreach_face(){
      g.x[] = 0;
      g.x[] += (a[-2,2]-27*a[-1,2]+27*a[0,2]-a[1,2] )*(-17.)/(5760.*24.*Delta);
      g.x[] += (a[-2,1]-27*a[-1,1]+27*a[0,1]-a[1,1] )*(308.)/(5760.*24.*Delta);      
      g.x[] += (a[-2,0]-27*a[-1,0]+27*a[0,0]-a[1,0] )*(5178.)/(5760.*24.*Delta);
      g.x[] += (a[-2,-1]-27*a[-1,-1]+27*a[0,-1]-a[1,-1] )*(308.)/(5760.*24.*Delta);
      g.x[] += (a[-2,-2]-27*a[-1,-2]+27*a[0,-2]-a[1,-2] )*(-17.)/(5760.*24.*Delta);
    }
   boundary_flux({g});
   foreach(){
       Lap[] = 0;
       foreach_dimension()
           Lap[] += (g.x[1]-g.x[])/Delta;
     }        

#else

   foreach(){
      Lap[] = 0;
      foreach_dimension(){
         Lap[] += ( (a[-1,2]-27*a[0,2]+27*a[1,2]-a[2,2]) - (a[-2,2]-27*a[-1,2]+27*a[0,2]-a[1,2]) )*(-17.)/(5760.*24.*sq(Delta));
         Lap[] += ( (a[-1,1]-27*a[0,1]+27*a[1,1]-a[2,1]) - (a[-2,1]-27*a[-1,1]+27*a[0,1]-a[1,1]) )*(308.)/(5760.*24.*sq(Delta));
         Lap[] += ( (a[-1,0]-27*a[0,0]+27*a[1,0]-a[2,0]) - (a[-2,0]-27*a[-1,0]+27*a[0,0]-a[1,0]) )*(5178.)/(5760.*24.*sq(Delta));
         Lap[] += ( (a[-1,-1]-27*a[0,-1]+27*a[1,-1]-a[2,-1]) - (a[-2,-1]-27*a[-1,-1]+27*a[0,-1]-a[1,-1]) )*(308.)/(5760.*24.*sq(Delta));
         Lap[] += ( (a[-1,-2]-27*a[0,-2]+27*a[1,-2]-a[2,-2]) - (a[-2,-2]-27*a[-1,-2]+27*a[0,-2]-a[1,-2]) )*(-17.)/(5760.*24.*sq(Delta));
       }
   }
#endif

boundary({Lap});
}