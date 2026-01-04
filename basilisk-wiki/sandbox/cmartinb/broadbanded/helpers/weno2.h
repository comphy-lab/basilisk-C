/**
#WENO2

Taken from Gibou reinitialization papers.

##Minmod limiter
*/


double my_minmod(double a,double b){
  if(a*b>0){
    if(fabs(a) < fabs(b))
      return a;
    if(fabs(a) > fabs(b))
      return b;
  }
  return 0;
}

foreach_dimension()
static inline double WENOdiff_x(Point point, scalar s, int i){
  double s1 = (s[2*i,0,0] + s[] - 2*s[i,0,0])/Delta; 
  double s2 = (s[1,0,0] + s[-1,0,0] - 2*s[])/Delta;
  return i*((s[i,0,0] - s[])/Delta -my_minmod(s1,s2)/2.);
}