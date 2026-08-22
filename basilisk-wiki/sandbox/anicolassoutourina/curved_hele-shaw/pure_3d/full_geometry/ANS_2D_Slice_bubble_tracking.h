/* ---------------------------------------------------------- */
/* -- ANS : Header with all custom functions necessary ------ */
/* -- to perform data extraction of a 2D slice for bubble --- */
/* -- center tracking. -------------------------------------- */

/* Same as facets function in regular Basilisk : https://basilisk.fr/src/geometry.h#facets */
/* Except we force it to operate as it does in 2D for a 3D simulation */
int Forced2D_facets (coord n, double alpha, coord p[2]){
  int i = 0;
  for (double s = -0.5; s<= 0.5; s+=1.){
    // foreach_dimension done manually to force it to be 2D in XY plane
    if (fabs (n.y) > 1e-4 && i < 2){
      double a = (alpha - s*n.x)/n.y ;
      if (a >= -0.5 && a <= 0.5){
        p[i].x = s ;
        p[i++].y = a ;
      }
    }
    if (fabs (n.x) > 1e-4 && i < 2){
      double a = (alpha - s*n.y)/n.x ;
      if (a >= -0.5 && a <= 0.5){
        p[i].y = s ;
        p[i++].x = a ;
      }
    }
  }
  return i ;
}

/* We force interface_normal : https://basilisk.fr/src/fractions.h#interface_normal to behave in 2D mode */

#define NOT_ZERO 1e-30

/*-----------------------------------------------------*
 *MYC - Mixed Youngs and Central Scheme (2D)           *
 *-----------------------------------------------------*/
coord Forced2D_mycs (Point point, scalar c)
{
  int ix;
  double c_t,c_b,c_r,c_l;
  double mx0,my0,mx1,my1,mm1,mm2;
  
  /* top, bottom, right and left sums of c values */
  c_t = c[-1,1] + c[0,1] + c[1,1];
  c_b = c[-1,-1] + c[0,-1] + c[1,-1];
  c_r = c[1,-1] + c[1,0] + c[1,1];
  c_l = c[-1,-1] + c[-1,0] + c[-1,1];

  /* consider two lines: sgn(my) Y =  mx0 X + alpha,
     and: sgn(mx) X =  my0 Y + alpha */ 
  mx0 = 0.5*(c_l - c_r);
  my0 = 0.5*(c_b - c_t);

  /* minimum coefficient between mx0 and my0 wins */
  if (fabs(mx0) <= fabs(my0)) {
    my0 = my0 > 0. ? 1. : -1.;
    ix = 1;
  }
  else {
    mx0 = mx0 > 0. ? 1. : -1.;
    ix = 0;
  }

  /* Youngs' normal to the interface */
  mm1 = c[-1,-1] + 2.0*c[-1,0] + c[-1,1];
  mm2 = c[1,-1] + 2.0*c[1,0] + c[1,1];
  mx1 = mm1 - mm2 + NOT_ZERO;
  mm1 = c[-1,-1] + 2.0*c[0,-1] + c[1,-1];
  mm2 = c[-1,1] + 2.0*c[0,1] + c[1,1];
  my1 = mm1 - mm2 + NOT_ZERO;

  /* choose between the best central and Youngs' scheme */ 
  if (ix) {
    mm1 = fabs(my1);
    mm1 = fabs(mx1)/mm1;
    if (mm1 > fabs(mx0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }
  else {
    mm1 = fabs(mx1);
    mm1 = fabs(my1)/mm1;
    if (mm1 > fabs(my0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }
	
  /* normalize the set (mx0,my0): |mx0|+|my0|=1 and
     write the two components of the normal vector  */
  mm1 = fabs(mx0) + fabs(my0);
  coord n = {mx0/mm1, my0/mm1, 0};
  
  return n;
}

auto macro coord Forced2D_interface_normal (Point p, scalar c) {
  return Forced2D_mycs (p, c);
}

/* We force plane_alpha : https://basilisk.fr/src/geometry.h#plane_alpha to behave in 2D mode */

double Forced2D_plane_alpha (double c, coord n)
{
  double alpha, n1, n2;
  
  n1 = fabs (n.x); n2 = fabs (n.y);
  if (n1 > n2)
    swap (double, n1, n2);

  c = clamp (c, 0., 1.);
  double v1 = n1/2.;
  if (c <= v1/n2)
    alpha = sqrt (2.*c*n1*n2);
  else if (c <= 1. - v1/n2)
    alpha = c*n2 + v1;
  else
    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;

  return alpha - (n.x + n.y)/2.;
}

#include "curvature.h" // necessary to get the bool interfacial function

double ANS_interface_x_center (scalar f, int maxlevel, double domainsize, double x_offset, double y_offset, double z_offset){
  // We try to use foreach_region to force to work within a slice or pseudo-slice at Z = 0 (with offsets), assuming a centered domain
  double sx = 0., sl = 0. ;
  coord box[2] = {{x_offset-domainsize/2,y_offset-domainsize/2,z_offset},{x_offset+domainsize/2,y_offset+domainsize/2,z_offset}} ; // slice at exactly z = z_offset
  coord n_box = {pow(2,maxlevel), pow(2,maxlevel)}; // nb of points
  coord p_box ;
  foreach_region(p_box, box, n_box, reduction(+:sx) reduction(+:sl)){
    if (interfacial(point, f)){
      coord n = Forced2D_interface_normal(point, f);
      double alpha = Forced2D_plane_alpha(f[],n);
        
      coord p[2];
      int m = Forced2D_facets(n, alpha, p);
      
      //if (m==dimension_adjustor*m_compare){ // FIXME : "the dimensional constraints below are not compatible" for L0, idk how to fix this. Problem only there for non-cubic 3D simulations, this line works for cubic 3D
      if(true){
        // foreach_dimension done manually
        p[0].x = x + p[0].x*Delta ;
        p[1].x = x + p[1].x*Delta ;
        
        p[0].y = y + p[0].y*Delta ;
        p[1].y = y + p[1].y*Delta ;
          
        double l = sqrt(sq(p[1].x - p[0].x) + sq(p[1].y - p[0].y));
          
        double xm = 0.5*(p[0].x + p[1].x) ;
          
        sx += xm*l ;
        sl += l;
      }
    }
  }
  return sl > 0. ? sx/sl : nodata ;
}

double ANS_interface_y_center (scalar f, int maxlevel, double domainsize, double x_offset, double y_offset, double z_offset){
  double sy = 0., sl = 0. ;
  coord box[2] = {{x_offset-domainsize/2,y_offset-domainsize/2,z_offset},{x_offset+domainsize/2,y_offset+domainsize/2,z_offset}} ;
  //coord box[2] = {{-domainsize/2,-domainsize/2,z_offset},{domainsize/2,domainsize/2,z_offset}} ;
  coord n_box = {pow(2,maxlevel), pow(2,maxlevel)}; // nb of points
  coord p_box ;
  foreach_region(p_box, box, n_box, reduction(+:sy) reduction(+:sl)){
    if (interfacial(point, f)){
      coord n = Forced2D_interface_normal(point, f);
      double alpha = Forced2D_plane_alpha(f[],n);
        
      coord p[2];
      int m = Forced2D_facets(n, alpha, p);
        
      //if (m==2){ // FIXME : "the dimensional constraints below are not compatible" for L0, idk how to fix this. Problem only there for non-cubic 3D simulations, this line works for cubic 3D
      if(true){
        // foreach_dimension done manually
        p[0].x = x + p[0].x*Delta ;
        p[1].x = x + p[1].x*Delta ;
        
        p[0].y = y + p[0].y*Delta ;
        p[1].y = y + p[1].y*Delta ;
          
        double l = sqrt(sq(p[1].x - p[0].x) + sq(p[1].y - p[0].y));
          
        double ym = 0.5*(p[0].y + p[1].y) ;
          
        sy += ym*l ;
        sl += l;
      }
    }
  }
  return sl > 0. ? sy/sl : nodata ;
}