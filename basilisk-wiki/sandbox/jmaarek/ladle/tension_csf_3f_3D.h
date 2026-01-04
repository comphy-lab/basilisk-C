#include "curvature.h"

  attribute {
  double sigma;
}

event stability (i++)
{
  double T1 = sqrt((rho1 + rho2)/2*pow(L0/(1 << grid->maxdepth), 3)/M_PI/(f1.sigma+f2.sigma));
  double T2 = sqrt((rho1 + rho3)/2*pow(L0/(1 << grid->maxdepth), 3)/M_PI/(f1.sigma+f3.sigma));
  double T3 = sqrt((rho2 + rho3)/2*pow(L0/(1 << grid->maxdepth), 3)/M_PI/(f2.sigma+f3.sigma));
  double dt = min(min(T1,T2), T3);
  if (dt < dtmax)
    dtmax = dt;
}

event defaults (i = 0) {  
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
  }
}

scalar triple_point_tag[];


foreach_dimension()
double isotropic_gradient_x(Point point, scalar input){

  double alpha = 0.245;
  double beta = 0.085;

return  (1*(input[1,0,0]-input[-1,0,0] )
	    + alpha*(input[1,0,1]-input[-1,0,1] + input[1,0,-1]-input[-1,0,-1] + input[1,1,0]-input[-1,1,0] + input[1,-1,0]-input[-1,-1,0])
	    + beta*(input[1,1,1]-input[-1,1,1] + input[1,1,-1]-input[-1,1,-1] + input[1,-1,1]-input[-1,-1,1] + input[1,-1,-1]-input[-1,-1,-1]))/(2+8*alpha+8*beta);
}

/*
foreach_dimension()
double isotropic_gradient_x(Point point, scalar input){
return  (input[1,0,0]-input[-1,0,0] )/2.0;
}*/


event acceleration (i++){
  /**
  We check for all VOF interfaces for which $\phi$ is allocated. The
  corresponding volume fraction fields will be stored in *list*. */

  scalar * list = NULL;


  /** smooth f */
  scalar smoothf[];
  // normal vector
  vector normfff[];

  scalar curv_c[];
  scalar curv_c_hf[];
  
  
  for (scalar f in interfaces){
      list = list_add (list, f);

      foreach()
        f[] = clamp (f[], 0., 1.);

  face vector ia = a;

  /** New method to compute curvature from the divergence of the interface normal. Normal is computed from smoothed volume fraction **/



  /** smooth the f based on isotropic direction*/

    foreach(){
     smoothf[] =   (64.*f[]+
               16.*(f[-1,0,0]+f[1,0,0]+f[0,-1,0]+f[0,1,0]+f[0,0,-1]+f[0,0,1])+
                4.*(f[-1,0,-1]+f[1,0,-1]+f[0,-1,-1]+f[0,1,-1]+f[-1,-1,0]+f[-1,1,0]+f[1,-1,0]+f[1,1,0]+f[-1,0,1]+f[1,0,1]+f[0,-1,1]+f[0,1,1])+
                1.*(f[-1,-1,-1]+f[-1,1,-1]+f[1,-1,-1]+f[1,1,-1]+f[-1,-1,1]+f[-1,1,1]+f[1,-1,1]+f[1,1,1]))/216.;
     }


  boundary({smoothf});

  /** evaluate the gradient of f $\nabla f$ at cell center by average*/
  foreach(){

	normfff.x[] = isotropic_gradient_x(point, smoothf);
        normfff.y[] = isotropic_gradient_y(point, smoothf);
        normfff.z[] = isotropic_gradient_z(point, smoothf);

      double molef=sqrt(normfff.x[]*normfff.x[]+normfff.y[]*normfff.y[] + normfff.z[]*normfff.z[]);

      normfff.x[]= molef > 0.0 ? -normfff.x[]/molef : 0.0;
      normfff.y[]= molef > 0.0 ? -normfff.y[]/molef : 0.0;
      normfff.z[]= molef > 0.0 ? -normfff.z[]/molef: 0.0;
      }
  boundary({normfff});
  // compute the curvature 

  foreach(){
    curv_c[] = (isotropic_gradient_x(point, normfff.x) + isotropic_gradient_y(point, normfff.y) + isotropic_gradient_z(point, normfff.z))/Delta ;
  }

  curvature(f, curv_c_hf, sigma = 1.0);

  #if TREE
    f.prolongation = p.prolongation;
    triple_point_tag.prolongation = refine_injection;
    f.dirty = true; // boundary conditions need to be updated
  #endif

  //tag a region near the triple point
    //tag a region near the triple point
  foreach(){
    double tag_temp[3];
    triple_point_tag[] = 0.0;

    for(int ii = -2; ii <= 0; ii++){
      for(int jj = -2; jj <= 0; jj++){
        for(int mm = -2; mm <= 0; mm++){
          tag_temp[0] = 0.0;
          tag_temp[1] = 0.0;
          tag_temp[2] = 0.0;

          for(int kk = 0; kk <= 2; kk++){
            for(int ll = 0; ll <= 2; ll++){
              for(int nn = 0; nn <= 2; nn++){
                if ((f1[ii+kk,jj+ll, mm+nn] > 0.0) && (f1[ii+kk,jj+ll, mm+nn] < 1.0) && (f2[ii+kk,jj+ll, mm+nn] > 0.0) && (f2[ii+kk,jj+ll, mm+nn] < 1.0))
                  tag_temp[0] = 1.0;
                if ((f1[ii+kk,jj+ll, mm+nn] > 0.0) && (f1[ii+kk,jj+ll, mm+nn] < 1.0) && (f3[ii+kk,jj+ll, mm+nn] > 0.0) && (f3[ii+kk,jj+ll, mm+nn] < 1.0))
                  tag_temp[1] = 1.0;
                if ((f2[ii+kk,jj+ll, mm+nn] > 0.0) && (f2[ii+kk,jj+ll, mm+nn] < 1.0) && (f3[ii+kk,jj+ll, mm+nn] > 0.0) && (f3[ii+kk,jj+ll, mm+nn] < 1.0))
                  tag_temp[2] = 1.0;
              }
            }
          }

          if((tag_temp[0]+tag_temp[1]+tag_temp[2]) > 2 && Delta == (L0/(1 << grid->maxdepth)))
            triple_point_tag[] = 1.0;
        }
      }
    }
  }

  foreach_face(){
      double curv_f = 0.0;
      //if we are near the triple point we compute curvature from the divergence of the normal (brackbrill scheme)
      if((triple_point_tag[-1] + triple_point_tag[0]) > 0.0){
        curv_f = (curv_c[-1]*(1.-2.*fabs((0.5-smoothf[-1]))) + curv_c[0]*(1.-2.*fabs((0.5-smoothf[0]))))/((1.-2.*fabs((0.5-smoothf[-1]))) + (1.-2.*fabs((0.5-smoothf[0])))+1e-32);
      }
      else{
        //if we are far from the triple point we compute curvature with height functions (basilisk scheme)
        curv_f =
          (curv_c_hf[] < nodata && curv_c_hf[-1] < nodata) ?
          (curv_c_hf[] + curv_c_hf[-1])/2. :
          curv_c_hf[] < nodata ? curv_c_hf[] :
          curv_c_hf[-1] < nodata ? curv_c_hf[-1] :
          0.;
      }

      ia.x[] += alpha.x[]/(fm.x[])*f.sigma*curv_f*(f[]-f[-1])/Delta;
  }
  boundary((scalar *){ia});


  #if TREE
  for (scalar f in list) {
    f.prolongation = fraction_refine;
    f.dirty = true; // boundary conditions need to be updated
  }
#endif
}

  free (list);
}