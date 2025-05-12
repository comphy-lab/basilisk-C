/**
#Curvature of a level_set function
*/

#include "alex_functions.h"
#include "fractions.h"
#include "curvature.h" // we use prolongation and restriction operator defined
// in curvature.H

struct Curvature_LS{
  scalar LS, curve_LS;
};

void curvature_LS(struct Curvature_LS p){

  scalar LS = p.LS, curve_LS = p.curve_LS;
  vector gr_LS[];
#if TREE
  curve_LS.refine = curvature_prolongation;
  curve_LS.restriction = curvature_restriction;
#endif

  scalar c1[];
  foreach(){
    double norm = 1.e-15;
    foreach_dimension(){
      double temp = (LS[1,0]-LS[-1,0])/(2.*Delta);
      gr_LS.x[] = temp;
      norm += sq(temp); 
    }
    norm = sqrt(norm);
    foreach_dimension(){
      gr_LS.x[] /= norm;
    }
  }
  boundary((scalar *){gr_LS});
  
#if dimension == 2
  scalar phixx[], phixy[], phiyy[];

  foreach(){
    phixx[] = (gr_LS.x[1,0] - gr_LS.x[-1,0])/(2.*Delta); 
    phiyy[] = (gr_LS.y[0,1] - gr_LS.y[0,-1])/(2.*Delta); 
    phixy[] = ((gr_LS.y[1,0] - gr_LS.y[-1,0])/(2.*Delta) +  
              (gr_LS.x[0,1] - gr_LS.x[0,-1])/(2.*Delta))/2.;
  }

  boundary({phixx, phixy, phiyy});
  restriction({phixx, phixy, phiyy});


/**
In 2D, we use the following formula for curvature calculation :
$$
\Kappa = \dfrac{\phi_y^2\phi_{xx} - 2\phi_x\phi_y\phi_{xy}+\phi_x^2\phi_{yy}}
{\left| \phi_x^2 + \phi_y^2\right|^{3/2}}
$$
*/
  foreach(){
    c1[]= (sq(gr_LS.y[])*phixx[] 
      - 2.*gr_LS.x[]*gr_LS.y[]*phixy[] 
      + sq(gr_LS.x[])*phiyy[])/
    powf(1.e-15 + (sq(gr_LS.x[]) + sq(gr_LS.y[])),1.5);
  }
  boundary({c1});
  restriction({c1});

#else
//FIXME : STILL BUGGED I GUESS
  foreach(){
    double sum = 0;
    foreach_dimension(){
      sum+=(gr_LS.x[1] - gr_LS.x[-1])/(2*Delta);
    }
    c1[] = sum; 
  }

#endif

  foreach (noauto) { // fixme: this is a hack
    if(interfacial(point,cs)){
      coord n, p;
      embed_geometry(point, &p, &n);
/**
The interpolation stencil is chosen according to the position of the face
centroid.
*/

      coord p_interp = {p.x, p.y, p.z};
#ifdef BICUBIC
      int Stencil[2] = {-1,-1};
      curve_LS[] = -bicubic( point , c1, Stencil, p_interp);
#elif QUADRATIC
#if dimension == 2
      curve_LS[] = -mybiquadratic( point , c1, p_interp,0);
#else
      curve_LS[] = -mytriquadratic( point , c1, p_interp);
#endif
#endif
    }
    else{
      curve_LS[] = nodata;
    }
  }
  curve_LS.dirty = true; // fixme: this goes with the hack above
  boundary({curve_LS});
  // restriction({curve_LS});
}
