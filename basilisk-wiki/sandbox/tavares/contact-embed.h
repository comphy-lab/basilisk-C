/**
# Contact angles on an embedded boundary

This header file implemented contact angles for [VOF
interfaces](/src/vof.h) on [embedded boundaries](/src/embed.h).

The contact angle is defined by this field and can be constant or
variable. */

(const) scalar contact_angle;
scalar tag[];


/**
This function returns the properly-oriented normal of an interface
touching the embedded boundary. `ns` is the normal to the embedded
boundary, `nf` the VOF normal, not taking into account the contact
angle, and `angle` the contact angle. */

static inline coord normal_contact (coord ns, coord nf, double angle)
{
  coord n;
  if (dimension == 2){
    if (- ns.x*nf.y + ns.y*nf.x > 0) { // 2D cross product
      n.x = - ns.x*cos(angle) + ns.y*sin(angle);
      n.y = - ns.x*sin(angle) - ns.y*cos(angle);
    }
    else {
      n.x = - ns.x*cos(angle) - ns.y*sin(angle);
      n.y =   ns.x*sin(angle) - ns.y*cos(angle);
    }
  }

  else { //3D
  /** Two versions of the 3D normal *n* calculation has been tested */
  
/** Here is a first version using spherical coordinates */
    foreach_dimension()  //Axis angle
      if (ns.x != 0) ns.x = -ns.x ; //We take the normal pointing toward the fluid
    double gamma = atan2(ns.z,ns.x); 
    double beta = asin(ns.y);  
    coord tau2 ={-sin(gamma), 0, cos(gamma)};
    coord tau1; 
    foreach_dimension() //cross product
      tau1.x = ns.y*tau2.z - ns.z*tau2.y; 
    double phi1 =0, phi2 =0;
  //double phi1, phi2;
    foreach_dimension(){
      phi1 += nf.x*tau1.x;
      phi2 += nf.x*tau2.x;
    }
    double phi = atan2(phi2,phi1); 

//general configurations using Axis angle with spherical coordinates
    n.x = sin(angle)*cos(phi)*sin(beta)*cos(gamma) + cos(angle)*cos(beta)*cos(gamma) - sin(angle)*sin(phi)*sin(gamma);
    n.y = cos(angle)*sin(beta) - sin(angle)*cos(phi)*cos(beta);
    n.z = sin(angle)*cos(phi)*sin(beta)*sin(gamma) + cos(angle)*cos(beta)*sin(gamma) + sin(angle)*sin(phi)*cos(gamma);
  
/** Another version using Quaternions */
     coord rot_axis; //Quaternion and rotation matrix 
     foreach_dimension() //rotation axis to define a quaternion
       rot_axis.x = ns.y*nf.z - ns.z*nf.y; 
     normalize(&rot_axis);
     coord p1 = {ns.x, ns.y, ns.z};  //original position of point

     double re = cos(angle/2.); //real part of quaternion
     double x = rot_axis.x*sin(angle/2.); //imaginary i part of quaternion
     double y = rot_axis.y*sin(angle/2.); //imaginary j part of quaternion
     double z = rot_axis.z*sin(angle/2.); //imaginary k part of quaternion

/** We use the rotation matrix defined by the quaternion components ($re$, $x$ ,$y$ $z$) */
     n.x = re*re*p1.x + 2*y*re*p1.z - 2*z*re*p1.y + x*x*p1.x + 2*y*x*p1.y + 2*z*x*p1.z - z*z*p1.x - y*y*p1.x;
     n.y = 2*x*y*p1.x + y*y*p1.y + 2*z*y*p1.z + 2*re*z*p1.x - z*z*p1.y + re*re*p1.y - 2*x*re*p1.z - x*x*p1.y;
     n.z = 2*x*z*p1.x + 2*y*z*p1.y + z*z*p1.z - 2*re*y*p1.x - y*y*p1.z + 2*re*x*p1.y - x*x*p1.z + re*re*p1.z;
  }

  return n;
}

/**
This function is an adaptation of the
[reconstruction()](/src/fractions.h#reconstruction) function which
takes into account the contact angle on embedded boundaries. */

void reconstruction_contact (scalar f, vector n, scalar alpha)
{

/**
We first reconstruct the (n, alpha) fields everywhere, using the
standard function. */

  reconstruction (f, n, alpha);
  vector gradf[];
  gradients ({f}, {gradf});

/**
In cells which contain an embedded boundary and an interface, we
modify the reconstruction to take the contact angle into account. */

  foreach(){
  	tag[]=0.;
    //if (cs[] < 1. && cs[] > 0. && f[] < 1.  && f[] > 0.) { //standard version
    //if ( interfacial(point, cs) && interfacial(point, cs) && cs[]>0. ) {
  	if ( interfacial(point, cs) && f[] < 1.  && f[] > 0. && cs[]>0. ) {
  	  tag[]=1;
	    coord ns = facet_normal (point, cs, fs);
	    normalize (&ns);
	    coord nf;
	    foreach_dimension()
	    nf.x = n.x[];

	    coord nc = normal_contact (ns, nf, contact_angle[]);
	    foreach_dimension()
	    n.x[] = nc.x;
	    alpha[] = line_alpha (f[], nc);
  	}
	}
  boundary ({n, alpha});  

}

/**
At every timestep, we modify the volume fraction values in the
embedded solid to take the contact angle into account. */

event contact (i++)
{
  vector n[];
  scalar alpha[];

/**
We first reconstruct (n,alpha) everywhere. Note that this is
necessary since we will use neighborhood stencils whose values may
be reconstructed by adaptive and/or parallel boundary conditions. */

  reconstruction_contact (f, n, alpha);

/**
We then look for "contact" cells in the neighborhood of each cell
entirely contained in the embedded solid. */

  foreach() {
    if (cs[] == 0.) {
      double fc = 0., sfc = 0.;
      coord o = {x, y, z};
      foreach_neighbor()
        //if (cs[] < 1. && cs[] > 0. && f[] < 1.  && f[] > 0.) {	// standard version  
        if (tag[]) {
/**
This is a contact cell. We compute a coefficient meant to
weight the estimated volume fraction according to how well
the contact point is defined. We assume here that a contact
point is better reconstructed if the values of `cs` and `f`
are both close to 0.5. */
          double coef = cs[]*(1. - cs[])*f[]*(1. - f[]);  
          if (coef == 0) {	// Case where the embedded boundary is strictly aligned with the mesh (coef=0)
          double mean = (f[]+cs[])/2.; // Approximation of the solid interface position to avoid division by zero
          coef= mean*(1. - mean)*f[]*(1. - f[]);
        	}
          sfc += coef;
/**
We then compute the volume fraction of the solid cell
(centered on `o`), using the extrapolation of the interface
reconstructed in the contact cell. */

          coord nf;
          foreach_dimension()
          nf.x = n.x[];
          coord a = {x, y, z}, b;
          foreach_dimension()
          a.x = (o.x - a.x)/Delta - 0.5, b.x = a.x + 1.;
          	fc += coef*rectangle_fraction (nf, alpha[], a, b);
        }
/**
The new volume fraction value of the solid cell is the weighted
average of the volume fractions reconstructed from all
neighboring contact cells. */

      if (sfc > 0.){
        f[] = fc/sfc;
      }
    }

  }
  boundary ({f});
}

// fixme: adaptive mesh 


