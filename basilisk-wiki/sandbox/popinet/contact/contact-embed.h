/**
# Contact angles on an embedded boundary

This header file implements contact angles for [VOF
interfaces](/src/vof.h) on [embedded boundaries](/src/embed.h).

The contact angle is defined by this field and can be constant or
variable. */

(const) scalar contact_angle;

/**
This function returns the properly-oriented normal of an interface
touching the embedded boundary. `ns` is the normal to the embedded
boundary, `nf` the VOF normal, not taking into account the contact
angle, and `angle` the contact angle. */

static inline coord normal_contact (coord ns, coord nf, double angle)
{
  assert (dimension == 2); // fixme: 2D only for the moment
  coord n;
  if (- ns.x*nf.y + ns.y*nf.x > 0) { // 2D cross product
    n.x = - ns.x*cos(angle) + ns.y*sin(angle);
    n.y = - ns.x*sin(angle) - ns.y*cos(angle);
  }
  else {
    n.x = - ns.x*cos(angle) - ns.y*sin(angle);
    n.y =   ns.x*sin(angle) - ns.y*cos(angle);
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

  /**
  In cells which contain an embedded boundary and an interface, we
  modify the reconstruction to take the contact angle into account. */
  
  foreach()
    if (cs[] < 1. && cs[] > 0. && f[] < 1.  && f[] > 0.) {
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
	if (cs[] < 1. && cs[] > 0. && f[] < 1.  && f[] > 0.) {	  

	  /**
	  This is a contact cell. We compute a coefficient meant to
	  weight the estimated volume fraction according to how well
	  the contact point is defined. We assume here that a contact
	  point is better reconstructed if the values of `cs` and `f`
	  are both close to 0.5. */

	  double coef = cs[]*(1. - cs[])*f[]*(1. - f[]);
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
      
      if (sfc > 0.)
	f[] = fc/sfc;
    }
  }
  boundary ({f});
}
