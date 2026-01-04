#include "fractions.h"
#include "geometry_ebm.h"

/** 
 ### Bug fix for fraction function (2D case) */

trace
void fractions_ebm (vertex scalar Phi, scalar c,
		face vector s = {0}, double val = 0.)
{
#if dimension > 1
  face vector as = automatic (s);
  
  /**
  We store the positions of the intersections of the surface with the
  edges of the cell in vector field `p`. In two dimensions, this field
  is just the transpose of the *line fractions* `s`, in 3D we need to
  allocate a new field. */
  
#if dimension == 3
  vector p[];
#else // dimension == 2
  vector p;
  p.x = as.y; p.y = as.x;
#endif
  
  /**
  ### Line fraction computation
  
  We start by computing the *line fractions* i.e. the (normalised)
  lengths of the edges of the cell within the surface. */

  foreach_edge() {

    /**
    If the values of $\Phi$ on the vertices of the edge have opposite
    signs, we know that the edge is cut by the interface. */

    if ((Phi[] - val)*(Phi[1] - val) < 0.) {

      /**
      In that case we can find an approximation of the interface position by
      simple linear interpolation. We also check the sign of one of the
      vertices to orient the interface properly. */

      p.x[] = (Phi[] - val)/(Phi[] - Phi[1]);
      if (Phi[] < val)
	p.x[] = 1. - p.x[];
    }

    /**
    If the values of $\Phi$ on the vertices of the edge have the same sign
    (or are zero), then the edge is either entirely outside or entirely
    inside the interface. We check the sign of both vertices to treat
    limit cases properly (when the interface intersects the edge exactly
    on one of the vertices). */

    else
      p.x[] = (Phi[] > val || Phi[1] > val);
  }

  /**
  ### Surface fraction computation 

  We can now compute the surface fractions. In 3D they will be
  computed for each face (in the z, x and y directions) and stored in
  the face field `s`. In 2D the surface fraction in the z-direction is
  the *volume fraction* `c`. */

#if dimension == 3

  /**
  In 3D we need to prevent boundary conditions, since this would
  impose vertex field BCs which are not (apparently) consistent for
  the edge intersection coordinates. This can probably be improved. */
  
  foreach_dimension()
    p.x.dirty = false;
  
  scalar s_x = as.x, s_y = as.y, s_z = as.z;
  foreach_face(z,x,y)
#else // dimension == 2
  scalar s_z = c;
  foreach()
#endif
  {

    /**
    We first compute the normal to the interface. This can be done easily
    using the line fractions. The idea is to compute the circulation of
    the normal along the boundary $\partial\Omega$ of the fraction of the
    cell $\Omega$ inside the interface. Since this is a closed curve, we
    have
    $$
    \oint_{\partial\Omega}\mathbf{n}\;dl = 0
    $$ 
    We can further decompose the integral into its parts along the edges
    of the square and the part along the interface. For the case pictured
    above, we get for one component (and similarly for the other)
    $$
    - s_x[] + \oint_{\Phi=0}n_x\;dl = 0
    $$
    If we now define the *average normal* to the interface as
    $$
    \overline{\mathbf{n}} = \oint_{\Phi=0}\mathbf{n}\;dl
    $$
    We have in the general case
    $$
    \overline{\mathbf{n}}_x = s_x[] - s_x[1,0]
    $$
    and
    $$
    |\overline{\mathbf{n}}| = \oint_{\Phi=0}\;dl
    $$ 
    Note also that this average normal is exact in the case of a linear
    interface. */

    coord n;
    double nn = 0.;
    foreach_dimension(2) {
      n.x = p.y[] - p.y[1];
      nn += fabs(n.x);
    }
    
    /**
    If the norm is zero, the cell is full or empty and the surface fraction
    is identical to one of the line fractions. */

    if (nn == 0.)
      s_z[] = p.x[];
    else {
    
      /**
      Otherwise we are in a cell containing the interface. We first
      normalise the normal. */

      foreach_dimension(2)
	n.x /= nn;

      /**
      To find the intercept $\alpha$, we look for edges which are cut by the
      interface, find the coordinate $a$ of the intersection and use it to
      derive $\alpha$. We take the average of $\alpha$ for all intersections. */
      
      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
	foreach_dimension(2)
	  if (p.x[0,i] > 0. && p.x[0,i] < 1.) {
	    double a = sign(Phi[0,i] - val)*(p.x[0,i] - 0.5);
	    alpha += n.x*a + n.y*(i - 0.5);
	    ni++;
	  }
    else if (p.x[0,i] == 0 && p.y[i,0] == 1) { //for interface crossing grid corner
      alpha += n.x*(i - 0.5) + n.y*(i - 0.5);
      ni++;
    }
    else if (p.x[0,i] == 0 && p.y[1-i,0] == 1) {
      alpha += n.x*(0.5 - i) + n.y*(i - 0.5);
      ni++;
    }

      /**
      Once we have $\mathbf{n}$ and $\alpha$, the (linear) interface
      is fully defined and we can compute the surface fraction using
      our pre-defined function. For marginal cases, the cell is full
      or empty (*ni == 0*) and we look at the line fractions to
      decide. */

      if (ni == 0)
	s_z[] = max (p.x[], p.y[]);
      else if (ni != 4)
	s_z[] = line_area (n.x, n.y, alpha/ni);
      else {
#if dimension == 3
	s_z[] = (p.x[] + p.x[0,1] + p.y[] + p.y[1] > 2.);
#else
	s_z[] = 0.;
#endif
      }
    }
  }
  
  /**
  ### Volume fraction computation

  To compute the volume fraction in 3D, we use the same approach. */
  
#if dimension == 3
  foreach() {

    /**
    Estimation of the average normal from the surface fractions. */
       
    coord n;
    double nn = 0.;
    foreach_dimension(3) {
      n.x = as.x[] - as.x[1];
      nn += fabs(n.x);
    }
    if (nn == 0.)
      c[] = as.x[];
    else {
      foreach_dimension(3)
	n.x /= nn;

      /**
      We compute the average value of *alpha* by looking at the
      intersections of the surface with the twelve edges of the
      cube. */
      
      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
	for (int j = 0; j <= 1; j++)
	  foreach_dimension(3)
	    if (p.x[0,i,j] > 0. && p.x[0,i,j] < 1.) {
	      double a = sign(Phi[0,i,j] - val)*(p.x[0,i,j] - 0.5);
	      alpha += n.x*a + n.y*(i - 0.5) + n.z*(j - 0.5);
	      ni++;
	    }

      /**
      Finally we compute the volume fraction. */

      if (ni == 0)
	c[] = as.x[];
      else if (ni < 3 || ni > 6)
	c[] = 0.; // this is important for robustness of embedded boundaries
      else
	c[] = plane_volume (n, alpha/ni);
    }
  }
#endif // dimension == 3
#else  // dimension == 1
  if (s.x.i)
    foreach_face()
      s.x[] = Phi[] > 0.;
  foreach()
    if ((Phi[] - val)*(Phi[1] - val) < 0.) {
      c[] = (Phi[] - val)/(Phi[] - Phi[1]);
      if (Phi[] < val)
	c[] = 1. - c[];
    }
    else
      c[] = (Phi[] > val || Phi[1] > val);
#endif
}

macro fraction_ebm (scalar f, double func)
{
  {
    vertex scalar phi[];
    foreach_vertex()
      phi[] = func;
    fractions_ebm (phi, f);
  }
}

macro solid_ebm (scalar cs, face vector fs, double func)
{
  {
    vertex scalar phi[];
    foreach_vertex()
      phi[] = func;
    fractions_ebm (phi, cs, fs);
  }
}

//
void reconstruction_cs (const scalar cs, face vector fs, vector n, scalar alpha)
{
  foreach() {
    if (cs[] <= 0. || cs[] >= 1.) {
      alpha[] = 0.;
      foreach_dimension()
      	n.x[] = 0.;
    }
    else {
      coord m = facet_normal (point, cs, fs);
      foreach_dimension()
	      n.x[] = m.x;
      alpha[] = plane_alpha (cs[], m);
    }
  }

#if TREE
  foreach_dimension()
    n.x.refine = n.x.prolongation = refine_injection;
  alpha.n = n;
  alpha.refine = alpha.prolongation = alpha_refine;
#endif
}

void output_facets_ebm (scalar c, scalar cs, scalar cet, FILE * fp = stdout, face vector s = {{-1}})
{
  foreach()
    if (c[] > 0. && c[] < cs[]) {
      coord n = interface_normal (point, c);
      double alpha = plane_alpha (cet[], n);
#if dimension == 1
      fprintf (fp, "%g\n", x + Delta*alpha/n.x);
#elif dimension == 2
      coord segment[2];
      if (facets_ebm (n, alpha, segment) == 2)
	fprintf (fp, "%g %g\n%g %g\n\n", 
		 x + segment[0].x*Delta, y + segment[0].y*Delta, 
		 x + segment[1].x*Delta, y + segment[1].y*Delta);
#else // dimension == 3
      coord v[12];
      int m = facets_ebm (n, alpha, v, 1.);
      for (int i = 0; i < m; i++)
	fprintf (fp, "%g %g %g\n",
		 x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
      if (m > 0)
	fputc ('\n', fp);
#endif
    }

  fflush (fp);
}

#if dimension == 2

void intersect_vof_solid (scalar c, scalar cs, face vector fs)
{
  scalar cfg = c.cfg;
  scalar cet = c.cet;
  if (iter == 0)
    foreach() {
      cfg[] = nodata;
      cet[] = nodata;
    }

  scalar cnew[];
  foreach()
    if (cs[] == 1.)
      cnew[] = c[];
    else if (cs[] == 0.)
      cnew[] = 0.;
    else { // 0 < cs < 1
      if (c[] == 0.) {
        cnew[] = 0.;
        if (iter == 0)
          cfg[] = -2.;
      }
      else if (c[] == 1.) {
        cnew[] = cs[];
        if (iter == 0)
          cfg[] = 2.;
      }
      else { // 0 < c < 1
        coord c_n = mycs (point, c);
        double c_alpha = plane_alpha (c[], c_n);
        coord segment_c[2];
        assert (facets_ebm (c_n, c_alpha, segment_c) == 2);

        coord cs_n = facet_normal (point, cs, fs);
        double cs_alpha = plane_alpha (cs[], cs_n);
        coord segment_cs[2];
        assert (facets_ebm (cs_n, cs_alpha, segment_cs) == 2);

        double sign_c1 = cs_alpha - cs_n.x*segment_c[0].x - cs_n.y*segment_c[0].y;
        double sign_c2 = cs_alpha - cs_n.x*segment_c[1].x - cs_n.y*segment_c[1].y;
        double sign_cs1 = c_alpha - c_n.x*segment_cs[0].x - c_n.y*segment_cs[0].y;
        double sign_cs2 = c_alpha - c_n.x*segment_cs[1].x - c_n.y*segment_cs[1].y;

        if (sign_cs1 <= 0. && sign_cs2 <= 0.) {
          if (sign_c1 >= 0. && sign_c2 >= 0.)
            cnew[] = c[];
          else {
            assert (sign_c1 <= 0. && sign_c2 <= 0.);
            cnew[] = 0.;
          }
          if (iter == 0)
            cfg[] = -2.;
        }
        else if (sign_cs1 >= 0. && sign_cs2 >= 0.) {
          if (sign_c1 >= 0. && sign_c2 >= 0.)
            cnew[] = c[] + cs[] - 1.;
          else {
            assert (sign_c1 <= 0. && sign_c2 <= 0.);
            cnew[] = cs[];
          }
          if (iter == 0)
            cfg[] = 2.;
        }
        else {
          assert ((sign_cs1 > 0. && sign_cs2 < 0.) || (sign_cs1 < 0. && sign_cs2 > 0.));
          assert ((sign_c1 > 0. && sign_c2 < 0.) || (sign_c1 < 0. && sign_c2 > 0.));
          coord start, end, cross; 
          coord cornor[4];

          if (sign_cs1 > 0)
            start = segment_cs[0];
          else
            start = segment_cs[1];
          
          if (sign_c1 < 0)
            end = segment_c[0];
          else
            end = segment_c[1];

          int i = 0;
          for (double sx = -0.5; sx <= 0.5; sx += 1.)
            for (double sy = -0.5; sy <= 0.5; sy += 1.)
              if (cs_alpha - cs_n.x*sx - cs_n.y*sy < 0 &&
                   c_alpha -  c_n.x*sx -  c_n.y*sy > 0) {
                cornor[i].x = sx;
                cornor[i].y = sy;
                i++;
              }
          assert (i < 4);

          assert (fabs(cs_n.y*c_n.x - cs_n.x*c_n.y) > 0);
          cross.x = (c_alpha*cs_n.y - cs_alpha*c_n.y)/(cs_n.y*c_n.x - cs_n.x*c_n.y);
          assert (fabs(c_n.y) > 0 || fabs(cs_n.y) > 0);
          if (fabs(c_n.y) > fabs(cs_n.y))
            cross.y = (c_alpha - c_n.x*cross.x)/c_n.y;
          else
            cross.y = (cs_alpha - cs_n.x*cross.x)/cs_n.y;

          coord polygon[7] = {{nodata}, {nodata}, {nodata}, {nodata}, {nodata}, {nodata}, {nodata}};
          polygon[0] = start;

          int num = 0;
          for (int k = 0; k < i; k++)
            for (int j = 0; j < i-k; j++)
              if (cornor[j].x == polygon[k].x || cornor[j].y == polygon[k].y) {
                polygon[k+1] = cornor[j];
                int m = j;
                while (m+1 < i-k) {
                  cornor[m] = cornor[m+1];
                  m++;
                }
                num++;
                break;
              }
          assert (num == i);

          polygon[i+1] = end;
          polygon[i+2] = cross;
          polygon[i+3] = start;

          double sum = 0;
          for (int k = 0; k < i+3; k++) {
            assert (polygon[k].x != nodata);
            sum += polygon[k].x*polygon[k+1].y - polygon[k].y*polygon[k+1].x;
          }
          double area = fabs(sum)/2.;
          cnew[] = c[] - area;

          if (iter == 0)
            cfg[] = 1.;
        }
      }
    }
  foreach()
    c[] = cnew[] < 1e-10 ? 0. : cnew[] > cs[] - 1e-10 ? cs[] : cnew[];
}

void reconstruction_emd (const scalar c, vector n, scalar alpha, const scalar cs, vector ncs, scalar alphacs)
{
  foreach() {
    assert (c[] >= 0. && c[] <= cs[]);
    if (c[] <= 0. || c[] >= cs[]) {
      alpha[] = 0.;
      foreach_dimension()
	      n.x[] = 0.;
    }
    else if (cs[] >= 1.) {
      coord m = interface_normal (point, c);
      foreach_dimension()
	      n.x[] = m.x;
      alpha[] = plane_alpha (c[], m);
    }
    else {
      coord m = interface_normal (point, c), mcs;
      foreach_dimension() {
	      n.x[] = m.x;
        mcs.x = ncs.x[];
      }
      alpha[] = line_alpha_ebm (c[], cs[], mcs, alphacs[], m);
    }
  }
}

static coord identify_cl_pos (coord c_n, double c_alpha, coord cs_n, double cs_alpha)
{
  coord cross;
  assert (fabs(cs_n.y*c_n.x - cs_n.x*c_n.y) > 0);
  cross.x = (c_alpha*cs_n.y - cs_alpha*c_n.y)/(cs_n.y*c_n.x - cs_n.x*c_n.y);
  assert (fabs(c_n.y) > 0 || fabs(cs_n.y) > 0);
  if (fabs(c_n.y) > fabs(cs_n.y))
    cross.y = (c_alpha - c_n.x*cross.x)/c_n.y;
  else
    cross.y = (cs_alpha - cs_n.x*cross.x)/cs_n.y;

  if (fabs(cross.x) <= 0.5 && fabs(cross.y) <= 0.5)
    return cross;
  else {
    coord segment_cs[2];
    assert (facets_ebm (cs_n, cs_alpha, segment_cs) == 2);
    double d[2];
    for (int k = 0; k < 2; k++)
      d[k] = fabs(segment_cs[k].x - cross.x) + fabs(segment_cs[k].y - cross.y);

    return (d[0] < d[1] ? segment_cs[0] : segment_cs[1]);
  }
}

void update_cet (const scalar c, const scalar cs, face vector fs, scalar cfg, scalar ceff, scalar angle)
{
  vector mcl = c.mcl;
  vector hp  = c.hp;
  foreach()
    foreach_dimension() {
      mcl.x[] = 0.;
      hp.x[]  = 0.;
    }

  foreach()
    ceff[] = nodata;

  foreach() {
    assert (c[] >= 0. && c[] <= cs[]);
    if (cs[] <= 0.)
      ceff[] = nodata;
    else if (c[] <= 0. || c[] >= cs[])
      ceff[] = (c[] <= 0. ? 0. : 1.);
    else if (cs[] >= 1.)
      ceff[] = c[];
    else {
      assert (cfg[] != nodata);
      if (cfg[] == 2.)
        ceff[] = c[] - cs[] + 1.;
      else if (cfg[] == -2.)
        ceff[] = c[];
      else if (cfg[] == 1.) {
        coord ncs = facet_normal (point, cs, fs);
        double alphacs = plane_alpha (cs[], ncs);
        
        coord n = mycs (point, c);
        coord ncl = normal_contact (ncs, n, angle[]); normalize (&ncl);
        foreach_dimension()
	        mcl.x[] = ncl.x;
        double alphacl = line_alpha_ebm (c[], cs[], ncs, alphacs, ncl);

        coord clp = identify_cl_pos (ncl, alphacl, ncs, alphacs);
        assert (fabs(clp.x) <= 0.5 && fabs(clp.y) <= 0.5);
        foreach_dimension()
	        hp.x[] = clp.x;

        ceff[] = plane_volume (ncl, alphacl);
      }
      else
        assert (0);
    }
  }
}

#endif