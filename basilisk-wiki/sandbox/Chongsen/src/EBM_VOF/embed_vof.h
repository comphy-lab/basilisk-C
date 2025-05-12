attribute {
  scalar * tracers, c;
  bool inverse;
}

#include "tmp_fraction_field.h"

extern scalar * interfaces;
extern scalar * tmp_interfaces;
extern scalar * interfaces_mark;
extern scalar * interfaces_contact_angle;
extern face vector uf;
extern double dt;


foreach_dimension()
static double vof_concentration_gradient_x (Point point, scalar c, scalar t)
{
  static const double cmin = 0.5;
  double cl = c[-1], cc = c[], cr = c[1];
  if (t.inverse)
    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  if (cc >= cmin && t.gradient != zero) {
    if (cr >= cmin) {
      if (cl >= cmin) {
	if (t.gradient)
	  return t.gradient (t[-1]/cl, t[]/cc, t[1]/cr)/Delta;
	else
	  return (t[1]/cr - t[-1]/cl)/(2.*Delta);
      }
      else
	return (t[1]/cr - t[]/cc)/Delta;
    }
    else if (cl >= cmin)
      return (t[]/cc - t[-1]/cl)/Delta;
  }
  return 0.;
}

#if TREE
static void vof_concentration_refine (Point point, scalar s)
{
  scalar f = s.c;
  if (cm[] == 0. || (!s.inverse && f[] <= 0.) || (s.inverse && f[] >= 1.))
    foreach_child()
      s[] = 0.;
  else {
    coord g;
    foreach_dimension()
      g.x = Delta*vof_concentration_gradient_x (point, f, s);
    double sc = s.inverse ? s[]/(1. - f[]) : s[]/f[], cmc = 4.*cm[];
    foreach_child() {
      s[] = sc;
      foreach_dimension()
	s[] += child.x*g.x*cm[-child.x]/cmc;
      s[] *= s.inverse ? 1. - f[] : f[];
    }
  }
}

event defaults (i = 0)
{

  for (scalar c in interfaces) {
    c.refine = c.prolongation = fraction_refine;
    c.dirty = true;
    scalar * tracers = c.tracers;
    for (scalar t in tracers) {
      t.restriction = restriction_volume_average;
      t.refine = t.prolongation = vof_concentration_refine;
      t.dirty = true;
      t.c = c;
    }
  }


}
#endif // TREE


event defaults (i = 0)
{

  for (scalar c in interfaces) {
    scalar * tracers = c.tracers;
    for (scalar t in tracers)
      t.depends = list_add (t.depends, c);
  }

}


event stability (i++) {

  double vmin = 1;
  for (scalar f in interfaces) {
    foreach(reduction(min:vmin)){
      foreach_neighbor(1){
        if (cs[]>0. && cs[]<1. && f[]>0.)
          vmin = min(vmin,cs[]);
      }
    }
  }

  if (CFL > 0.5*vmin)
    CFL = 0.5*vmin;
  
}


foreach_dimension()
static void sweep_x (scalar c, scalar tmp_c, scalar mark, scalar contact_angle, scalar cc, int dime, scalar * tcl)
{ 


  scalar flux[];

  double cfl = 0.;

  scalar * tracers = c.tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    for (scalar t in tracers) {
      scalar gf = new scalar, flux = new scalar;
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }

    foreach() {
      scalar t, gf;
      for (t,gf in tracers,gfl)
	gf[] = vof_concentration_gradient_x (point, c, t);
    }
  }

    vector ms[], tmp_mc[], mc[];

    scalar alphacs[], nalphac[];
    
    sort_cell(c, cs, contact_angle, mark);

    reconstruction_cs(cs, fs, ms, alphacs);

    reconstruction_mc(c, tmp_mc);
    
    reconstruction_tmp_embed_fraction_field (c, cs, tmp_mc, ms, alphacs, contact_angle, tmp_c, mark);

    boundary({tmp_c});

    
    reconstruction_test (c, tmp_c, mc, nalphac);

    foreach_face(x, reduction (max:cfl)) {

      double s = sign(uf.x[]), tfm = 1.;
      int i = -(s + 1.)/2.;
      
      
      if ((cs[i]>0. && cs[i]<1. && fs.x[]>0)){
        double tun = fabs(uf.x[]*dt/(Delta + SEPS));

        coord pp1[5]={{10,10},{10,10},{10,10},{10,10},{10,10}};
        coord p_mof1[2]={{10,10},{10,10}};

        double alphac1 = 0;
        if (fabs(tun)<1e-14)
          alphac1 = -0.5+tun;
        else {
          alphac1 = polygon_alpha (tun, (coord){1, 0}, (coord){-s*ms.x[i], ms.y[i]}, alphacs[i], p_mof1, pp1);
        }


        double tfmun = fabs(-0.5-alphac1);

        tfm = tun/(tfmun+SEPS);

      }

      
      double un = uf.x[]*dt/(Delta*tfm + SEPS);

      if (cs[] >= 1.){
        if (un*fm.x[]*s/(cm[] + SEPS) > cfl)
          cfl = un*fm.x[]*s/(cm[] + SEPS);
      }


      coord pp[5]={{10,10},{10,10},{10,10},{10,10},{10,10}};
      coord p_mof[2]={{10,10},{10,10}};

      double cf=0;

      if (cs[i] >= 1){
      
            cf = (c[i] <= 0. || c[i] >= 1.) ? c[i] :
        rectangle_fraction ((coord){-s*mc.x[i], mc.y[i], mc.z[i]}, nalphac[i], (coord){-0.5, -0.5, -0.5}, (coord){s*un - 0.5, 0.5, 0.5}); 
        
        }

      else if (cs[i] < 1. && cs[i] > 0.){  //  embed cells
      
          if (c[i] <= 0.)    
            cf = 0;

          
          else if (c[i] >= cs[i]){
              if (un != 0){
                coord tms={ms.x[i], ms.y[i]};

                double alphac = plane_alpha (c[i], tms);
                
                      cf = (c[i] <= 0. || c[i] >= 1.) ? c[i] :
                  rectangle_fraction ((coord){-s*ms.x[i], ms.y[i], ms.z[i]}, alphac, (coord){-0.5, -0.5, -0.5}, (coord){s*un - 0.5, 0.5, 0.5}); 
              }
              else
                cf = 0;
          }


          else if (c[i]>0. && c[i]<cs[i]){

              if (un != 0){
                
                if (mark[i]==4 || mark[i]==5){
                  double alphac = polygon_alpha (c[i], (coord){mc.x[i], mc.y[i], mc.z[i]}, (coord){ms.x[i], ms.y[i], ms.z[i]}, alphacs[i], p_mof, pp);
                
                  cf = polygon_fraction(pp, un, s, p_mof, (coord){mc.x[i], mc.y[i], mc.z[i]}, alphac)/fabs((un));
                }
                else if (mark[i]==3){
                  coord tms={ms.x[i], ms.y[i]};

                  double alphac = plane_alpha (c[i], tms);

                      cf = (c[i] <= 0. || c[i] >= 1.) ? c[i] :
                  rectangle_fraction ((coord){-s*ms.x[i], ms.y[i], ms.z[i]}, alphac, (coord){-0.5, -0.5, -0.5}, (coord){s*un - 0.5, 0.5, 0.5}); 
                }
                else if (mark[i]==2){
                  cf = 0;
                }

                }

              else
              cf = 0;

          }
                
      }
      else
        cf = 0;

      flux[] = cf*uf.x[]/(tfm+SEPS);

        scalar t, gf, tflux;
    for (t,gf,tflux in tracers,gfl,tfluxl) {
      double cf1 = cf, ci = c[i];
      if (t.inverse)
	  cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
	  double ff = t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2.;
	  tflux[] = ff*cf1*uf.x[];
      }
      else
   	tflux[] = 0.;
    }
    
  }
  
    
      delete (gfl); free (gfl);

  if (cfl > 0.5 + 1e-6)
  fprintf (ferr, 
	     "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n", 
	     cfl - 0.5), fflush (ferr);

  foreach(){
    if (cs[] > 0.) {
      
      c[] += dt*(flux[] - flux[1] + cc[]*(uf.x[1] - uf.x[]))/Delta;
      
      scalar t, tc, tflux;
      for (t, tc, tflux in tracers, tcl, tfluxl)

      t[] += dt*(tflux[] - tflux[1] + tc[]*(uf.x[1] - uf.x[]))/Delta;
    }
    
    
  }


  delete (tfluxl); free (tfluxl);

}

/**
## Multi-dimensional advection

The multi-dimensional advection is performed by the event below. */
trace
void vof_advection (scalar * interfaces, scalar * tmp_interfaces, scalar * interfaces_mark, scalar * interfaces_contact_angle, int i)
{

  for (scalar contact_angle in interfaces_contact_angle)
  for (scalar mark in interfaces_mark)
  for (scalar c in interfaces) 
  for (scalar tmp_c in tmp_interfaces) {

    scalar cc[], * tcl = NULL, * tracers = c.tracers;    
    for (scalar t in tracers) {
      scalar tc = new scalar;
      tcl = list_append (tcl, tc);
#if TREE
      if (t.refine != vof_concentration_refine) {
	t.refine = t.prolongation = vof_concentration_refine;
	t.restriction = restriction_volume_average;
	t.dirty = true;
	t.c = c;
      }
#endif // TREE
    }

    foreach() {

      double vmin = cs[];

      cc[] = (c[] > 0.5*vmin);

      scalar t, tc;
      for (t, tc in tracers, tcl) {
	if (t.inverse)
	  tc[] = c[] < 0.5 ? t[]/(1. - c[]) : 0.;
	else
	  tc[] = c[] > 0.5 ? t[]/c[] : 0.;
      }
    }
    /**
    We then apply the one-dimensional advection scheme along each
    dimension. To try to minimise phase errors, we alternate dimensions
    according to the parity of the iteration index `i`. */

    void (* sweep[dimension]) (scalar, scalar, scalar, scalar, scalar, int, scalar *);
    int d = 0;
    foreach_dimension()
      sweep[d++] = sweep_x;
    for (d = 0; d < dimension; d++){

      sweep[(i + d) % dimension] (c, tmp_c, mark, contact_angle, cc, (i + d) % dimension, tcl);

    }


    reconstruction_tmp_c(c, contact_angle, cs, fs, tmp_c, mark);
    delete (tcl), free (tcl);
    
  }
  
}

event vof (i++){


  vof_advection (interfaces, tmp_interfaces, interfaces_mark, interfaces_contact_angle, i) ; 


}
/**
## References

~~~bib
~~~
*/