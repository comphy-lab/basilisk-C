/**
# EHD stresses (alternative scheme) */

event defaults (i = 0) {
  if (is_constant(a.x))
    a = new face vector;
}

/**
Center values of the permittivity are needed. For the moment we
have decided to define it outside of the library. Also,
note that the metric factors are not considered in *epsilonc[]*. */

//scalar epsilonc[];

event acceleration (i++) {
  symmetric tensor M[];

  foreach()
    foreach_dimension()
       M.x.x[] = epsilonc[]/2.*cm[]*
       (sq(phi[1] - phi[-1]) - sq(phi[0,1] - phi[0,-1]))/sq(2.*Delta);

  foreach_vertex()
    M.x.y[] =
    (epsilonc[] + epsilonc[-1] + epsilonc[0,-1] + epsilonc[-1,-1])/4.*
    (cm[] + cm[-1] + cm[0,-1] + cm[-1,-1])/4.*
    (phi[] - phi[0,-1] + phi[-1] - phi[-1,-1])*
    (phi[] - phi[-1] + phi[0,-1] - phi[-1,-1])/sq(2.*Delta);

  boundary ({M.x.x, M.y.y});
  
  face vector av = a;
  foreach_face()
    av.x[] += (fm.x[] < 1.1e-20 ? 0 : alpha.x[]/fm.x[]* 
	       (M.x.x[] - M.x.x[-1,0] +
		M.x.y[0,1] - M.x.y[])/(fm.x[]*Delta));

#if AXI
  /**
  There is a round-off problem at the axis of revolution. We set the
  acceleration directly to zero. */

  foreach_face(y)
    av.y[] += (y == 0 ? 0 : (epsilonc[] + epsilonc[0,-1])/4.
               *alpha.y[]/fm.y[]*
               (sq(phi[] - phi[0,-1]) + 
                sq(phi[1,0] + phi[1,-1] -
                   phi[-1,0]-phi[-1,-1])/16.)/(fm.y[]*sq(Delta)));
#endif
}
