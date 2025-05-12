#if AXI
# define gamma (1./y)
#else
# define gamma (1.)
#endif

void dissipation (scalar dis, vector u,  (const) face vector mu
#if COMPRESSIBLE
  , (const) face vector lambdav
#endif
)
{

#if AXI
  scalar ur = u.y;
#endif
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */

  foreach () {
#if AXI
    dis[] = ((mu.x[] + mu.x[1] + mu.y[] + mu.y[0,1])/2.*gamma*sq(ur[]/y)
#if COMPRESSIBLE
	     + (lambdav.x[] + lambdav.x[1] +lambdav.y[]+lambdav.y[0,1])/4.*
	     (1 + (u.x[1] -u.x[-1] + u.y[0,1] -u.y[0,-1])/2./Delta))*ur[]/y
#endif
	     );
#else
    dis[] = 0.;
#endif
  }

  foreach_dimension() {
    face vector taux[];
#if COMPRESSIBLE
    face vector tauc[];
#if AXI
    face vector axic[];
#endif
#endif
    foreach_face(x) {
#if COMPRESSIBLE
#if AXI
      axic.x[] = lambdav.x[]*(ur[]+ur[-1])
	*(u.x[] - u.x[-1])/Delta/2.;
#endif
      tauc.x[] = lambdav.x[]*(u.x[] - u.x[-1]
     #if dimension > 1
			      + (u.y[0,1] + u.y[-1,1])/4
			      - (u.y[0,-1] + u.y[-1,-1])/4.
    #endif
    #if dimension > 2
			      + (u.z[0,0,1] + u.z[-1,0,1])/4
			      - (u.z[0,0,-1] + u.z[-1,0,-1])/4.
    #endif
			      )*(u.x[] - u.x[-1])/sq(Delta);
#endif
      taux.x[] =  2.*mu.x[]*sq((u.x[] - u.x[-1])/Delta);
    }
    #if dimension > 1
      foreach_face(y)
	taux.y[] = mu.y[]*(u.x[] - u.x[0,-1] +
			   (u.y[1,-1] + u.y[1,0])/4. -
			   (u.y[-1,-1] + u.y[-1,0])/4.)
	*(u.x[] - u.x[0,-1])/sq(Delta);
    #endif
    #if dimension > 2
      foreach_face(z)
	taux.z[] = mu.z[]*(u.x[] - u.x[0,0,-1] +
			   (u.z[1,0,-1] + u.z[1,0,0])/4. -
			   (u.z[-1,0,-1] + u.z[-1,0,0])/4.)
	*(u.x[] - u.x[0,0,-1])/sq(Delta);
    #endif  
      foreach () {
        double d = 0;
	foreach_dimension()
	  d += (taux.x[1] + taux.x[])*gamma;
	dis[] += (d
#if COMPRESSIBLE 
               + tauc.x[1] + tauc.x[]
#if AXI
               + (axic.x[1] + axic.x[])/y
#endif
#endif
	       )/2.;
      }
  }
#else
/*   /\* "naive" discretisation (only 1st order on trees) *\/ */
  foreach () {
#if AXI
    dis[] = ((mu.x[] + mu.x[1] + mu.y[] + mu.y[0,1])/2.*gamma*sq(ur[]/y)
#if COMPRESSIBLE
	     + (lambdav.x[] + lambdav.x[1] +lambdav.y[]+lambdav.y[0,1])/4.*
	     (ur[]/y + (u.x[1] -u.x[-1] + u.y[0,1] -u.y[0,-1])/2./Delta)*ur[]/y
#endif
	     );
#else
    dis[] = 0.;
#endif

    foreach_dimension() 
	     dis[] += ((mu.x[1,0]*sq(u.x[1] - u.x[])
		+ mu.x[]*sq(u.x[] - u.x[-1])
      #if dimension > 1
		+ mu.y[0,1]*(u.x[0,1] - u.x[] +
		(u.y[1,0] + u.y[1,1])/4. -
		(u.y[-1,0] + u.y[-1,1])/4.)*(u.x[0,1] - u.x[])/2.
		+ mu.y[]*(u.x[] - u.x[0,-1] +
		(u.y[1,-1] + u.y[1,0])/4. -
		(u.y[-1,-1] + u.y[-1,0])/4.)*(u.x[] - u.x[0,-1])/2.
      #endif
      #if dimension > 2
		+ mu.z[0,0,1]*(u.x[0,0,1] - u.x[] +
			       (u.z[1,0,0] + u.z[1,0,1])/4. -
			       (u.z[-1,0,0] + u.z[-1,0,1])/4.)*(u.x[0,0,1] - u.x[])/2.
		- mu.z[]*(u.x[] - u.x[0,0,-1] +
			  (u.z[1,0,-1] + u.z[1,0,0])/4. -
			  (u.z[-1,0,-1] + u.z[-1,0,0])/4.)*(u.x[] - u.x[0,0,-1])/2.
       #endif
		)*gamma
#if COMPRESSIBLE
	     + lambdav.x[1]*sq(u.x[1] - u.x[])/2.
	     + lambdav.x[]*sq(u.x[] - u.x[-1])/2.
               #if dimension > 1
	     + lambdav.x[1]*((u.y[1,1] + u.y[0,1])/4 -
	                     (u.y[1,-1] + u.y[0,-1])/4.)*(u.x[1] - u.x[])/2.
	     + lambdav.x[]*((u.y[0,1] + u.y[-1,1])/4 -
			    (u.y[0,-1] + u.y[-1,-1])/4.)*(u.x[] - u.x[-1])/2.
#if AXI
	     + lambdav.x[1]*(ur[1] + ur[])*(u.x[1] - u.x[])/4./y*Delta
	     + lambdav.x[]*(ur[-1] + ur[])*(u.x[] - u.x[-1])/4./y*Delta
#endif
               #endif
               #if dimension > 2
		  + lambdav.x[1]*((u.z[1,0,1] + u.z[0,0,1])/4 -
				  (u.z[1,0,-1] + u.z[0,0,-1])/4.)*(u.x[1] - u.x[])/2.
		  + lambdav.x[]*((u.z[0,0,1] + u.z[-1,0,1])/4 -
				    (u.z[0,0,-1] + u.z[-1,0,-1])/4.)*(u.x[] - u.x[-1])/2.
              #endif
#endif
	     )/sq(Delta);
}
#endif
}	     
  
#undef gamma
