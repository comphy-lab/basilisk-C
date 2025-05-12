#ifndef thrsh
#define thrsh 0.15
#endif

#include "functionsHF.h"

extern scalar * interfaces;

attribute {
  scalar sigmaf;
}

event defaults (i = 0) {
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() {
      a.x[] = 0.;
      dimensional (a.x[] == Delta/sq(DT));
    }
  }
}

event stability (i++)
{
  double sigma = 0.;
  for (scalar d in interfaces)
    if (is_constant (d.sigmaf))
      sigma += constant (d.sigmaf);
  double sigmamax = sigma;
  
  double amin = HUGE, amax = -HUGE, dmin = HUGE;
  foreach_face (reduction(min:amin) reduction(max:amax) reduction(min:dmin)
		reduction(max:sigmamax))
    if (fm.x[] > 0.) {
      if (alpha.x[]/fm.x[] > amax) amax = alpha.x[]/fm.x[];
      if (alpha.x[]/fm.x[] < amin) amin = alpha.x[]/fm.x[];
      if (Delta < dmin) dmin = Delta;

      double sigmai = sigma;
      for (scalar d in interfaces)
	if (!is_constant (d.sigmaf) && d[] < 0.999999 && d[] > 0.000001) {
	  scalar sigmaf = d.sigmaf;
	  sigmai += sigmaf[];
	}
      if (sigmai > sigmamax)
	sigmamax = sigmai;
    }
  double rhom = (1./amin + 1./amax)/2.;

  if (sigmamax) {
    double dt = sqrt (rhom*cube(dmin)/(pi*sigmamax));
    if (dt < dtmax)
      dtmax = dt;
  }
}

#if AXI
# include "fractions.h"
#endif

vector nrm[];

foreach_dimension()
  double getdiagterm_y(Point point, scalar d,scalar sigma, double kappahf){
  double termyy=0.;
  for (int i = -1; i <= 1; i += 2)
    if (d[]*(d[] + d[i]) < 0.) {
      double xi = d[]/(d[] - d[i] + SEPS);
      double nx = 2.*xi*(nrm.x[] + nrm.x[i])/2. + (1. - 2.*xi)*(nrm.x[]);
      double sigmai = sigma[] + xi*(sigma[i] - sigma[]);
      double ki = kappahf;
      termyy += sigmai*(fabs(nx)/Delta - sign(d[])*ki*(0.5 - xi));
    }
  return termyy;
}

foreach_dimension()
  double getnondiagterm_y(Point point, scalar d,scalar sigma){
  double termxy = 0.;
  if ((d[] + d[0,-1])*(d[-1] + d[-1,-1]) > 0.)
    return termxy;
  else {
    double xi = (d[-1] + d[-1,-1])/(d[-1] + d[-1,-1] - d[] - d[0,-1]);
    double ny = (xi*(nrm.y[] + nrm.y[0,-1])/2. + (1. - xi)*(nrm.y[-1,0] + nrm.y[-1,-1])/2.);
    double sigmai = (sigma[-1] + sigma[-1,-1] +
		     xi*(sigma[] + sigma[0,-1] - sigma[-1] -  sigma[-1,-1]))/2.;
    termxy = - sigmai*sign(d[] + d[0,-1])*ny/Delta;
  }
  return termxy;
}

bool interfacialvertex(Point point, scalar f){

  for(int ii = -1; ii <=0; ii++)
    for(int jj = -1; jj <=0; jj++)
      if(f[ii,jj] > 0 && f[ii,jj] < 1.)
	return true;

  return false;	
}


scalar kappa[];
event acceleration (i++)
{

  for (scalar d in interfaces)
    if (d.sigmaf.i) {
      (const) scalar sigma = d.sigmaf;
      vector hf[];
      scalar dx[],dy[];
      scalar sd[];

      foreach(){
	d[] = clamp(d[],0.,1.);
	sd[] = (4.*d[] + 
		2.*(d[0,1] + d[0,-1] + d[1,0] + d[-1,0]) +
		d[-1,-1] + d[1,-1] + d[1,1] + d[-1,1])/16.;
      }
      
      heights(d,hf);
      /* curvature(sd,kappa); */
      
      foreach(){
	kappa[] =  height_curvature_fit(point,f,hf);
	//kappa[] = height_curvatureV2(point,sd,hf);
	coord nrmf = height_normal(point,f,hf);
	foreach_dimension()
	  nrm.x[] = -nrmf.x; // to match with d[1] - d[-1], I put -ve sign
      }

      foreach(){
	if(hf.x[] < nodata)
	  dx[] = hf.x[] > HSHIFT/2. ? HSHIFT - hf.x[] : hf.x[];
	else
	  dx[] = nodata;
	if(hf.y[] < nodata)
	  dy[] = hf.y[] > HSHIFT/2. ? HSHIFT - hf.y[] : hf.y[];
	else
	  dy[] = nodata;

      }
      boundary({dx,dy});
      
      scalar Sxx[], Sxy[], Syy[], Syx[];
      tensor S; S.x.x = Sxx, S.x.y = Sxy, S.y.y = Syy, S.y.x = Syx;
      
      foreach(){
	coord nrf;
	double mag = 0.;
	double wx = 0, wy = 0;
	foreach_dimension(){
	  nrf.x = sd[1] - sd[-1];
	  mag += sq(sd[1] - sd[-1]);
	}

	if(mag > 0.){
	  foreach_dimension()
	    nrf.x /= sqrt(mag);

	  if(sq(nrf.x) < thrsh){
	    wx = 0.;
	    wy = 1.;
	  }
	  else if(sq(nrf.x) > (1. - thrsh)){
	    wx = 1.;
	    wy = 0.;
	  }
	  else{
	    wx = sq(nrf.x);
	    wy = sq(nrf.y);
	  }
	}
	
	foreach_dimension() {
	  S.y.y[] = 0.;
	  if(interfacial(point,d)){
	    double syyhx = getdiagterm_y(point,dx,sigma,kappa[]);
	    double syyhy = getdiagterm_y(point,dy,sigma,kappa[]);
	    double syy = wx*syyhx + wy*syyhy;
	    S.y.y[] += syy;
	  }
	}
      }
      
      /**
      We compute the off-diagonal components of the tensor.  */
      
      foreach_vertex(){
	coord nrf;
	double mag = 0.;
	double wx = 0, wy = 0;
	
	foreach_dimension(){
	  nrf.x = sd[1] + sd[0,-1] - sd[-1] - sd[-1,-1];
	  mag += sq(sd[1] + sd[0,-1] - sd[-1] - sd[-1,-1]);
	}
	if(mag > 0.){
	  foreach_dimension()
	    nrf.x /= sqrt(mag);
	  
	  if(sq(nrf.x) < thrsh){
	    wx = 0.;
	    wy = 1.;
	  }
	  else if(sq(nrf.x) > (1. - thrsh)){
	    wx = 1.;
	    wy = 0.;
	  }
	  else{
	    wx = sq(nrf.x);
	    wy = sq(nrf.y);
	  }
	}
	foreach_dimension(){
	  if(interfacialvertex(point,d)){
	    double sxyhx = getnondiagterm_y(point,dx,sigma);
	    double sxyhy = getnondiagterm_y(point,dy,sigma);
	    S.x.y[] = wx*sxyhx + wy*sxyhy;
	  }
	  else
	    S.x.y[] = 0.;
	}
      }
      
      face vector av = a;
      foreach_face() {
	av.x[] += alpha.x[]/(fm.x[] + SEPS)*(S.x.x[] - S.x.x[-1] + S.x.y[0,1] - S.x.y[])/Delta;

#if AXI
	coord n = {
	  (nrm.x[] + nrm.x[-1])/2.,
	  (nrm.y[] + nrm.y[-1])/2.
	};	
	extern scalar f;
	av.x[] -= alpha.x[]/sq(fm.x[] + SEPS)*(sigma[] + sigma[-1])/2.*n.y*(f[] - f[-1])/Delta;

#endif // AXI

      }
    }
}