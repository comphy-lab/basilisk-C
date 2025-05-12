/**
   Implementation of the Spalart-Allmaras model */

/**
The equation for the model is 
$$
\nu_t=\hat{\nu}f_{v1}
$$
with $f_{v1}=\frac{\chi^3}{\chi^3+c_{v1}^3}$ where $\chi=\frac{\hat{\nu}}{\nu}$ and $\hat{\nu}$ verify the equation (using Einstein notation)
$$
\frac{\partial\hat{\nu}}{\partial t}+u_j\frac{\partial\hat{\nu}}{\partial x_j}=
c_{b1}\left(1-f_{t2}\right)\hat{S}\hat{\nu}-\left(c_{w1}f_w-\frac{c_{b1}}{\kappa^2}f_{t2}\right) \left(\frac{\hat{\nu}}{d}\right)^2+
\frac{1}{\sigma}\left(\frac{\partial}{\partial x_j}\left(\left(\nu+\hat{\nu}\right)\frac{\partial\hat{\nu}}{\partial x_j}\right)+c_{b2}\frac{\hat{nu}}{\partial x_i}\frac{\partial\hat{\nu}}{\partial x_i}\right)
$$
where
$$
\hat{S}=\Omega+\frac{\hat{\nu}}{\kappa^2 d^2}f_{v2} ~~~~~   \Omega=\sqrt{2W_{ij}W_{ij}} ~~~~~~  \text{with} ~~~\bold{W}=\frac{1}{2}\left(\nabla \bold{u}-\nabla\bold{u}^t\right) ~~~~ \text{the rotation tensor}
$$
d is the distance to the nearest wall
$$
f_{v2}=1-\frac{\chi}{1+\chi f_{v1}}
$$
$$
f_{w}=g\left(\frac{1+c_{w3}^6}{g^6+c_{w3}^6}\right)^{\frac{1}{6}}  ~~~~ g=r+c_{w2}\left(r^6-r\right) ~~~  \text{with} ~~~~ r=min\left(\frac{\hat{\nu}}{\hat{S}\kappa^2 d^2},10\right)
$$
and
$$
f_{t2}=c_{t3}e^{-c_{t4}\chi^2}
$$
The boundary condition is $\hat{\nu}_{Wall}=0$

The value for the constants is
$$
c_{b1}=0.1355,~ \sigma=\frac{2}{3},~ c_{b2}=0.622,~ \kappa=0.41,~ c_{w2}=0.3,~ c_{w3}=2,~ c_{v1}=7.1, ~c_{t3}=1.2,~ c_{t4}=0.5 ~\text{and} ~c_{w1}=\frac{c_{b1}}{\kappa^2}+\frac{1+c_{b2}}{\sigma}
$$

For stability reason, we use
$$
\hat{S}=\Omega+\bar{S}~~~ \text{when}~~~\bar{S}\geq-c_2\Omega
$$ 
$$
\hat{S}=\Omega+\frac{\Omega\left(c_2^2\Omega+c_3\bar{S}\right)}{\left(c_3-2c_2\right)\Omega-\bar{S}} ~~~\text{when}~~~\bar{S}<-c_2\Omega
$$
with $c_2=0.7$ and $c_3=0.9$

*/


scalar muhat;   //field that the Spalart-Allmaras equation propagated
face vector mut;    //turbulent viscosity

scalar source;   //source term


/**
   Definition of all the constant for the model */
double cb1,sigma,cb2,kappa,cw2,cw3,cv1,ct3,ct4,cw1,c2,c3;
double molvis;    //molecular viscosity


event defaults(i=0) {
  cb1=0.1355;
  sigma=2./3.;
  cb2=0.622;
  kappa=0.41;
  cw2=0.3;
  cw3=2.;
  cv1=7.1;
  ct3=1.2;
  ct4=0.5;
  cw1=cb1/(kappa*kappa)+(1+cb2)/sigma;
  c2=0.7;
  c3=0.9;
  molvis=1.;
  muhat=new scalar[nl];
  mut=new face vector[nl];
  source=new scalar[nl];
  foreach() {
    foreach_layer() {
      dimensional (source[]=Delta*molvis/t);
    }
  }
  nu=mut;

}

/**
   Definition of the different functions used by the model */

double fv1(double chi) {
  return pow(chi,3.)/(pow(chi,3.)+pow(cv1,3.));
}

double fv2(double chi) {
  return 1-chi/(1+chi*fv1(chi));
}

double fw(double r) {
  double g=r+cw2*(pow(r,6.)-r);
  return g*pow((1+pow(cw3,6.))/(pow(g,6.)+pow(cw3,6.)),1./6.);
}

double ft2(double chi) {
  return ct3*exp(-ct4*pow(chi,2.));
}

/**
   We first advect the viscosity */

void advection_source() {

  foreach_face() {
    foreach_layer() {
      hf.x[]=(h[]+h[-1])/2.;
      hu.x[]=(h[]*u.x[]+h[-1]*u.x[-1])/2.;
    }
  }
  advect_tracer((scalar *) {muhat},hu,hf,dt,(scalar *) {source});

}

/**
   We take care of the diffusion with a semi-implicit method */


/**
   We need to add a correction to the viscosity computed to account for the term in the right hand side */

void correction_nu(double dt) {
  foreach() {
    foreach_layer() {
      muhat[]+=dt*source[]/h[];
    }
  }
}

event face_fields(i++) {
  foreach() {
    foreach_layer() {
      double dt=sq(Delta)/(muhat[]+1e-10);
      if (10*dt<dtmax) {
	dtmax=10*dt;
      }
    }
  }
}

event viscous_term(i++) {
  horizontal_diffusion({u,w},nu,dt);
}


event reaction_diffusion(i++,last) {
  
  advection_source();

  //correction_nu(dt);

  face vector muhat2;
  face vector muhat3;
  muhat2=new face vector[nl];
  muhat3=new face vector[nl];
  foreach_face() {
    foreach_layer() {
      muhat2.x[]=1/sigma*(molvis+(muhat[]+muhat[-1])/2.);
      muhat3.x[]=1/sigma*(molvis+muhat[]);
    }
  }
  
  
  foreach() {
    vertical_diffusion(point,x,h,muhat,dt,muhat3,0,0.*molvis,0,true,false,3.*molvis,molvis/sigma);
  }

  
  foreach_face() {
    foreach_layer() {
      muhat2.x[]=1/sigma*(molvis+(muhat[]+muhat[-1])/2.);
      muhat3.x[]=1/sigma*(molvis+muhat[]/2.);
    }
  }
  
  horizontal_diffusion((scalar *) {muhat},muhat2,dt);
  
  delete ((scalar *) {muhat2,muhat3});
  
  //correction_nu(-dt);

  
  trash({source});
  
  double chi;
  
  double Omega2;                 //norm of the rotation tensor
  double S=0.;
  double d=0.;              //distance to the nearest wall
  
    
  if (i>=100) {
    foreach() {
      double xpos=x;
      double zlayer=zb[];
      foreach_layer() {
	double sum=0.;
	foreach_dimension() {
	  if (point.l==0) {
	    if (xpos<=1.) {
	      sum+=((u.x[0,0,1]-u.x[0,0,0])/(3./2.*h[]+h[0,0,1]/2.)-(w[1]-w[-1])/(2.*Delta))*((u.x[0,0,1]-u.x[0,0,0])/(3./2.*h[]+h[0,0,1]/2.)-(w[1]-w[-1])/(2.*Delta));
	    }
	    else {
	      sum+=((u.x[0,0,1]+u.x[0,0,0])/(3./2.*h[]+h[0,0,1]/2.)-(w[1]-w[-1])/(2.*Delta))*((u.x[0,0,1]+u.x[0,0,0])/(3./2.*h[]+h[0,0,1]/2.)-(w[1]-w[-1])/(2.*Delta));
	    }
	  }
	  else if (point.l==nl-1) {
	   sum+=((u.x[0,0,0]-u.x[0,0,-1])/(3./2.*h[]+h[0,0,-1]/2.)-(w[1]-w[-1])/(2.*Delta))*((u.x[0,0,0]-u.x[0,0,-1])/(3./2.*h[]+h[0,0,-1]/2.)-(w[1]-w[-1])/(2.*Delta));
	  }
	  else {
	    sum+=((u.x[0,0,1]-u.x[0,0,-1])/(h[0,0,-1]/2.+h[0,0,0]+h[0,0,1]/2.)-(w[1]-w[-1])/(2.*Delta))*((u.x[0,0,1]-u.x[0,0,-1])/(h[0,0,-1]/2.+h[0,0,0]+h[0,0,1]/2.)-(w[1]-w[-1])/(2.*Delta));
	  }
	  //new version
	  //if (point.l==0) {
	  //  if (xpos<=1.) {
	  //    sum+=((u.x[0,0,1]-u.x[])/(h[]+h[0,0,1])-(w[1]-w[-1])/(2.*Delta))*((u.x[0,0,1]-u.x[])/(h[]+h[0,0,1])-(w[1]-w[-1])/(2.*Delta));
	  //  }
	  //  else {
	  //    sum+=((u.x[0,0,1]-u.x[])/(h[]+h[0,0,1])+u.x[]/h[]-(w[1]-w[-1])/(2.*Delta))*((u.x[0,0,1]-u.x[])/(h[]+h[0,0,1])+u.x[]/h[]-(w[1]-w[-1])/(2.*Delta));
	  //  }
	  //}
	  //else if (point.l==nl-1) {
	  //  sum+=((u.x[]-u.x[0,0,-1])/(h[]+h[0,0,-1])-(w[1]-w[-1])/(2.*Delta))*((u.x[]-u.x[0,0,-1])/(h[]+h[0,0,-1])-(w[1]-w[-1])/(2.*Delta));
	  //}
	  //else {
	  //  sum+=((u.x[0,0,1]-u.x[])/(h[0,0,1]+h[])+(u.x[]-u.x[0,0,-1])/(h[]+h[0,0,-1])-(w[1]-w[-1])/(2.*Delta))*((u.x[0,0,1]-u.x[])/(h[0,0,1]+h[])+(u.x[]-u.x[0,0,-1])/(h[]+h[0,0,-1])-(w[1]-w[-1])/(2.*Delta));
	  //}
	  //#if dimension==3
	  //sum+=(uf.x[0,0,1]-uf.x[])*(uf.x[0,0,1]-uf.x[]-uf.z[1,0,0]+uf.z[])/Delta/Delta;
	  //#endif      
	}
	Omega2=max(sum,0.);
	chi=muhat[]/molvis;
#if WALL
	d=distance_to_wall(x,zlayer+h[]/2.);
	zlayer+=h[];
#else
	d=HUGE;
#endif
	//S=max(sqrt(Omega2[])+muhat[]/(pow(kappa*d,2.)+1e-10)*fv2(chi),0.3*sqrt(Omega2[]));
	S=muhat[]/(pow(kappa*d,2.)+1e-10)*fv2(chi);
	if (S>=-c2*sqrt(Omega2)) {
	  S+=sqrt(Omega2);
	}
	else {
	  S=sqrt(Omega2)+sqrt(Omega2)*(pow(c2,2.)*sqrt(Omega2)+c3*S)/((c3-2*c2)*sqrt(Omega2)-S);
	}
	source[]=h[]*cb1*(1-ft2(chi))*S*muhat[];
	
     
	double r=min(muhat[]/(S*pow(kappa*d,2.)+1e-10),10);
	source[]+=-h[]*(cw1*fw(r)-cb1/(pow(kappa,2.))*ft2(chi))*pow(muhat[]/(d+1e-10),2.);
	
	foreach_dimension() {
	  source[]+=h[]*cb2/sigma*sq((muhat[1]-muhat[-1])/(2*Delta));
	}
	
	if (point.l==0) {
	  if (x<=1.) {
	    source[]+=h[]*cb2/sigma*sq((muhat[0,0,1]-muhat[0,0,0])/(3./2.*h[]+h[0,0,1]/2.));
	  }
	  else {
	    source[]+=h[]*cb2/sigma*sq((muhat[0,0,1]+muhat[0,0,0])/(3./2.*h[]+h[0,0,1]/2.));
	  }
	}
	else if (point.l<nl-1) {
	  source[]+=h[]*cb2/sigma*sq((muhat[0,0,1]-muhat[0,0,-1])/(h[0,0,-1]/2.+h[]+h[0,0,1]/2.+dry));
	}
	else if (point.l==nl-1) {
	  source[]+=h[]*cb2/sigma*sq((muhat[0,0,0]-muhat[0,0,-1])/(3./2.*h[]+h[0,0,-1]/2.+dry));
	}
	  //new version
	//if (point.l==0) {
	//  if (xpos<=1.) {
	//    source[]+=h[]*cb2/sigma*sq((muhat[0,0,1]-muhat[])/(h[]+h[0,0,1]));
	//  }
	//  else {
	//    source[]+=h[]*cb2/sigma*sq((muhat[0,0,1]-muhat[])/(h[]+h[0,0,1])+muhat[]/h[]);
	//  }
	//}
	//else if (point.l<nl-1) {
	//  source[]+=h[]*cb2/sigma*sq((muhat[0,0,1]-muhat[])/(h[]+h[0,0,1]+dry)+(muhat[]-muhat[0,0,-1])/(h[]+h[0,0,-1]));
	//}
	//else if (point.l==nl-1) {
	//  source[]+=h[]*cb2/sigma*sq((muhat[]-muhat[0,0,-1])/(h[]+h[0,0,-1]+dry));
	//}
      }
    }
  }
  
  
  correction_nu(dt);

  
  foreach() {
    foreach_layer() {
      muhat[]=max(muhat[],0.);
    }
  }
  
  foreach_face() {
    foreach_layer() {
      chi=(muhat[]+muhat[-1])/2./molvis;
      mut.x[]=fm.x[]*(molvis+fv1(chi)*molvis*chi);
    }
  }
  
  delete ((scalar *) {hu});
  delete ((scalar *) {hf});
}

event cleanup (t = end, last)
{
  delete ((scalar *) {mut});
  delete ((scalar *) {source,muhat});
}