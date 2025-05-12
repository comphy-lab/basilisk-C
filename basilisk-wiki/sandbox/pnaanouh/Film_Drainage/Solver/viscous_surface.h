/**
# Surface stress continuity and viscosity

The complete surface stress continuity of a free surface (without surface
 tension) writes
 $$
(\mathbf{T_{liquid}} - \mathbf{T_{gas}}) \cdot \mathbf{n} = 0
$$

Projecting this equation along normal and tangential directions gives (in 2D)
$$
\begin{aligned}
  &\phi|_{top} = - \frac{2\,\nu}{1 + \eta_{x}^2} (\partial_x u|_{top} 
      (1 - \eta_{x}^2) - \eta_{x} (\partial_z u + \partial_x w))_{top} \\
  &\frac{\rho\,\nu}{1 + \eta_{x}^2}(1-\eta_{x}^2) (\partial_z u + \partial_x w)_{top} 
    = 4\,\eta_{x}\,\partial_x u|_{top}
\end{aligned}
$$
where &\phi_{nu}$ is the non-hydrostratique pressure. 

In the limits of small slopes, we get the simplified equations
$$
\begin{aligned}
  &\phi|_{top} = - 2\,\nu\,\partial_x u|_{top} \\
  &\partial_z u|_{top} = - \partial_x w|_{top} 
\end{aligned}
$$

These terms are added respecively in the momentum equation and in the
 vertical viscosity diffusion. */

bool visc_activate = true;
 
/**
## Normal stress continuity
The pressure contribution obtained in the normal stress continuty can be added
 to the momentum equation of the [multilayer solver](hydro.h) as
$$
\begin{aligned}
  \partial_t \left( h \mathbf{u} \right)_k + \mathbf{{\nabla}} \cdot \left(
  h \mathbf{u}  \mathbf{u} \right)_k & = - gh_k  \mathbf{{\nabla}} (\eta)
  {\color{blue} - \mathbf{{\nabla}} (h \phi_{\nu})_k + \left[ \phi_{\nu} 
  \mathbf{{\nabla}} z \right]_k}
\end{aligned}
$$
where the terms in blue have been added and $\phi_{\nu}$ is the viscous
pressure deviation due to velocity variations on the surface.

These terms are added to the acceleration of the [multilayer solver](hydro.h).*/


event acceleration (i++)
{
  /**
  The hydrostatic pressure deviation $\phi_{\nu}$ is stored on the interfaces
  between layers (consistently with the Keller box scheme
  discretisation, see [Popinet,
  2020](/Bibliography#popinet2020)). This gives the following vertical
  discrete integration scheme. */
  if (visc_activate){
    scalar phiNu = new scalar[nl];
    foreach(){
      double phiNu0 = 0.;
      foreach_dimension(){
	double etax=(eta[1]-eta[-1])/(2.*Delta);
        phiNu0 -= nu*2.*(1.+sq(etax))/(1.-sq(etax))*(u.x[1,0,nl-1] - u.x[-1,0,nl-1]+h[1,0,nl-1]/2*dut.x[1,0]-h[-1,0,nl-1]/2*dut.x[-1,0])/(2.*Delta);
	// fixme dx_s and dx are not the same 
      }
      foreach_layer()
        phiNu[] = phiNu0;
    }
    boundary ({phiNu});

    /**
    Once the pressure deviation is known, the terms in blue above are
    added to the face acceleration field `ha`, using the
    [pressure-gradient macro](hydro.h#horizontal-pressure-gradient). */
    
    foreach_face() {
      hpg_2(pg, phiNu, 0, ha.x[] += pg);
    }
    boundary ((scalar *){ha});
  
    delete ({phiNu});
  }
}


/**
## Tangential stress continuity 
The viscous contribution is added as a boundary condition in the
 vertical viscosity solver, an auxilliary field $du_{\nu}$ is needed to
 contain this condition.
*/

vector du_nu[];
#if NH

event viscous_term (i++) 
{
  if (visc_activate){
    double maxdiff=1.;
    int iter=0;

    scalar wt[];
    scalar eta_star[];
    vector du_nu0[];
    
    foreach()
      foreach_dimension (){
        du_nu.x[]=0;
	du_nu0.x[]=0;
    }
    
        foreach(){
      double b=zb[];
      foreach_layer(){
	  zl[]=b;
	  b+=h[];
      }
    }

    foreach(){
      eta_star[] = 0; 
      foreach_layer(){
	eta_star[]+=h[];
      }
    }
    eta_star[left]=eta[left];
    eta_star[right]=eta[right];
    
    
    while(iter<110){
      maxdiff=0.;

      foreach(){
	wt[]= 0.;
	if(nl>1){
	  wt[] += (h[0,0,nl-2]*(3.*h[0,0,nl-1]+2.*h[0,0,nl-2])*w[0,0,nl-1]+sq(h[0,0,nl-1])*w[0,0,nl-2])/sq(h[0,0,nl-1]+h[0,0,nl-2]);
	}
	else{
	  wt[] += 0;
	}
	foreach_dimension(){
	  wt[] -= (h[1,0,nl-1]*u.x[1,0,nl-1]-h[-1,0,nl-1]*u.x[-1,0,nl-1])/(2.*Delta);
	  if(nl>1){
	    wt[] += (((3.*sq(h[0,0,nl-1])+3.*h[0,0,nl-1]*h[0,0,nl-2]+sq(h[0,0,nl-2]))*u.x[0,0,nl-1]-sq(h[0,0,nl-1])*u.x[0,0,nl-2])/((h[0,0,nl-1]+h[0,0,nl-2])*(2.*h[0,0,nl-1]+h[0,0,nl-2]))+(1./2.)*du_nu0.x[0,0]*h[0,0,nl-1]*(h[0,0,nl-1]+h[0,0,nl-2])/(2.*h[0,0,nl-1]+h[0,0,nl-2]))*(zl[1,0,nl-1]+h[1,0,nl-1]-zl[-1,0,nl-1]-h[-1,0,nl-1])/(2.*Delta);
	    wt[] -= ((h[0,0,nl-2]*(3.*h[0,0,nl-1]+h[0,0,nl-2])*u.x[0,0,nl-1]+2.*sq(h[0,0,nl-1])*u.x[0,0,nl-2])/((h[0,0,nl-1]+h[0,0,nl-2])*(2.*h[0,0,nl-1]+h[0,0,nl-2]))-(1./2.)*du_nu0.x[0,0]*h[0,0,nl-1]*h[0,0,nl-2]/(2.*h[0,0,nl-1]+h[0,0,nl-2]))*(zl[1,0,nl-1]-zl[-1,0,nl-1])/(2.*Delta);
	  }
	  else{
	    if (symmetric_bathymetry){
	      wt[] += (u.x[0,0,nl-1]+h[0,0,nl-1]/6.*(2.*du_nu0.x[0,0]+dub.x[0,0]))*(zl[1,0,nl-1]+h[1,0,nl-1]-zl[-1,0,nl-1]-h[-1,0,nl-1])/(2.*Delta);
	      wt[] -= (u.x[0,0,nl-1]-h[0,0,nl-1]/6.*(du_nu0.x[0,0]+2.*dub.x[0,0]))*(zl[1,0,nl-1]-zl[-1,0,nl-1])/(2.*Delta);
	    }
	    else{
	      wt[] += ((3.*u.x[0,0,nl-1]-u_b.x[0,0])/2.+h[0,0,nl-1]*du_nu0.x[0,0]/4.)*(zl[1,0,nl-1]+h[1,0,nl-1]-zl[-1,0,nl-1]-h[-1,0,nl-1])/(2.*Delta);
	      wt[] -= u_b.x[0,0,nl-1]*(zl[1,0,nl-1]-zl[-1,0,nl-1])/(2.*Delta);
	    }
	  }
	}
	if(nl>1){
	  wt[] /= (2.*h[0,0,nl-2]+h[0,0,nl-1])/(h[0,0,nl-1]+h[0,0,nl-2]);
	}
      }
      
      foreach()
	foreach_dimension ()
	{
	  double etax=(eta_star[1]-eta_star[-1])/(2.*Delta);
	  du_nu.x[] = (- (wt[1,0,nl-1] - wt[-1,0,nl-1])/(2.*Delta) + 4.*(u.x[1,0,nl-1]+ h[1,0,nl-1]/2*du_nu0.x[1,0] - u.x[-1,0,nl-1]- h[-1,0,nl-1]/2*du_nu0.x[-1,0])/(2.*Delta)*etax/(1.-etax*etax))/(1.+etax*(4.*etax/(1-sq(etax))-etax));
	}
      iter++;
      foreach(reduction(max:maxdiff)){
	foreach_dimension(){
	  if(fabs((du_nu.x[]-du_nu0.x[])/(du_nu0.x[]+1E-13))>maxdiff){
	    maxdiff=fabs((du_nu.x[]-du_nu0.x[])/(du_nu0.x[]+1E-13));
	  }
	  du_nu0.x[]=du_nu.x[];
	  
	}
      }
      if (maxdiff<1E-5){
	iter=1000;
      }
    }
    boundary ((scalar *){du_nu});
    dut = du_nu;
  }
}

#else

event viscous_term (i++) 
{
  if (visc_activate){
    double etax=0;
    double maxdiff=1.;
    int iter=0;

    scalar wt[];
    scalar eta_star[];
    vector du_nu0[];

    scalar QED[];
    vertical_velocity(QED,hu,hf);
    
    foreach()
      foreach_dimension (){
      du_nu.x[]=0;
      du_nu0.x[]=0;
    }
    
    foreach(){
      double b=zb[];
      foreach_layer(){
	zl[]=b;
	b+=h[];
      }
    }

    foreach(){
      eta_star[] = 0; 
      foreach_layer(){
	eta_star[]+=h[];
      }
    }
    eta_star[left]=eta[left];
    eta_star[right]=eta[right];
    
    
    while(iter<110){
      maxdiff=0.;

      foreach(){
	wt[]= 0.;
	if(nl>1){
	  wt[] += (h[0,0,nl-2]*(3.*h[0,0,nl-1]+2.*h[0,0,nl-2])*QED[0,0,nl-1]+sq(h[0,0,nl-1])*QED[0,0,nl-2])/sq(h[0,0,nl-1]+h[0,0,nl-2]);
	}
	else{
	  wt[] += 0;
	}
	foreach_dimension(){
	  wt[] -= (h[1,0,nl-1]*u.x[1,0,nl-1]-h[-1,0,nl-1]*u.x[-1,0,nl-1])/(2.*Delta);
	  if(nl>1){
	    wt[] += (((3.*sq(h[0,0,nl-1])+3.*h[0,0,nl-1]*h[0,0,nl-2]+sq(h[0,0,nl-2]))*u.x[0,0,nl-1]-sq(h[0,0,nl-1])*u.x[0,0,nl-2])/((h[0,0,nl-1]+h[0,0,nl-2])*(2.*h[0,0,nl-1]+h[0,0,nl-2]))+(1./2.)*du_nu0.x[0,0]*h[0,0,nl-1]*(h[0,0,nl-1]+h[0,0,nl-2])/(2.*h[0,0,nl-1]+h[0,0,nl-2]))*(zl[1,0,nl-1]+h[1,0,nl-1]-zl[-1,0,nl-1]-h[-1,0,nl-1])/(2.*Delta);
	    wt[] -= ((h[0,0,nl-2]*(3.*h[0,0,nl-1]+h[0,0,nl-2])*u.x[0,0,nl-1]+2.*sq(h[0,0,nl-1])*u.x[0,0,nl-2])/((h[0,0,nl-1]+h[0,0,nl-2])*(2.*h[0,0,nl-1]+h[0,0,nl-2]))-(1./2.)*du_nu0.x[0,0]*h[0,0,nl-1]*h[0,0,nl-2]/(2.*h[0,0,nl-1]+h[0,0,nl-2]))*(zl[1,0,nl-1]-zl[-1,0,nl-1])/(2.*Delta);
	  }
	  else{
	    if (symmetric_bathymetry){
	      wt[] += (u.x[0,0,nl-1]+h[0,0,nl-1]/6.*(2.*du_nu0.x[0,0]+dub.x[0,0]))*(zl[1,0,nl-1]+h[1,0,nl-1]-zl[-1,0,nl-1]-h[-1,0,nl-1])/(2.*Delta);
	      wt[] -= (u.x[0,0,nl-1]-h[0,0,nl-1]/6.*(du_nu0.x[0,0]+2.*dub.x[0,0]))*(zl[1,0,nl-1]-zl[-1,0,nl-1])/(2.*Delta);
	    }
	    else{
	      wt[] += ((3.*u.x[0,0,nl-1]-u_b.x[0,0])/2.+h[0,0,nl-1]*du_nu0.x[0,0]/4.)*(zl[1,0,nl-1]+h[1,0,nl-1]-zl[-1,0,nl-1]-h[-1,0,nl-1])/(2.*Delta);
	      wt[] -= u_b.x[0,0,nl-1]*(zl[1,0,nl-1]-zl[-1,0,nl-1])/(2.*Delta);
	    }
	  }
	}
	if(nl>1){
	  wt[] /= (2.*h[0,0,nl-2]+h[0,0,nl-1])/(h[0,0,nl-1]+h[0,0,nl-2]);
	}
      }
      
      foreach()
	foreach_dimension ()
	{
	  etax=(eta_star[1]-eta_star[-1])/(2.*Delta);
	  du_nu.x[] = (- (wt[1,0,nl-1] - wt[-1,0,nl-1])/(2.*Delta) + 4.*(u.x[1,0,nl-1]+ h[1,0,nl-1]/2*du_nu0.x[1,0] - u.x[-1,0,nl-1]- h[-1,0,nl-1]/2*du_nu0.x[-1,0])/(2.*Delta)*etax/(1.-etax*etax))/(1.+etax*(4.*etax/(1-sq(etax))-etax));
	}
      iter++;
      foreach(){
	foreach_dimension(){
	  if(fabs((du_nu.x[]-du_nu0.x[])/(du_nu0.x[]+1E-13))>maxdiff){
	    maxdiff=fabs((du_nu.x[]-du_nu0.x[])/(du_nu0.x[]+1E-13));
	  }
	  du_nu0.x[]=du_nu.x[];
	  
	}
      }
      if (maxdiff<1E-5){
	iter=1000;
      }
    }
    boundary ((scalar *){du_nu});
    dut = du_nu;
  }
}

#endif
