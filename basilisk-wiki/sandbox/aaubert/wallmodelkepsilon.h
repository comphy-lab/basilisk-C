/**
   Wall function used for the k-epsilon model */

#include "functionkepsilon.h"

attribute {
  scalar kwall;
  scalar epsilonwall;
  scalar utau;
  scalar utIP;
  scalar duF;
}

event defaults(i=0) {
  d_IP=2.*Deltamin;//sqrt(2.)*Deltamin;

  scalar kwall=new scalar;
  mutke.kwall=kwall;
  scalar epsilonwall=new scalar;
  mutke.epsilonwall=epsilonwall;
  scalar utau=new scalar;
  mutke.utau=utau;
  scalar utIP=new scalar;
  mutke.utIP=utIP;
  scalar duF=new scalar;
  mutke.duF=duF;

}
#if 1
/**
   Condition for the velocity on the solid */
double wall_condition_velocity(double x2, double y2,bool xny, scalar s) {

  
  bool reverse=false;
  double utau=0.;
  double d=0.;
  foreach_point(x2,y2) {
    d=distance_to_wall(x,y);
    if (fabs(d)<Delta) {
      utau=s.utau[];
      if (utau<0) {  //we take the information if we need to reverse the tangential vector
	      reverse=true;
	      utau=-utau;
      }
    }
  }
  if (utau==0.) {
    return 0.;
  }

  double ut;

#if LINEARISATION  //linearise version
  foreach_point(x2,y2) {
    ut=s.utIP[]-d_IP*s.duF[];
  }
#else //not linearise
  ut=0.;
#endif

  coord n,t;
  double dx,dy;
  foreach_point(x2,y2) {
    dx=(distance_to_wall(x+Delta,y)-distance_to_wall(x-Delta,y))/(2.*Delta);
    dy=(distance_to_wall(x,y+Delta)-distance_to_wall(x,y-Delta))/(2.*Delta);
	  
    n.x=dx/sqrt(sq(dx)+sq(dy)+1e-10);
    n.y=dy/sqrt(sq(dx)+sq(dy)+1e-10);

    t.x=n.y;
    t.y=-n.x;
    if (reverse) {
      t.x=-t.x;
      t.y=-t.y;
    }
  }

  return xny ? t.x*ut : t.y*ut;
}


double wall_condition_k(double x2, double y2, scalar s) {
  double kvalue=0.;
  foreach_point(x2,y2) {
    kvalue=s.kwall[];
  }
  return kvalue;
}

double wall_condition_epsilon(double x2, double y2, scalar s) {
  double epsilonvalue=0.;
  foreach_point(x2,y2) {
    epsilonvalue=s.epsilonwall[];
  }
  return epsilonvalue;
}

#endif


/**
   We applied the wall model in this event */
event wall_model(i++) {

  scalar ux[];
  scalar uy[];
  scalar pvalue[];
  scalar kvalue[];
  scalar epsilonvalue[];

  if (i>100) {

    foreach() {
      double dstar=distance_to_wall(x,y);
      if ((dstar>-Deltamin)&&(dstar<d_IP)&&(x>leading)&&(x<trailing)) {
      
        double nx,ny;
#if EMBED        
        if ((cs[]>0.)&&(cs[]<1.)) {
          coord n,b;
          embed_geometry(point,&b,&n);
          nx=-n.x;
          ny=-n.y;
        }
        else {
#endif          
          double dx=(distance_to_wall(x+Delta,y)-distance_to_wall(x-Delta,y))/(2.*Delta);
          double dy=(distance_to_wall(x,y+Delta)-distance_to_wall(x,y-Delta))/(2.*Delta);
          nx=dx/sqrt(sq(dx)+sq(dy)+1e-20);
          ny=dy/sqrt(sq(dx)+sq(dy)+1e-20);
#if EMBED          
        }
#endif        

        double tx=ny;
        double ty=-nx;

        double alpha=d_IP-dstar;
        double x_IP=x+alpha*nx;
        double y_IP=y+alpha*ny;

        double ux_IP=interpolate(u.x,x_IP,y_IP);
        double uy_IP=interpolate(u.y,x_IP,y_IP);
        double p_IP=interpolate(p,x_IP,y_IP);
        double ut_IP=ux_IP*tx+uy_IP*ty;
        bool reverse=false;
        if (ut_IP<0.) {
	        reverse=true;
	        tx=-tx;
	        ty=-ty;
	        ut_IP=-ut_IP;
        }
        double un_IP=ux_IP*nx+uy_IP*ny;

        double utau=newton_utau(0.1,ut_IP,d_IP,molvis[]);

      
        kvalue[]=max(sq(utau)/sqrt(C_mu),pow(10.,-10.));
        epsilonvalue[]=rho[]*C_mu*sq(kvalue[])/(rho[]*molvis[]*viscosity_profile(utau,d_IP,ut_IP));

        mutke.kwall[]=kvalue[];
        mutke.epsilonwall[]=epsilonvalue[];

#if LINEARISATION
      
        mutke.utIP[]=ut_IP;
        mutke.duF[]=1./(rho[]*molvis[]+rho[]*C_mu*sq(kvalue[])/epsilonvalue[])*rho[]*sq(utau);
        mutke.utau[]=utau;

        double ut=ut_IP+mutke.duF[]*(dstar-d_IP);
        double un=-un_IP*sq(dstar/d_IP)+2.*un_IP*dstar/d_IP;

      
#else
        double ut=utau*fwallSA(dstar*utau/molvis[]);
        double un=un_IP*dstar/d_IP;
#endif
        ux[]=un*nx+ut*tx;
        uy[]=un*ny+ut*ty;
        pvalue[]=p_IP;
      }
    }

  

    foreach() {
      double dstar=distance_to_wall(x,y);
      if ((dstar>-Deltamin)&&(dstar<d_IP)&&(x>leading)&&(x<trailing)) {
        u.x[]=ux[];
        u.y[]=uy[];
        p[]=pvalue[];
        k[]=kvalue[];
        epsilon[]=epsilonvalue[];
      }
    }
    foreach() {
      mutke[]=rho[]*C_mu*sq(k[])/epsilon[];
    }

    foreach_face() {
      mut.x[]=fm.x[]*face_value(rho,0)*face_value(molvis,0)+fm.x[]*face_value(mutke,0);
    }
  }
}

/**
   Boundary condition need to be change. This event allows to hook up other boundary conditions that need to be modified (for example if a solid is a frontier of the domain) */
event boundary_condition(i++) {
#if EMBED
    k[embed]=dirichlet(wall_condition_k(x,y,mutke));
    epsilon[embed]=dirichlet(wall_condition_epsilon(x,y,mutke));
    u.n[embed]=dirichlet(wall_condition_velocity(x,y,true,mutke));
    u.t[embed]=dirichlet(wall_condition_velocity(x,y,false,mutke));
#endif
  boundary({k,epsilon});
  foreach() {
    mutke[]=rho[]*C_mu*sq(k[])/epsilon[];
  }
  foreach_face() {
    mut.x[]=fm.x[]*face_value(rho,0)*face_value(molvis,0)+fm.x[]*face_value(mutke,0);
  }
}