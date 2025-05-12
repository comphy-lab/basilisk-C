/**
   Wall function used for the Spalart-Allmaras model 
   
   We use a classical wall function given by Spalart and Allmaras in [here](https://www.iccfd.org/iccfd7/assets/pdf/papers/ICCFD7-1902_paper.pdf). This function is an approximation of a real solution of their model. We also try to implement a wall model proposed by [Berger and Aftosmis](https://www.nas.nasa.gov/publications/software/docs/cart3d/pages/publications/AIAA_2017-0528_berger_aftosmis.pdf) based on the solution of a boundary value problem to better take into account the pressure gradient near the solid. This wall model is currently not working*/

#include "functionSA.h"

#if !defined(BVP4)           //version with a resolution of the boundary layer equations
#define BVP4 0
#endif

#if !defined(LINEARISATION)      //use non linearize version for the wall function if not state otherwise
#define LINEARISATION 0
#endif

attribute {
  scalar utau;
  scalar duF;
  scalar utIP;
  scalar Cplus;
  scalar uboundary;
  scalar dnuhatF;
  scalar utauprev;
  scalar nuhatIP;
  double uxIPedge;
  double uyIPedge;
  double molvisIPedge;
  scalar uyboundary;
  double pIPedge;
  double pfIPedge;
  scalar unIPprevious;
  scalar dunIPdt;
  double dunIPdtedge;
  double unIPpreviousedge;
  scalar pIP;
  scalar pfIP;
  scalar viscSAdiss;
}

typedef struct{
  int npoint;
  int nitermax,nitermin;
  double niteravg;
  double norm2max, norm2avg;
}BVPstatsglobal;

BVPstatsglobal  statssolver;  //global information on the BVP solver

scalar pavg3[]; //average of the pressure

#if BVP4

#define EULER 0
#define RK4 0
#define MIRK221L 0
#define MIRK6 1

#if !defined(NSUBPOINT)
#define NSUBPOINT 51
#endif

#define NEQUATION 4
#define NBOUNDARY0 2

int Nsubpoint;  //number of point in the mesh
int Nequation;  //number of first order equation to solve
int Nboundary0; //number of boundary condition at $y=0$

scalar * smesh;
scalar * ssol;

double U0; //free stream speed

scalar updatedbvp[]; //mark the cells where the BVP have been solved

/**
   The BVP solver needs function to compute the right hand side, the Jacobian and the boundary conditions. This function are defined here */


typedef struct {
  double dp, adv, yF, uF, nuhatF, utau, rho, mu;
} Param;


//Function for blending of the advection term
double blending(double y, double y_IP, double utau, double nu) {
  return fwallSA(y*utau/nu)/fwallSA(y_IP*utau/nu);
}

//RHS of the differential equation
void rhs(double *res,double y, double *w, int niter, int Nitermax, Param param) {
  
  //Version originale
  double nu=param.mu/param.rho;
  res[0]=w[2]/(param.mu+param.rho*w[1]*fv1(w[1]/nu));
  res[1]=w[3]/(nu+w[1]);
  res[2]=param.dp+blending(y,param.yF,param.utau,nu)*param.rho*param.adv;
  double Omega=fabs(w[2]/(param.mu+param.rho*w[1]*fv1(w[1]/nu)));
  double Shat=w[1]*fv2(w[1]/nu)/sq(kappa*y+1e-20);
  if (w[1]>=0.) {
    if (Shat>-c2*Omega) {
      Shat+=Omega;
    }
    else {
      Shat=Omega+Omega*(sq(c2)*Omega+c3*Shat)/(Omega*(c3-2.*c2)-Shat+1e-30);
    }
  }
  else { //$\tilde{\nu}<0$
    Shat=(1.-ct3)*Omega;
  }
  double r,fwv;
  if (w[1]>=0.) {
    r=fabs(w[1]/(Shat*sq(kappa*y+1e-20)+1e-30))<10. ? min(w[1]/(Shat*sq(kappa*y+1e-20)+1e-30),10.) : 10.;
    fwv=fw(r);
  }
  else {  //$\tilde{\nu}<0$
    r=0.;
    fwv=-1;
  }
  Shat*=pow(2.,niter-Nitermax+1);
  res[3]=-cb2*sq(w[3]/(nu+w[1]))-sigma*(cb1*Shat*w[1]-cw1*fwv*sq(w[1]/(y+1e-20)));
}

//derivative of the RHS of the differential equation
void drhs(int Nequation, double dres[Nequation][Nequation], double y, double *w, int niter, int Nitermax, Param param) {
  
  // Version originale
  double nu=param.mu/param.rho;
  double fv1prime=3*sq(w[1])*pow(cv1*nu,3.)/sq(pow(w[1],3.)+pow(cv1*nu,3.));
  double fv2prime=-(nu-sq(w[1])*fv1prime)/sq(nu+w[1]*fv1(w[1]/nu));
  double S=w[1]*fv2(w[1]/nu)/sq(kappa*y+1e-20);
  double Omega=fabs(w[2]/(param.mu+param.rho*w[1]*fv1(w[1]/nu)));
  double dSd1=(fv2(w[1]/nu)+w[1]*fv2prime)/sq(kappa*y+1e-20);
  double dOmegad1=w[2]>0. ? -w[2]*(param.rho*fv1(w[1]/nu)+param.rho*w[1]*fv1prime)/sq(param.mu+param.rho*w[1]*fv1(w[1]/nu)) : w[2]*(param.rho*fv1(w[1]/nu)+param.rho*w[1]*fv1prime)/sq(param.mu+param.rho*w[1]*fv1(w[1]/nu));
  double dOmegad2=w[2]>0. ? 1./(param.mu+param.rho*w[1]*fv1(w[1]/nu)) : -1./(param.mu+param.rho*w[1]*fv1(w[1]/nu));
  double Shat, dShatd1, dShatd2;
  if (w[1]>=0.) {
    if (S>-c2*Omega) {
      Shat=S+Omega;
      dShatd1=dSd1+dOmegad1;
      dShatd2=dOmegad2;
    }
    else {
      Shat=Omega+Omega*(sq(c2)*Omega+c3*S)/(Omega*(c3-2.*c2)-S+1e-30);
      dShatd1=dOmegad1+(-c3*sq(S)*dOmegad1+Omega*(sq(c2)*dOmegad1+c3*dSd1)*Omega*(c3-2.*c2)+sq(c2)*sq(Omega)*dSd1-2.*sq(c2)*S*Omega*dOmegad1)/sq(Omega*(c3-2.*c2)-S+1e-30);
      dShatd2=dOmegad2+(-c3*sq(S)*dOmegad2+Omega*sq(c2)*dOmegad2*Omega*(c3-2.*c2)-2.*sq(c2)*S*Omega*dOmegad2)/sq(Omega*(c3-2.*c2)-S+1e-30);
    }
  }
  else { //$\tilde{\nu}<0$
    Shat=(1.-ct3)*Omega;
    dShatd1=(1.-ct3)*dOmegad1;
    dShatd2=(1.-ct3)*dOmegad2;
  }
  double r, drd1, drd2;
  double g, gprime;
  double fw, fwprime;
  if (w[1]>=0.) {
    if (fabs(w[1]/(Shat*sq(kappa*y+1e-20)+1e-30))<=10.) {
      r=w[1]/(Shat*sq(kappa*y+1e-20)+1e-30);
      drd1=(Shat-w[1]*dShatd1)/sq(Shat*kappa*y+1e-20);
      drd2=-w[1]*dShatd2/sq(Shat*kappa*y+1e-20);
    }
    else {
      r=10.;
      drd1=0.;
      drd2=0.;
    }
    g=r+cw2*(pow(r,6.)-r);
    gprime=1.+cw2*(6.*pow(r,5.)-1.);
    fw=g*pow((1+pow(cw3,6.))/(pow(g,6.)+pow(cw3,6.)),1./6.);
    fwprime=pow((1.+pow(cw3,6.))/(pow(g,6.)+pow(cw3,6.)),1./6.)-pow(g,6.)*(1.+pow(cw3,6.))/sq(pow(g,6.)+pow(cw3,6.))*pow((1.+pow(cw3,6.))/(pow(g,6.)+pow(cw3,6.)),-5./6.);
  }
  else {//$\tilde{\nu}<0$
    r=1;
    drd1=0.;
    drd2=0.;
    g=0.;
    gprime=0.;
    fw=-1;
    fwprime=0.;
  }
  Shat*=pow(2.,niter-Nitermax+1);
  dShatd1*=pow(2.,niter-Nitermax+1);
  dShatd2*=pow(2.,niter-Nitermax+1);
  
  dres[0][0]=0.;
  dres[0][1]=-w[2]*(param.rho*fv1(w[1]/nu)+param.rho*w[1]*fv1prime)/sq(param.mu+param.rho*w[1]*fv1(w[1]/nu));
  dres[0][2]=1./(param.mu+param.rho*w[1]*fv1(w[1]/nu));
  dres[0][3]=0.;
  dres[1][0]=0.;
  dres[1][1]=-w[3]/sq(nu+w[1]);
  dres[1][2]=0.;
  dres[1][3]=1./(nu+w[1]);
  dres[2][0]=0.;
  dres[2][1]=0.;
  dres[2][2]=0.;
  dres[2][3]=0.;
  dres[3][0]=0.;
  dres[3][1]=2.*cb2*sq(w[3])/pow(nu+w[1],3.)-sigma*(cb1*(Shat+w[1]*dShatd1)-cw1*(2.*w[1]/sq(y+1e-20)*fw+sq(w[1]/(y+1e-20))*fwprime*gprime*drd1));
  dres[3][2]=-sigma*(cb1*w[1]*dShatd2-cw1*sq(w[1]/(y+1e-20))*fwprime*gprime*drd2);
  dres[3][3]=-2.*cb2*w[3]/sq(nu+w[1]);
}

//boundary condition at 0
void boundary0(double *boundary, double *w, Param param) {
  boundary[0]=w[0];
  boundary[1]=w[1];
}

//boundary condition at F
void boundaryF(double *boundary, double *w, Param param) {
  boundary[0]=w[0]-param.uF;
  boundary[1]=w[1]-param.nuhatF;///param.mu*param.rho;
}

//derivative of the boundary condition at 0
void dboundary0(int Nequation,int Nboundary0, double dboundary[Nboundary0][Nequation], double *w, Param param) {
  dboundary[0][0]=1.;
  dboundary[0][1]=0.;
  dboundary[0][2]=0.;
  dboundary[0][3]=0.;
  dboundary[1][0]=0.;
  dboundary[1][1]=1.;
  dboundary[1][2]=0.;
  dboundary[1][3]=0.;
}

//derivative of the boundary condition at F
void dboundaryF(int Nequation, int Nboundary0, double dboundary[Nequation-Nboundary0][Nequation], double *w, Param param) {
  dboundary[0][0]=1.;
  dboundary[0][1]=0.;
  dboundary[0][2]=0.;
  dboundary[0][3]=0.;
  dboundary[1][0]=0.;
  dboundary[1][1]=1.;
  dboundary[1][2]=0.;
  dboundary[1][3]=0.;
}


#include "BVPsolver.h"

BVPstats statssolverpart = {0};

#endif



event defaults (i=0) {  //Deltamin need to be initialised in the main code
  
  d_IP=2.*Deltamin;  //distance of image point equals 2*Delta
  scalar duF=new scalar;
  muhat.duF=duF;
  scalar utau=new scalar;
  muhat.utau=utau;
  scalar utIP=new scalar;
  muhat.utIP=utIP;
  scalar Cplus=new scalar;
  muhat.Cplus=Cplus;
  scalar uboundary=new scalar;
  muhat.uboundary=uboundary;
  scalar dnuhatF=new scalar;
  muhat.dnuhatF=dnuhatF;
  scalar utauprev=new scalar;
  muhat.utauprev=utauprev;
  scalar nuhatIP=new scalar;
  muhat.nuhatIP=nuhatIP;
  scalar uyboundary=new scalar;
  muhat.uyboundary=uyboundary;
  scalar unIPprevious=new scalar;
  muhat.unIPprevious=unIPprevious;
  scalar dunIPdt=new scalar;
  muhat.dunIPdt=dunIPdt;
  scalar pIP=new scalar;
  muhat.pIP=pIP;
  scalar pfIP=new scalar;
  muhat.pfIP=pfIP;
  scalar viscSAdissinit=new scalar;
  muhat.viscSAdiss=viscSAdissinit;

  foreach() {
    muhat.duF[]=0.;
    muhat.utau[]=0.;
    muhat.Cplus[]=0.;
    muhat.uboundary[]=0.;
    muhat.uyboundary[]=0.;
    muhat.dnuhatF[]=0.;
    muhat.nuhatIP[]=0.;
    muhat.viscSAdiss[]=0.;
  }
  
#if BVP4

  Nsubpoint=NSUBPOINT;
  Nequation=NEQUATION;
  Nboundary0=NBOUNDARY0;
  statssolverpart=(BVPstats){0};

  U0=1.;  //free stream speed equals 1, can be overrided by the user

  ssol=malloc(1.*sizeof(scalar));
  for (int kk=0;kk<Nequation*Nsubpoint;kk++) {
    scalar s = new scalar;
    qrealloc (ssol, kk+2, scalar);  //fixme : same things as list_append but initial malloc gives a wrong size that needs to be compensate by the real length
    ssol[kk]=s;
    ssol[kk+1].i=-1;
    //ssol=list_append (ssol,s);
  }

  smesh=malloc(1.*sizeof(scalar));
  for (int kk=0;kk<Nsubpoint;kk++) {
    scalar s = new scalar;
    qrealloc (smesh, kk+2, scalar);  //fixme : same things as list_append but initial malloc gives a wrong size that needs to be compensate by the real length
    smesh[kk]=s;
    smesh[kk+1].i=-1;
    //ssol=list_append (ssol,s);
  }


#endif
}

#if BVP4
event init (i=0) {
  foreach() {
    updatedbvp[]=0.;
    double dstar=distance_to_wall(x,y);
    if ((dstar>-Deltamin)&&(dstar<d_IP)) {
      int kkm=0;
      for (scalar s in smesh) {
	s[]=kkm*d_IP/(1.*Nsubpoint-1.);  //initialisation of the mesh : equirepartirtion of the point between the wall and the image point
	kkm+=1;
      }
      int kks=0;
      for (scalar s in ssol) {  //we start with linear profile for $u$ and $\tilde{\nu}$ to join the wall with the free stream flow
	if (kks%4==0) {  //$u$
	  s[]=(kks/4)*U0/10./(1.*Nsubpoint-1.);
	}
	else if (kks%4==1) {  //$\tilde{\nu}$
	  s[]=(kks/4)*3.*molvis[]/(1.*Nsubpoint-1.);  //turbulent viscosity equals 3 times the molecular viscosity in the free stream, its chi not the turbulent viscosity
	}
	else if (kks%4==2) {  //$(\mu_t+\mu)*\frac{\partial u}{\partial y}
	  s[]=(molvis[]*rho[]+rho[]*(kks/4)*3.*molvis[]/(1.*Nsubpoint-1.)*fv1((kks/4)*3./(1.*Nsubpoint-1.)))*U0/10./(d_IP*(1.*Nsubpoint-1.));
	}
	else {  //$(\nu+\tilde{\nu})\frac{\partial \tilde{\nu}}{\partial y}
	  s[]=(molvis[]+(kks/4)*3.*molvis[]/(1.*Nsubpoint-1.))*3.*molvis[]/(d_IP*(1.*Nsubpoint-1.));
	}
	kks+=1;
      }
    }
  }
}
#endif
	  
    

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



/**
   Condition for the turbulent viscosity  on the solid */
double wall_condition_viscosity(double x2,double y2, scalar s) {

  double cplus;
  double mol_vis;
  double duF;
  double utIP;
  double uboundary;
  double dnuhatF;

  foreach_point(x2,y2) {
    cplus=s.Cplus[];
    mol_vis=molvis[];
    duF=s.duF[];
    utIP=s.utIP[];
    uboundary=s.uboundary[];
    dnuhatF=s.dnuhatF[];
  }
  
  if (fabs(cplus)<1e-30) {
    return 0.;
  }



  double sigmaplus=cplus*pow(cv1*mol_vis,3.);
  double Rea=0.;
  double Reb=cplus/4.;
  double Rem=(Rea+Reb)/2.;
  double resm;
  while (fabs(Rea-Reb)/Reb>1e-5) {
    Rem=(Rea+Reb)/2.;
    resm=function_Relambda(Rem,cplus,sigmaplus);
    if (resm>0.) {
      Reb=Rem;
    }
    else {
      Rea=Rem;
    }
  }
  double discriminant=sq(cplus-2.*Rem)*(1.+8.*Rem/(cplus-4.*Rem));

  
  double viscC=1./2.*(cplus-2.*Rem+sqrt(discriminant));

  return viscC;
}

double wall_condition_viscosity_SA(double x2,double y2, scalar s) {

  double dnuhatF;
  double nuhatIP;

  foreach_point(x2,y2) {
    dnuhatF=s.dnuhatF[];
    nuhatIP=s.nuhatIP[];
  }
  return nuhatIP+dnuhatF*(-d_IP);
}

double wall_condition_viscosity_SA_diss(double x2,double y2, scalar s) {

  double viscSAvalue=0.;

  foreach_point(x2,y2) {
    viscSAvalue=s.viscSAdiss[];
  }

  return viscSAvalue;
}


/**
   We applied the wall model in this event */
event wall_model(i++) {

  if (i>100) {  //(i>100)

    int compteurerror=0;
    int compteurpoint=0;
    scalar pavg[];
    scalar pavg2[];
    scalar pfavg2[];
    vector uavg[];
    vector uavg2[];
    scalar muhatSAavg[];
    scalar muhatSAavg2[];
//We average in time and space the pressure to obtain a less noisy pressure gradient
    foreach() {
      if ((x-Delta>0.)&&(x+Delta<L0)&&(y-Delta>0.)&&(y+Delta<L0)&&(1<0)) {
        pavg[]=(cs[-2,2]*p[-2,2]+4.*cs[-1,2]*p[-1,2]+7.*cs[0,2]*p[0,2]+4.*cs[1,2]*p[1,2]+cs[2,2]*p[2,2]+ \
            4.*cs[-2,1]*p[-2,1]+16.*cs[-1,1]*p[-1,1]+26.*cs[0,1]*p[0,1]+16.*cs[1,1]*p[1,1]+4.*cs[2,1]*p[2,1]+ \
            7.*cs[-2,0]*p[-2,0]+26.*cs[-1,0]*p[-1,0]+41.*cs[0,0]*p[0,0]+26.*cs[1,0]*p[1,0]+7.*cs[2,0]*p[2,0]+ \
            4.*cs[-2,-1]*p[-2,-1]+16.*cs[-1,-1]*p[-1,-1]+26.*cs[0,-1]*p[0,-1]+16.*cs[1,-1]*p[1,-1]+4.*cs[2,-1]*p[2,-1]+ \
            1.*cs[-2,-2]*p[-2,-2]+4.*cs[-1,-2]*p[-1,-2]+7.*cs[0,-2]*p[0,-2]+4.*cs[1,-2]*p[1,-2]+1.*cs[2,-2]*p[2,-2]+SEPS*p[])/ \
          (cs[-2,2]+4.*cs[-1,2]+7.*cs[0,2]+4.*cs[1,2]+cs[2,2]+ \
            4.*cs[-2,1]+16.*cs[-1,1]+26.*cs[0,1]+16.*cs[1,1]+4.*cs[2,1]+ \
            7.*cs[-2,0]+26.*cs[-1,0]+41.*cs[0,0]+26.*cs[1,0]+7.*cs[2,0]+ \
            4.*cs[-2,-1]+16.*cs[-1,-1]+26.*cs[0,-1]+16.*cs[1,-1]+4.*cs[2,-1]+ \
            cs[-2,-2]+4.*cs[-1,-2]+7.*cs[0,-2]+4.*cs[1,-2]+cs[2,-2]+SEPS);//(4.*cs[]*p[]+cs[1]*p[1]+cs[-1]*p[-1]+cs[0,1]*p[0,1]+cs[0,-1]*p[0,-1])/(4.*cs[]+cs[1]+cs[-1]+cs[0,1]+cs[0,-1]+SEPS);
      }
          

      else {
	        pavg[]=(1.*cs[-1,1]*p[-1,1]+2.*cs[0,1]*p[0,1]+1.*cs[1,1]*p[1,1]+ \
            2.*cs[-1,0]*p[-1,0]+4.*cs[0,0]*p[0,0]+2.*cs[1,0]*p[1,0]+ \
            1.*cs[-1,-1]*p[-1,-1]+2.*cs[0,-1]*p[0,-1]+1.*cs[1,-1]*p[1,-1]+SEPS*p[])/ \
          (cs[-1,1]+2.*cs[0,1]+cs[1,1]+ \
            2.*cs[-1,0]+4.*cs[0,0]+2.*cs[1,0]+ \
            cs[-1,-1]+2.*cs[0,-1]+cs[1,-1]+SEPS);//(4.*cs[]*p[]+cs[1]*p[1]+cs[-1]*p[-1]+cs[0,1]*p[0,1]+cs[0,-1]*p[0,-1])/(4.*cs[]+cs[1]+cs[-1]+cs[0,1]+cs[0,-1]+SEPS); //put gaussian filter here too
          
      }
    }
    foreach() {
      pavg2[]=pavg[];
    }
    for (int iii=0;iii<10;iii++) {
      foreach() {
        if ((x-Delta>0.)&&(x+Delta<L0)&&(y-Delta>0.)&&(y+Delta<L0)&&(1<0)) {
          pavg[]=(cs[-2,2]*pavg2[-2,2]+4.*cs[-1,2]*pavg2[-1,2]+7.*cs[0,2]*pavg2[0,2]+4.*cs[1,2]*pavg2[1,2]+cs[2,2]*pavg2[2,2]+ \
        4.*cs[-2,1]*pavg2[-2,1]+16.*cs[-1,1]*pavg2[-1,1]+26.*cs[0,1]*pavg2[0,1]+16.*cs[1,1]*pavg2[1,1]+4.*cs[2,1]*pavg2[2,1]+ \
        7.*cs[-2,0]*pavg2[-2,0]+26.*cs[-1,0]*pavg2[-1,0]+41.*cs[0,0]*pavg2[0,0]+26.*cs[1,0]*pavg2[1,0]+7.*cs[2,0]*pavg2[2,0]+ \
        4.*cs[-2,-1]*pavg2[-2,-1]+16.*cs[-1,-1]*pavg2[-1,-1]+26.*cs[0,-1]*pavg2[0,-1]+16.*cs[1,-1]*pavg2[1,-1]+4.*cs[2,-1]*pavg2[2,-1]+ \
        1.*cs[-2,-2]*pavg2[-2,-2]+4.*cs[-1,-2]*pavg2[-1,-2]+7.*cs[0,-2]*pavg2[0,-2]+4.*cs[1,-2]*pavg2[1,-2]+1.*cs[2,-2]*pavg2[2,-2]+SEPS*pavg2[])/ \
      (cs[-2,2]+4.*cs[-1,2]+7.*cs[0,2]+4.*cs[1,2]+cs[2,2]+ \
        4.*cs[-2,1]+16.*cs[-1,1]+26.*cs[0,1]+16.*cs[1,1]+4.*cs[2,1]+ \
        7.*cs[-2,0]+26.*cs[-1,0]+41.*cs[0,0]+26.*cs[1,0]+7.*cs[2,0]+ \
        4.*cs[-2,-1]+16.*cs[-1,-1]+26.*cs[0,-1]+16.*cs[1,-1]+4.*cs[2,-1]+ \
        cs[-2,-2]+4.*cs[-1,-2]+7.*cs[0,-2]+4.*cs[1,-2]+cs[2,-2]+SEPS);//(4.*cs[]*p[]+cs[1]*p[1]+cs[-1]*p[-1]+cs[0,1]*p[0,1]+cs[0,-1]*p[0,-1])/(4.*cs[]+cs[1]+cs[-1]+cs[0,1]+cs[0,-1]+SEPS);
        }
	    
        else {
	    pavg[]=(1.*cs[-1,1]*pavg2[-1,1]+2.*cs[0,1]*pavg2[0,1]+1.*cs[1,1]*pavg2[1,1]+ \
        2.*cs[-1,0]*pavg2[-1,0]+4.*cs[0,0]*pavg2[0,0]+2.*cs[1,0]*pavg2[1,0]+ \
        1.*cs[-1,-1]*pavg2[-1,-1]+2.*cs[0,-1]*pavg2[0,-1]+1.*cs[1,-1]*pavg2[1,-1]+SEPS*pavg2[])/ \
      (cs[-1,1]+2.*cs[0,1]+cs[1,1]+ \
        2.*cs[-1,0]+4.*cs[0,0]+2.*cs[1,0]+ \
        cs[-1,-1]+2.*cs[0,-1]+cs[1,-1]+SEPS);//(4.*cs[]*pavg2[]+cs[1]*pavg2[1]+cs[-1]*pavg2[-1]+cs[0,1]*pavg2[0,1]+cs[0,-1]*pavg2[0,-1])/(4.*cs[]+cs[1]+cs[-1]+cs[0,1]+cs[0,-1]+SEPS);
      
        }
      }
      foreach() {
        pavg2[]=pavg[];
      }
    }
      
    if (i<=900) {
      foreach() {
        pavg3[]=pavg[];
      }
    }
    else {
      foreach() {
        pavg3[]=pavg3[]*0.9+pavg[]*0.1;
      }
    }
    
#if BVP4
    //We reinitialize the information on the BVP solver
    statssolver.npoint=0;
    statssolver.nitermax=0;
    statssolver.nitermin=100;
    statssolver.niteravg=0.;
    statssolver.norm2max=0.;
    statssolver.norm2avg=0.;
#endif

  
    vector gu[],gv[],gp[];
    foreach() {
      foreach_dimension() {
        gu.x[]=center_gradient(u.x);
        gv.x[]=center_gradient(u.y);
        gp.x[]=center_gradient(p);
      }
    }

    scalar visc[];
    scalar ux[];
    scalar uy[];
    scalar viscSA[];
    scalar viscSAdiss[];

    //data of the previous point to see if we can use the previous solution to start the Newton method
    double dpprec=1.;  
    double advprec=1.;

#if BVP4
    bool initial=false;
    //initial solution and mesh at the point
    double mesh[Nsubpoint];
    double sol[Nsubpoint*Nequation];
#endif

    foreach() {
      double dstar=distance_to_wall(x,y);
      if ((dstar>-Deltamin)&&(dstar<d_IP+Deltamin)&&(x>leading)&&(x<trailing)) {
        compteurpoint+=1;
      
        double nx,ny;
        if ((cs[]>1e-3)&&(cs[]<1-1e-3)) {
          coord n,b;
          embed_geometry(point,&b,&n);
          nx=-n.x;
          ny=-n.y;
        }
        else {
          double dx=(distance_to_wall(x+Delta,y)-distance_to_wall(x-Delta,y))/(2.*Delta);
          double dy=(distance_to_wall(x,y+Delta)-distance_to_wall(x,y-Delta))/(2.*Delta);
          nx=dx/sqrt(sq(dx)+sq(dy)+1e-20);
          ny=dy/sqrt(sq(dx)+sq(dy)+1e-20);
        }

        double tx=ny;
        double ty=-nx;

        if (dstar>Deltamin/sqrt(2.)) {  //close to the wall we seek a better approximation of the wall normal as the results are extremely sensitive to the value at the image point
          double nxp,nyp,nxm,nym;
          nxp=nx;
          nxm=nx;
          nyp=ny;
          nym=ny;
          int i2=sign2(tx);
          int j2=sign2(ty);
          double xcell=x;
          double ycell=y;

          foreach_neighbor(1) {
            if ((fabs(x-(xcell+i2*Delta))<Delta/2.)&&(fabs(y-(ycell+j2*Delta))<Delta/2.)) {
              if ((cs[]>0.)&&(cs[]<1.)) {
                coord n,b;
                embed_geometry(point,&b,&n);
                nxp=-n.x;
                nyp=-n.y;
              }
              else {
                double dxp=(distance_to_wall(xcell+i2*Delta+Delta,ycell+j2*Delta)-distance_to_wall(xcell+i2*Delta-Delta,ycell+j2*Delta))/(2.*Delta);
                double dyp=(distance_to_wall(xcell+i2*Delta,y+j2*Delta+Delta)-distance_to_wall(xcell+i2*Delta,ycell+j2*Delta-Delta))/(2.*Delta);
                nxp=dxp/sqrt(sq(dxp)+sq(dyp)+1e-20);
                nyp=dyp/sqrt(sq(dxp)+sq(dyp)+1e-20);
              }
            }

            if ((fabs(x-(xcell-i2*Delta))<Delta/2.)&&(fabs(y-(ycell-j2*Delta))<Delta/2.)) {
              if ((cs[]>0.)&&(cs[]<1.)) {
                coord n,b;
                embed_geometry(point,&b,&n);
                nxm=-n.x;
                nym=-n.y;
              }
              else {
                double dxm=(distance_to_wall(xcell-i2*Delta+Delta,ycell-j2*Delta)-distance_to_wall(xcell-i2*Delta-Delta,ycell-j2*Delta))/(2.*Delta);
                double dym=(distance_to_wall(xcell-i2*Delta,ycell-j2*Delta+Delta)-distance_to_wall(xcell-i2*Delta,ycell-j2*Delta-Delta))/(2.*Delta);
                nxm=dxm/sqrt(sq(dxm)+sq(dym)+1e-20);
                nym=dym/sqrt(sq(dxm)+sq(dym)+1e-20);
              }
            }
          }

          nx=(nxm+2.*nx+nxp)/4.;
          ny=(nym+2.*ny+nyp)/4.;

          nx=nx/sqrt(sq(nx)+sq(ny));
          ny=ny/sqrt(sq(nx)+sq(ny));

          tx=ny;
          ty=-nx;
        }

        double alpha=d_IP-dstar;
        double x_IP=x+alpha*nx;
        double y_IP=y+alpha*ny;

        double ux_IP=interpolate(u.x,x_IP,y_IP);
        double uy_IP=interpolate(u.y,x_IP,y_IP);
        double muhat_IP=interpolate(muhatSA,x_IP,y_IP);
        double rho_IP=1.;//interpolate(rho,x_IP,y_IP);
        double p_IP=interpolate(p,x_IP,y_IP);
        double pf_IP=interpolate(pf,x_IP,y_IP);
        double ut_IP=ux_IP*tx+uy_IP*ty;
        bool reverse=false;
        if (ut_IP<0.) {
          reverse=true;
          tx=-tx;
          ty=-ty;
          ut_IP=-ut_IP;
        }
        double un_IP=ux_IP*nx+uy_IP*ny;
        double dun_IP=tx*(interpolate(gu.x,x_IP,y_IP)*nx+interpolate(gv.y,x_IP,y_IP)*ny)+ty*(interpolate(gv.x,x_IP,y_IP)*nx+interpolate(gu.y,x_IP,y_IP)*ny);
        double dut_IP=tx*(interpolate(gu.x,x_IP,y_IP)*tx+interpolate(gv.y,x_IP,y_IP)*ty)+ty*(interpolate(gv.x,x_IP,y_IP)*tx+interpolate(gu.y,x_IP,y_IP)*ty);

        //computation of the tangential pressure gradient : pressure is noisy and the value of the gradient can change totally the solution of the BVP so filtering is applied to try to reduce noise. We use a Savitzky-Golay filter of order 4
        //double dp_IP=interpolate(gp.x,x_IP,y_IP)*tx+interpolate(gp.y,x_IP,y_IP)*ty;
      
        double p_IP_1=interpolate(pavg3,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta);
        double p_IP_m1=interpolate(pavg3,x_IP-tx*4.*Delta,y_IP-ty*4.*Delta);
        double p_IP_2=interpolate(pavg3,x_IP+2.*tx*4.*Delta,y_IP+2.*ty*4.*Delta);
        double p_IP_m2=interpolate(pavg3,x_IP-2.*tx*4.*Delta,y_IP-2.*ty*4.*Delta);
        double p_IP_3=interpolate(pavg3,x_IP+3.*tx*4.*Delta,y_IP+3.*ty*4.*Delta);
        double p_IP_m3=interpolate(pavg3,x_IP-3.*tx*4.*Delta,y_IP-3.*ty*4.*Delta);
        double p_IP_4=interpolate(pavg3,x_IP+4.*tx*4.*Delta,y_IP+4.*ty*4.*Delta);
        double p_IP_m4=interpolate(pavg3,x_IP-4.*tx*4.*Delta,y_IP-4.*ty*4.*Delta);

        double dp_IP;
      
        if ((x_IP-4.*tx*4.*Delta>0.)&&(x_IP-4.*tx*4.*Delta<L0)&&(x_IP+4.*tx*4.*Delta>0.)&&(x_IP+4.*tx*4.*Delta<L0)&&(y_IP-4.*ty*4.*Delta>0.)&&(y_IP-4.*ty*4.*Delta<L0)&&(y_IP+4.*ty*4.*Delta>0.)&&(y_IP+4.*ty*4.*Delta<L0)) { //case everything is in the simulation
          if ((interpolate(cs,x_IP+4.*tx*4.*Delta,y_IP+4.*ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-4.*tx*4.*Delta,y_IP-4.*ty*4.*Delta)>0.)) {  //value not in the solid
            dp_IP=(86.*p_IP_m4-142.*p_IP_m3-193.*p_IP_m2-126.*p_IP_m1+126.*p_IP_1+193.*p_IP_2+142.*p_IP_3-86.*p_IP_4)/(1188.*4.*Delta);  //cs_IP=1.
    //p_IP=(15.*p_IP_m4-55.*p_IP_m3+30.*p_IP_m2+135*p_IP_m1+179.*p_IP+135.*p_IP_1+30.*p_IP_2-55.*p_IP_3+15.*p_IP_4)/429.;
          }
          else if ((interpolate(cs,x_IP+3.*tx*4.*Delta,y_IP+3.*ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-3.*tx*4.*Delta,y_IP-3.*ty*4.*Delta)>0.)) {
            dp_IP=(22.*p_IP_m3-67.*p_IP_m2-58.*p_IP_m1+58.*p_IP_1+67.*p_IP_2-22.*p_IP_3)/(252.*4.*Delta);
    //p_IP=(5.*p_IP_m3-30.*p_IP_m2+75.*p_IP_m1+131.*p_IP+75.*p_IP_1-30.*p_IP_2+5.*p_IP_3)/231.;
          }
          else if ((interpolate(cs,x_IP+2.*tx*4.*Delta,y_IP+2.*ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-2.*tx*4.*Delta,y_IP-2.*ty*4.*Delta)>0.)) {
            dp_IP=(p_IP_m2-8.*p_IP_m1+8.*p_IP_1-p_IP_2)/(12.*4.*Delta);
    //p_IP=(-3.*p_IP_m2+12.*p_IP_m1+17.*p_IP+12.*p_IP_1-3.*p_IP_2)/35.;
          }
          else if ((interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-tx*4.*Delta,y_IP-ty*4.*Delta)>0.)) {  //if 2 or -2 value in the solid, switch to a second order gradient
            dp_IP=(p_IP_1-p_IP_m1)/(2.*4.*Delta);
    //p_IP=(p_IP_m1+2.*p_IP+p_IP_1)/4.;
          }
          else if (interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.) {  //gradient biased towards point in the fluid
            dp_IP=(p_IP_1-p_IP)/(4.*Delta);
    //p_IP=(p_IP_1+p_IP)/2.;
          }
          else {
            dp_IP=(p_IP-p_IP_m1)/(4.*Delta);
    //p_IP=(p_IP+p_IP_m1)/2.;
          }
        }
        else if ((x_IP-3.*tx*4.*Delta>0.)&&(x_IP-3.*tx*4.*Delta<L0)&&(x_IP+3.*tx*4.*Delta>0.)&&(x_IP+3.*tx*4.*Delta<L0)&&(y_IP-3.*ty*4.*Delta>0.)&&(y_IP-3.*ty*4.*Delta<L0)&&(y_IP+3.*ty*4.*Delta>0.)&&(y_IP+3.*ty*4.*Delta<L0)) { //case only 3 points are in the simulation
          if ((interpolate(cs,x_IP+3.*tx*4.*Delta,y_IP+3.*ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-3.*tx*4.*Delta,y_IP-3.*ty*4.*Delta)>0.)) {
            dp_IP=(22.*p_IP_m3-67.*p_IP_m2-58.*p_IP_m1+58.*p_IP_1+67.*p_IP_2-22.*p_IP_3)/(252.*4.*Delta);
            //p_IP=(5.*p_IP_m3-30.*p_IP_m2+75.*p_IP_m1+131.*p_IP+75.*p_IP_1-30.*p_IP_2+5.*p_IP_3)/231.;
          }
          else if ((interpolate(cs,x_IP+2.*tx*4.*Delta,y_IP+2.*ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-2.*tx*4.*Delta,y_IP-2.*ty*4.*Delta)>0.)) {
            dp_IP=(p_IP_m2-8.*p_IP_m1+8.*p_IP_1-p_IP_2)/(12.*4.*Delta);
            //p_IP=(-3.*p_IP_m2+12.*p_IP_m1+17.*p_IP+12.*p_IP_1-3.*p_IP_2)/35.;
          }
          else if ((interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-tx*4.*Delta,y_IP-ty*4.*Delta)>0.)) {  //if 2 or -2 value in the solid, switch to a second order gradient
            dp_IP=(p_IP_1-p_IP_m1)/(2.*4.*Delta);
            //p_IP=(p_IP_m1+2.*p_IP+p_IP_1)/4.;
          }
          else if (interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.) {  //gradient biased towards point in the fluid
            dp_IP=(p_IP_1-p_IP)/(4.*Delta);
            //p_IP=(p_IP_1+p_IP)/2.;
          }
          else {
            dp_IP=(p_IP-p_IP_m1)/(4.*Delta);
            //p_IP=(p_IP+p_IP_m1)/2.;
          }
        }
        else if ((x_IP-2.*tx*4.*Delta>0.)&&(x_IP-2.*tx*4.*Delta<L0)&&(x_IP+2.*tx*4.*Delta>0.)&&(x_IP+2.*tx*4.*Delta<L0)&&(y_IP-2.*ty*4.*Delta>0.)&&(y_IP-2.*ty*4.*Delta<L0)&&(y_IP+2.*ty*4.*Delta>0.)&&(y_IP+2.*ty*4.*Delta<L0)) { //case only 2 points are in the simulation
          if ((interpolate(cs,x_IP+2.*tx*4.*Delta,y_IP+2.*ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-2.*tx*4.*Delta,y_IP-2.*ty*4.*Delta)>0.)) {  //value not in the solid
            dp_IP=(p_IP_m2-8.*p_IP_m1+8.*p_IP_1-p_IP_2)/(12.*4.*Delta);
            //p_IP=(-3.*p_IP_m2+12.*p_IP_m1+17.*p_IP+12.*p_IP_1-3.*p_IP_2)/35.;
          }
          else if ((interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-tx*4.*Delta,y_IP-ty*4.*Delta)>0.)) {  //if 2 or -2 value in the solid, switch to a second order gradient
            dp_IP=(p_IP_1-p_IP_m1)/(2.*4.*Delta);
            //p_IP=(p_IP_m1+2.*p_IP+p_IP_1)/4.;
          }
          else if (interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.) {  //gradient biased towards point in the fluid
            dp_IP=(p_IP_1-p_IP)/(4.*Delta);
            //p_IP=(p_IP_1+p_IP)/2.;
          }
          else {
            dp_IP=(p_IP-p_IP_m1)/(4.*Delta);
    //p_IP=(p_IP+p_IP_m1)/2.;
          }
        }
        else {
          if ((interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-tx*4.*Delta,y_IP-ty*4.*Delta)>0.)) {  //if 2 or -2 not in the simulation box, switch to a second order gradient
            dp_IP=(p_IP_1-p_IP_m1)/(2.*4.*Delta);
    //p_IP=(p_IP_m1+2.*p_IP+p_IP_1)/4.;
          }
          else if (interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.) {  //gradient biased towards point in the fluid
            dp_IP=(p_IP_1-p_IP)/(4.*Delta);
    //p_IP=(p_IP_1+p_IP)/2.;
          }
          else {
            dp_IP=(p_IP-p_IP_m1)/(4.*Delta);
    //p_IP=(p_IP+p_IP_m1)/2.;
          }
        }


        double ut_IP_1=interpolate(uavg.x,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)*tx+interpolate(uavg.y,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)*ty;
double ut_IP_m1=interpolate(uavg.x,x_IP-tx*4.*Delta,y_IP-ty*4.*Delta)*tx+interpolate(uavg.y,x_IP-tx*4.*Delta,y_IP-ty*4.*Delta)*ty;
        double ut_IP_2=interpolate(uavg.x,x_IP+2.*tx*4.*Delta,y_IP+2.*ty*4.*Delta)*tx+interpolate(uavg.y,x_IP+2.*tx*4.*Delta,y_IP+2.*ty*4.*Delta)*ty;
        double ut_IP_m2=interpolate(uavg.x,x_IP-2.*tx*4.*Delta,y_IP-2.*ty*4.*Delta)*tx+interpolate(uavg.y,x_IP-2.*tx*4.*Delta,y_IP-2.*ty*4.*Delta)*ty;
        double ut_IP_3=interpolate(uavg.x,x_IP+3.*tx*4.*Delta,y_IP+3.*ty*4.*Delta)*tx+interpolate(uavg.y,x_IP+3.*tx*4.*Delta,y_IP+3.*ty*4.*Delta)*ty;
        double ut_IP_m3=interpolate(uavg.x,x_IP-3.*tx*4.*Delta,y_IP-3.*ty*4.*Delta)*tx+interpolate(uavg.y,x_IP-3.*tx*4.*Delta,y_IP-3.*ty*4.*Delta)*ty;
        double ut_IP_4=interpolate(uavg.x,x_IP+4.*tx*4.*Delta,y_IP+4.*ty*4.*Delta)*tx+interpolate(uavg.y,x_IP+4.*tx*4.*Delta,y_IP+4.*ty*4.*Delta)*ty;
        double ut_IP_m4=interpolate(uavg.x,x_IP-4.*tx*4.*Delta,y_IP-4.*ty*4.*Delta)*tx+interpolate(uavg.y,x_IP-4.*tx*4.*Delta,y_IP-4.*ty*4.*Delta)*ty;
      
        if ((x_IP-4.*tx*4.*Delta>0.)&&(x_IP-4.*tx*4.*Delta<L0)&&(x_IP+4.*tx*4.*Delta>0.)&&(x_IP+4.*tx*4.*Delta<L0)&&(y_IP-4.*ty*4.*Delta>0.)&&(y_IP-4.*ty*4.*Delta<L0)&&(y_IP+4.*ty*4.*Delta>0.)&&(y_IP+4.*ty*4.*Delta<L0)) { //case everything is in the simulation
          if ((interpolate(cs,x_IP+4.*tx*4.*Delta,y_IP+4.*ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-4.*tx*4.*Delta,y_IP-4.*ty*4.*Delta)>0.)) {  //value not in the solid
            dut_IP=(86.*ut_IP_m4-142.*ut_IP_m3-193.*ut_IP_m2-126.*ut_IP_m1+126.*ut_IP_1+193.*ut_IP_2+142.*ut_IP_3-86.*ut_IP_4)/(1188.*4.*Delta);  //cs_IP=1.
            //ut_IP=(15.*ut_IP_m4-55.*ut_IP_m3+30.*ut_IP_m2+135.*ut_IP_m1+179.*ut_IP+135.*ut_IP_1+30.*ut_IP_2-55.*ut_IP_3+15.*ut_IP_4)/429.;
          }
          else if ((interpolate(cs,x_IP+3.*tx*4.*Delta,y_IP+3.*ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-3.*tx*4.*Delta,y_IP-3.*ty*4.*Delta)>0.)) {
            dut_IP=(22.*ut_IP_m3-67.*ut_IP_m2-58.*ut_IP_m1+58.*ut_IP_1+67.*ut_IP_2-22.*ut_IP_3)/(252.*4.*Delta);
            //ut_IP=(5.*ut_IP_m3-30.*ut_IP_m2+75.*ut_IP_m1+131.*ut_IP+75.*ut_IP_1-30.*ut_IP_2+5.*ut_IP_3)/231.;
          }
          else if ((interpolate(cs,x_IP+2.*tx*4.*Delta,y_IP+2.*ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-2.*tx*4.*Delta,y_IP-2.*ty*4.*Delta)>0.)) {
            dut_IP=(ut_IP_m2-8.*ut_IP_m1+8.*ut_IP_1-ut_IP_2)/(12.*4.*Delta);
    //ut_IP=(-3.*ut_IP_m2+12.*ut_IP_m1+17.*ut_IP+12.*ut_IP_1-3.*ut_IP_2)/35.;
          }
          else if ((interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-tx*4.*Delta,y_IP-ty*4.*Delta)>0.)) {  //if 2 or -2 value in the solid, switch to a second order gradient
            dut_IP=(ut_IP_1-ut_IP_m1)/(2.*4.*Delta);
            //ut_IP=(ut_IP_m1+2.*ut_IP+ut_IP_1)/4.;
          }
          else if (interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.) {  //gradient biased towards point in the fluid
            dut_IP=(ut_IP_1-ut_IP)/(4.*Delta);
            //ut_IP=(ut_IP_1+ut_IP)/2.;
          }
          else {
            dut_IP=(ut_IP-ut_IP_m1)/(4.*Delta);
            //ut_IP=(ut_IP+ut_IP_m1)/2.;
          }
        }
        else if ((x_IP-3.*tx*4.*Delta>0.)&&(x_IP-3.*tx*4.*Delta<L0)&&(x_IP+3.*tx*4.*Delta>0.)&&(x_IP+3.*tx*4.*Delta<L0)&&(y_IP-3.*ty*4.*Delta>0.)&&(y_IP-3.*ty*4.*Delta<L0)&&(y_IP+3.*ty*4.*Delta>0.)&&(y_IP+3.*ty*4.*Delta<L0)) { //case only 3 points are in the simulation
          if ((interpolate(cs,x_IP+3.*tx*4.*Delta,y_IP+3.*ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-3.*tx*4.*Delta,y_IP-3.*ty*4.*Delta)>0.)) {
            dut_IP=(22.*ut_IP_m3-67.*ut_IP_m2-58.*ut_IP_m1+58.*ut_IP_1+67.*ut_IP_2-22.*ut_IP_3)/(252.*4.*Delta);
            //ut_IP=(5.*ut_IP_m3-30.*ut_IP_m2+75.*ut_IP_m1+131.*ut_IP+75.*ut_IP_1-30.*ut_IP_2+5.*ut_IP_3)/231.;
          }
          else if ((interpolate(cs,x_IP+2.*tx*4.*Delta,y_IP+2.*ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-2.*tx*4.*Delta,y_IP-2.*ty*4.*Delta)>0.)) {
            dut_IP=(ut_IP_m2-8.*ut_IP_m1+8.*ut_IP_1-ut_IP_2)/(12.*4.*Delta);
            //ut_IP=(-3.*ut_IP_m2+12.*ut_IP_m1+17.*ut_IP+12.*ut_IP_1-3.*ut_IP_2)/35.;
          }
          else if ((interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-tx*4.*Delta,y_IP-ty*4.*Delta)>0.)) {  //if 2 or -2 value in the solid, switch to a second order gradient
            dut_IP=(ut_IP_1-ut_IP_m1)/(2.*4.*Delta);
            //ut_IP=(ut_IP_m1+2.*ut_IP+ut_IP_1)/4.;
          }
          else if (interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.) {  //gradient biased towards point in the fluid
            dut_IP=(ut_IP_1-ut_IP)/(4.*Delta);
            //ut_IP=(ut_IP_1+ut_IP)/2.;
          }
          else {
            dut_IP=(ut_IP-ut_IP_m1)/(4.*Delta);
            //ut_IP=(ut_IP+ut_IP_m1)/2.;
          }
        }
        else if ((x_IP-2.*tx*4.*Delta>0.)&&(x_IP-2.*tx*4.*Delta<L0)&&(x_IP+2.*tx*4.*Delta>0.)&&(x_IP+2.*tx*4.*Delta<L0)&&(y_IP-2.*ty*4.*Delta>0.)&&(y_IP-2.*ty*4.*Delta<L0)&&(y_IP+2.*ty*4.*Delta>0.)&&(y_IP+2.*ty*4.*Delta<L0)) { //case only 2 points are in the simulation
          if ((interpolate(cs,x_IP+2.*tx*4.*Delta,y_IP+2.*ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-2.*tx*4.*Delta,y_IP-2.*ty*4.*Delta)>0.)) {  //value not in the solid
            dut_IP=(ut_IP_m2-8.*ut_IP_m1+8.*ut_IP_1-ut_IP_2)/(12.*4.*Delta);
            //ut_IP=(-3.*ut_IP_m2+12.*ut_IP_m1+17.*ut_IP+12.*ut_IP_1-3.*ut_IP_2)/35.;
          }
          else if ((interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-tx*4.*Delta,y_IP-ty*4.*Delta)>0.)) {  //if 2 or -2 value in the solid, switch to a second order gradient
            dut_IP=(ut_IP_1-ut_IP_m1)/(2.*4.*Delta);
            //ut_IP=(ut_IP_m1+2.*ut_IP+ut_IP_1)/4.;
          }
          else if (interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.) {  //gradient biased towards point in the fluid
            dut_IP=(ut_IP_1-ut_IP)/(4.*Delta);
    //ut_IP=(ut_IP_1+ut_IP)/2.;
          }
          else {
            dut_IP=(ut_IP-ut_IP_m1)/(4.*Delta);
    //ut_IP=(ut_IP+ut_IP_m1)/2.;
          }
        }
        else {
          if ((interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-tx*4.*Delta,y_IP-ty*4.*Delta)>0.)) {  //if 2 or -2 not in the simulation box, switch to a second order gradient
            dut_IP=(ut_IP_1-ut_IP_m1)/(2.*4.*Delta);
            //ut_IP=(ut_IP_m1+2.*ut_IP+ut_IP_1)/4.;
          }
          else if (interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.) {  //gradient biased towards point in the fluid
            dut_IP=(ut_IP_1-ut_IP)/(4.*Delta);
            //ut_IP=(ut_IP_1+ut_IP)/2.;
          }
          else {
            dut_IP=(ut_IP-ut_IP_m1)/(4.*Delta);
            //ut_IP=(ut_IP+ut_IP_m1)/2.;
          }
        } 

        double un_IP_1=interpolate(uavg.x,x_IP+nx*Delta,y_IP+ny*Delta)*tx+interpolate(uavg.y,x_IP+nx*Delta,y_IP+ny*Delta)*ty;
        double un_IP_m1=interpolate(uavg.x,x_IP-nx*Delta,y_IP-ny*Delta)*tx+interpolate(uavg.y,x_IP-nx*Delta,y_IP-ny*Delta)*ty;
        double un_IP_2=interpolate(uavg.x,x_IP+2.*nx*Delta,y_IP+2.*ny*Delta)*tx+interpolate(uavg.y,x_IP+2.*nx*Delta,y_IP+2.*ny*Delta)*ty;
        double un_IP_m2=interpolate(uavg.x,x_IP-2.*nx*Delta,y_IP-2.*ny*Delta)*tx+interpolate(uavg.y,x_IP-2.*nx*Delta,y_IP-2.*ny*Delta)*ty;
        double un_IP_3=interpolate(uavg.x,x_IP+3.*nx*Delta,y_IP+3.*ny*Delta)*tx+interpolate(uavg.y,x_IP+3.*nx*Delta,y_IP+3.*ny*Delta)*ty;
        double un_IP_m3=interpolate(uavg.x,x_IP-3.*nx*Delta,y_IP-3.*ny*Delta)*tx+interpolate(uavg.y,x_IP-3.*nx*Delta,y_IP-3.*ny*Delta)*ty;
        double un_IP_4=interpolate(uavg.x,x_IP+4.*nx*Delta,y_IP+4.*ny*Delta)*tx+interpolate(uavg.y,x_IP+4.*nx*Delta,y_IP+4.*ny*Delta)*ty;
        double un_IP_m4=interpolate(uavg.x,x_IP-4.*nx*Delta,y_IP-4.*ny*Delta)*tx+interpolate(uavg.y,x_IP-4.*nx*Delta,y_IP-4.*ny*Delta)*ty;
      
        if ((x_IP-4.*nx*Delta>0.)&&(x_IP-4.*nx*Delta<L0)&&(x_IP+4.*nx*Delta>0.)&&(x_IP+4.*nx*Delta<L0)&&(y_IP-4.*ny*Delta>0.)&&(y_IP-4.*ny*Delta<L0)&&(y_IP+4.*ny*Delta>0.)&&(y_IP+4.*ny*Delta<L0)) { //case everything is in the simulation
          if ((interpolate(cs,x_IP+4.*nx*Delta,y_IP+4.*ny*Delta)>0.)&&(interpolate(cs,x_IP-4.*nx*Delta,y_IP-4.*ny*Delta)>0.)) {  //value not in the solid
            dun_IP=(86.*un_IP_m4-142.*un_IP_m3-193.*un_IP_m2-126.*un_IP_m1+126.*un_IP_1+193.*un_IP_2+142.*un_IP_3-86.*un_IP_4)/(1188.*Delta);  //cs_IP=1.
          }
          else if ((interpolate(cs,x_IP+3.*nx*Delta,y_IP+3.*ny*Delta)>0.)&&(interpolate(cs,x_IP-3.*nx*Delta,y_IP-3.*ny*Delta)>0.)) {
            dun_IP=(22.*un_IP_m3-67.*un_IP_m2-58.*un_IP_m1+58.*un_IP_1+67.*un_IP_2-22.*un_IP_3)/(252.*Delta);
          }
          else if ((interpolate(cs,x_IP+2.*nx*Delta,y_IP+2.*ny*Delta)>0.)&&(interpolate(cs,x_IP-2.*nx*Delta,y_IP-2.*ny*Delta)>0.)) {
            dun_IP=(un_IP_m2-8.*un_IP_m1+8.*un_IP_1-un_IP_2)/(12.*Delta);
          }
          else if ((interpolate(cs,x_IP+nx*Delta,y_IP+ny*Delta)>0.)&&(interpolate(cs,x_IP-nx*Delta,y_IP-ny*Delta)>0.)) {  //if 2 or -2 value in the solid, switch to a second order gradient
            dun_IP=(un_IP_1-un_IP_m1)/(2.*Delta);
          }
          else if (interpolate(cs,x_IP+nx*Delta,y_IP+ny*Delta)>0.) {  //gradient biased towards point in the fluid
            dun_IP=(un_IP_1-un_IP)/(Delta);
          }
          else {
            dun_IP=(un_IP-un_IP_m1)/(Delta);
          }
        }
        else if ((x_IP-3.*nx*Delta>0.)&&(x_IP-3.*nx*Delta<L0)&&(x_IP+3.*nx*Delta>0.)&&(x_IP+3.*nx*Delta<L0)&&(y_IP-3.*ny*Delta>0.)&&(y_IP-3.*ny*Delta<L0)&&(y_IP+3.*ny*Delta>0.)&&(y_IP+3.*ny*Delta<L0)) { //case only 3 points are in the simulation
          if ((interpolate(cs,x_IP+3.*nx*Delta,y_IP+3.*ny*Delta)>0.)&&(interpolate(cs,x_IP-3.*nx*Delta,y_IP-3.*ny*Delta)>0.)) {
            dun_IP=(22.*un_IP_m3-67.*un_IP_m2-58.*un_IP_m1+58.*un_IP_1+67.*un_IP_2-22.*un_IP_3)/(252.*Delta);
          }
          else if ((interpolate(cs,x_IP+2.*nx*Delta,y_IP+2.*ny*Delta)>0.)&&(interpolate(cs,x_IP-2.*nx*Delta,y_IP-2.*ny*Delta)>0.)) {
            dun_IP=(un_IP_m2-8.*un_IP_m1+8.*un_IP_1-un_IP_2)/(12.*Delta);
          }
          else if ((interpolate(cs,x_IP+nx*Delta,y_IP+ny*Delta)>0.)&&(interpolate(cs,x_IP-nx*Delta,y_IP-ny*Delta)>0.)) {  //if 2 or -2 value in the solid, switch to a second order gradient
            dun_IP=(un_IP_1-un_IP_m1)/(2.*Delta);
          }
          else if (interpolate(cs,x_IP+nx*Delta,y_IP+ny*Delta)>0.) {  //gradient biased towards point in the fluid
            dun_IP=(un_IP_1-un_IP)/(Delta);
          }
          else {
            dun_IP=(un_IP-un_IP_m1)/(Delta);
          }
        }
        else if ((x_IP-2.*nx*Delta>0.)&&(x_IP-2.*nx*Delta<L0)&&(x_IP+2.*nx*Delta>0.)&&(x_IP+2.*nx*Delta<L0)&&(y_IP-2.*ny*Delta>0.)&&(y_IP-2.*ny*Delta<L0)&&(y_IP+2.*ny*Delta>0.)&&(y_IP+2.*ny*Delta<L0)) { //case only 2 points are in the simulation
          if ((interpolate(cs,x_IP+2.*nx*Delta,y_IP+2.*ny*Delta)>0.)&&(interpolate(cs,x_IP-2.*nx*Delta,y_IP-2.*ny*Delta)>0.)) {  //value not in the solid
            dun_IP=(un_IP_m2-8.*un_IP_m1+8.*un_IP_1-un_IP_2)/(12.*Delta);
          }
          else if ((interpolate(cs,x_IP+nx*Delta,y_IP+ny*Delta)>0.)&&(interpolate(cs,x_IP-nx*Delta,y_IP-ny*Delta)>0.)) {  //if 2 or -2 value in the solid, switch to a second order gradient
            dun_IP=(un_IP_1-un_IP_m1)/(2.*Delta);
          }
          else if (interpolate(cs,x_IP+nx*Delta,y_IP+ny*Delta)>0.) {  //gradient biased towards point in the fluid
            dun_IP=(un_IP_1-un_IP)/(Delta);
          }
          else {
            dun_IP=(un_IP-un_IP_m1)/(Delta);
          }
        }
        else {
          if ((interpolate(cs,x_IP+nx*Delta,y_IP+ny*Delta)>0.)&&(interpolate(cs,x_IP-nx*Delta,y_IP-ny*Delta)>0.)) {  //if 2 or -2 not in the simulation box, switch to a second order gradient
            dun_IP=(un_IP_1-un_IP_m1)/(2.*Delta);
          }
          else if (interpolate(cs,x_IP+nx*Delta,y_IP+ny*Delta)>0.) {  //gradient biased towards point in the fluid
            dun_IP=(un_IP_1-un_IP)/(Delta);
          }
          else {
            dun_IP=(un_IP-un_IP_m1)/(Delta);
          }
        }
/** 
    double muhat_IP_1=interpolate(muhatSA,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta);
    double muhat_IP_m1=interpolate(muhatSA,x_IP-tx*4.*Delta,y_IP-ty*4.*Delta);
    double muhat_IP_2=interpolate(muhatSA,x_IP+2.*tx*4.*Delta,y_IP+2.*ty*4.*Delta);
    double muhat_IP_m2=interpolate(muhatSA,x_IP-2.*tx*4.*Delta,y_IP-2.*ty*4.*Delta);
    double muhat_IP_3=interpolate(muhatSA,x_IP+3.*tx*4.*Delta,y_IP+3.*ty*4.*Delta);
    double muhat_IP_m3=interpolate(muhatSA,x_IP-3.*tx*4.*Delta,y_IP-3.*ty*4.*Delta);
    double muhat_IP_4=interpolate(muhatSA,x_IP+4.*tx*4.*Delta,y_IP+4.*ty*4.*Delta);
    double muhat_IP_m4=interpolate(muhatSA,x_IP-4.*tx*4.*Delta,y_IP-4.*ty*4.*Delta);
      
    if ((x_IP-4.*tx*4.*Delta>0.)&&(x_IP-4.*tx*4.*Delta<L0)&&(x_IP+4.*tx*4.*Delta>0.)&&(x_IP+4.*tx*4.*Delta<L0)&&(y_IP-4.*ty*4.*Delta>0.)&&(y_IP-4.*ty*4.*Delta<L0)&&(y_IP+4.*ty*4.*Delta>0.)&&(y_IP+4.*ty*4.*Delta<L0)) { //case everything is in the simulation
      if ((interpolate(cs,x_IP+4.*tx*4.*Delta,y_IP+4.*ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-4.*tx*4.*Delta,y_IP-4.*ty*4.*Delta)>0.)) {  //value not in the solid
        muhat_IP=(15.*muhat_IP_m4-55.*muhat_IP_m3+30.*muhat_IP_m2+135.*muhat_IP_m1+179.*muhat_IP+135.*muhat_IP_1+30.*muhat_IP_2-55.*muhat_IP_3+15.*muhat_IP_4)/429.;
      }
      else if ((interpolate(cs,x_IP+3.*tx*4.*Delta,y_IP+3.*ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-3.*tx*4.*Delta,y_IP-3.*ty*4.*Delta)>0.)) {
        muhat_IP=(5.*muhat_IP_m3-30.*muhat_IP_m2+75.*muhat_IP_m1+131.*muhat_IP+75.*muhat_IP_1-30.*muhat_IP_2+5.*muhat_IP_3)/231.;
      }
      else if ((interpolate(cs,x_IP+2.*tx*4.*Delta,y_IP+2.*ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-2.*tx*4.*Delta,y_IP-2.*ty*4.*Delta)>0.)) {
        muhat_IP=(-3.*muhat_IP_m2+12.*muhat_IP_m1+17.*muhat_IP+12.*muhat_IP_1-3.*muhat_IP_2)/35.;
      }
      else if ((interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-tx*4.*Delta,y_IP-ty*4.*Delta)>0.)) {  //if 2 or -2 value in the solid, switch to a second order gradient
        muhat_IP=(muhat_IP_m1+2.*muhat_IP+muhat_IP_1)/4.;
      }
      else if (interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.) {  //gradient biased towards point in the fluid
        muhat_IP=(muhat_IP_1+muhat_IP)/2.;
      }
      else {
        muhat_IP=(muhat_IP+muhat_IP_m1)/2.;
      }
    }
    else if ((x_IP-3.*tx*4.*Delta>0.)&&(x_IP-3.*tx*4.*Delta<L0)&&(x_IP+3.*tx*4.*Delta>0.)&&(x_IP+3.*tx*4.*Delta<L0)&&(y_IP-3.*ty*4.*Delta>0.)&&(y_IP-3.*ty*4.*Delta<L0)&&(y_IP+3.*ty*4.*Delta>0.)&&(y_IP+3.*ty*4.*Delta<L0)) { //case only 3 points are in the simulation
      if ((interpolate(cs,x_IP+3.*tx*4.*Delta,y_IP+3.*ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-3.*tx*4.*Delta,y_IP-3.*ty*4.*Delta)>0.)) {
        muhat_IP=(5.*muhat_IP_m3-30.*muhat_IP_m2+75.*muhat_IP_m1+131.*muhat_IP+75.*muhat_IP_1-30.*muhat_IP_2+5.*muhat_IP_3)/231.;
      }
      else if ((interpolate(cs,x_IP+2.*tx*4.*Delta,y_IP+2.*ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-2.*tx*4.*Delta,y_IP-2.*ty*4.*Delta)>0.)) {
        muhat_IP=(-3.*muhat_IP_m2+12.*muhat_IP_m1+17.*muhat_IP+12.*muhat_IP_1-3.*muhat_IP_2)/35.;
      }
      else if ((interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-tx*4.*Delta,y_IP-ty*4.*Delta)>0.)) {  //if 2 or -2 value in the solid, switch to a second order gradient
        muhat_IP=(muhat_IP_m1+2.*muhat_IP+muhat_IP_1)/4.;
      }
      else if (interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.) {  //gradient biased towards point in the fluid
        muhat_IP=(muhat_IP_1+muhat_IP)/2.;
      }
      else {
        muhat_IP=(muhat_IP+muhat_IP_m1)/2.;
      }
    }
    else if ((x_IP-2.*tx*4.*Delta>0.)&&(x_IP-2.*tx*4.*Delta<L0)&&(x_IP+2.*tx*4.*Delta>0.)&&(x_IP+2.*tx*4.*Delta<L0)&&(y_IP-2.*ty*4.*Delta>0.)&&(y_IP-2.*ty*4.*Delta<L0)&&(y_IP+2.*ty*4.*Delta>0.)&&(y_IP+2.*ty*4.*Delta<L0)) { //case only 2 points are in the simulation
      if ((interpolate(cs,x_IP+2.*tx*4.*Delta,y_IP+2.*ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-2.*tx*4.*Delta,y_IP-2.*ty*4.*Delta)>0.)) {  //value not in the solid
        muhat_IP=(-3.*muhat_IP_m2+12.*muhat_IP_m1+17.*muhat_IP+12.*muhat_IP_1-3.*muhat_IP_2)/35.;
      }
      else if ((interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-tx*4.*Delta,y_IP-ty*4.*Delta)>0.)) {  //if 2 or -2 value in the solid, switch to a second order gradient
        muhat_IP=(muhat_IP_m1+2.*muhat_IP+muhat_IP_1)/4.;
      }
      else if (interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.) {  //gradient biased towards point in the fluid
        muhat_IP=(muhat_IP_1+muhat_IP)/2.;
      }
      else {
        muhat_IP=(muhat_IP+muhat_IP_m1)/2.;
      }
    }
    else {
      if ((interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.)&&(interpolate(cs,x_IP-tx*4.*Delta,y_IP-ty*4.*Delta)>0.)) {  //if 2 or -2 not in the simulation box, switch to a second order gradient
        muhat_IP=(muhat_IP_m1+2.*muhat_IP+muhat_IP_1)/4.;
      }
      else if (interpolate(cs,x_IP+tx*4.*Delta,y_IP+ty*4.*Delta)>0.) {  //gradient biased towards point in the fluid
        muhat_IP=(muhat_IP_1+muhat_IP)/2.;
      }
      else {
        muhat_IP=(muhat_IP+muhat_IP_m1)/2.;
      }
    } */



        //Initial guess for utau, needed for the blending function
double utau=newton_utau(muhat.utau[] ? muhat.utau[] : 0.1,ut_IP,d_IP,molvis[]);

        muhat.utIP[]=ut_IP;

#if BVP4   

        //parameter for the function
        Param param;
        dp_IP=fabs(dp_IP)>100. ? 0. : dp_IP;
        param.dp=dp_IP;
        param.adv=ut_IP*dut_IP+un_IP*dun_IP;
        param.yF=d_IP;
        param.uF=ut_IP;
        param.nuhatF=muhat_IP/rho_IP;
        param.utau=utau;
        param.rho=1.;//rho[];
        param.mu=molvis[];//mu[];
        muhat.utauprev[]=utau;
      
        bool findneighbor=true;  //needed if the Newton did not converged

        if ((i>900)&&(x>leading+Deltamin)) {
	  
          if (i<=901) {
            if (initial) {
              double meshinit[Nsubpoint];
              double solinit[Nsubpoint*Nequation];
              bool initneighbor=false;
              foreach_neighbor(1) {
                if (updatedbvp[]>=0.5) {
                  initneighbor=true;
                  int kk1=0;
                  for (scalar s in smesh) {
                    meshinit[kk1]=s[];
                    kk1+=1;
                  }
                  int kk2=0;
                  for (scalar s in ssol) {
                    solinit[kk2]=s[];
                    kk2+=1; 
                  }
                }
              }
              double normphiinit=1e10;
              if (initneighbor) {
                double Phiinit[Nequation*Nsubpoint];
                constPhi(Nsubpoint,Nequation,Nboundary0,1,Phiinit,meshinit,solinit,0,param);
                normphiinit=0.;
                for (int pp=0;pp<Nequation*Nsubpoint;pp++) {
                  normphiinit+=sq(Phiinit[pp]);
                }
                normphiinit=sqrt(normphiinit);
              }
	    
              double meshpart[Nsubpoint];
              double solpart[Nsubpoint*Nequation];
              int kk1=0;
              for (scalar s in smesh) {
                if (kk1<25) {
                  meshpart[kk1]=d_IP*kk1/(15.*Nsubpoint);
                }
                else {
                  meshpart[kk1]=d_IP*24./(15.*Nsubpoint)+(d_IP-d_IP*24./(15.*Nsubpoint))/(1.*Nsubpoint-26.)*(kk1-25.);
                }

                kk1+=1;
              }
              int kk2=0;
              for (scalar s in ssol) {

                if (kk2%Nequation==0) {
                  solpart[kk2]=ut_IP*meshpart[kk2/Nequation]/d_IP;
                }
                else if (kk2%Nequation==1) {
                  solpart[kk2]=muhat_IP*meshpart[kk2/Nequation]/d_IP;
                }
                else if (kk2%Nequation==2) {
                  solpart[kk2]=(param.mu+param.rho*muhat_IP*meshpart[kk2/Nequation]/d_IP*fv1(muhat_IP*meshpart[kk2/Nequation]/d_IP/molvis[]))*ut_IP/d_IP;
                }
                else {
                  solpart[kk2]=(param.mu/param.rho+muhat_IP*meshpart[kk2/Nequation]/d_IP)*muhat_IP/d_IP;
                }
                kk2+=1; 
              }
	    
              double Phipart[Nequation*Nsubpoint];
              constPhi(Nsubpoint,Nequation,Nboundary0,1,Phipart,meshpart,solpart,10,param);
              double normphipart=0.;
              for (int pp=0;pp<Nequation*Nsubpoint;pp++) {
                normphipart+=sq(Phipart[pp]);
              }
              normphipart=sqrt(normphipart);
	    
              if (normphiinit>normphipart) {
                for (int kk1=0;kk1<Nsubpoint;kk1++) {
                  mesh[kk1]=meshpart[kk1];
                }
                for (int kk2=0;kk2<Nequation*Nsubpoint;kk2++) {
                  sol[kk2]=solpart[kk2];
                }
                initial=false;
              }
              else {
                for (int kk1=0;kk1<Nsubpoint;kk1++) {
                  mesh[kk1]=meshinit[kk1];
                }
                for (int kk2=0;kk2<Nequation*Nsubpoint;kk2++) {
                  sol[kk2]=solinit[kk2];
                }
                initial=true;
              }
            }
	  
            if (!initial) {
              int kk1=0;
              for (scalar s in smesh) {
                if (kk1<25) {
                  mesh[kk1]=d_IP*kk1/(15.*Nsubpoint);
                }
                else {
                  mesh[kk1]=d_IP*24./(15.*Nsubpoint)+(d_IP-d_IP*24./(15.*Nsubpoint))/(1.*Nsubpoint-26.)*(kk1-25.);
                }

                kk1+=1;
              }
              int kk2=0;
              for (scalar s in ssol) {

                if (kk2%Nequation==0) {
                  sol[kk2]=ut_IP*mesh[kk2/Nequation]/d_IP;
                }
                else if (kk2%Nequation==1) {
                  sol[kk2]=muhat_IP*mesh[kk2/Nequation]/d_IP;
                }
                else if (kk2%Nequation==2) {
                  sol[kk2]=(param.mu+param.rho*muhat_IP*mesh[kk2/Nequation]/d_IP*fv1(muhat_IP*mesh[kk2/Nequation]/d_IP/molvis[]))*ut_IP/d_IP;
                }
                else {
                  sol[kk2]=(param.mu/param.rho+muhat_IP*mesh[kk2/Nequation]/d_IP)*muhat_IP/d_IP;
                }

                kk2+=1; 
              }
            }
          }
          else {
	  
            double meshinit[Nsubpoint];
            double solinit[Nsubpoint*Nequation];
            bool initneighbor=false;
            foreach_neighbor(1) {
              if (updatedbvp[]>=0.5) {
                initneighbor=true;
                int kk1=0;
                for (scalar s in smesh) {
                  meshinit[kk1]=s[];
                  kk1+=1;
                }
                int kk2=0;
                for (scalar s in ssol) {
                  solinit[kk2]=s[];

                  kk2+=1; 
                }

              }
            }
            double normphiinit=1e10;
            if (initneighbor) {
              double Phiinit[Nequation*Nsubpoint];
              constPhi(Nsubpoint,Nequation,Nboundary0,1,Phiinit,meshinit,solinit,0,param);
              normphiinit=0.;
              for (int pp=0;pp<Nequation*Nsubpoint;pp++) {
                normphiinit+=sq(Phiinit[pp]);
              }
              normphiinit=sqrt(normphiinit);
            }
	  
            double meshpart[Nsubpoint];
            double solpart[Nsubpoint*Nequation];
            int kk1=0;
            for (scalar s in smesh) {
              if (kk1<25) {
                meshpart[kk1]=d_IP*kk1/(15.*Nsubpoint);
              }
              else {
                meshpart[kk1]=d_IP*24./(15.*Nsubpoint)+(d_IP-d_IP*24./(15.*Nsubpoint))/(1.*Nsubpoint-26.)*(kk1-25.);
              }
              kk1+=1;
            }
            int kk2=0;
            for (scalar s in ssol) {
              if (kk2%Nequation==0) {
                solpart[kk2]=ut_IP*meshpart[kk2/Nequation]/d_IP;
              }
              else if (kk2%Nequation==1) {
                solpart[kk2]=muhat_IP*meshpart[kk2/Nequation]/d_IP;
              }
              else if (kk2%Nequation==2) {
                solpart[kk2]=(param.mu+param.rho*muhat_IP*meshpart[kk2/Nequation]/d_IP*fv1(muhat_IP*meshpart[kk2/Nequation]/d_IP/molvis[]))*ut_IP/d_IP;
              }
              else {
                solpart[kk2]=(param.mu/param.rho+muhat_IP*meshpart[kk2/Nequation]/d_IP)*muhat_IP/d_IP;
              }
              kk2+=1; 
            }
	  
            double Phipart[Nequation*Nsubpoint];
            constPhi(Nsubpoint,Nequation,Nboundary0,1,Phipart,meshpart,solpart,10,param);
            double normphipart=0.;
            for (int pp=0;pp<Nequation*Nsubpoint;pp++) {
              normphipart+=sq(Phipart[pp]);
            }
            normphipart=sqrt(normphipart);

            double meshprev[Nsubpoint];
            double solprev[Nsubpoint*Nequation];
            kk1=0;
            for (scalar s in smesh) {
              meshprev[kk1]=s[];
              kk1+=1;
            }
            kk2=0;
            for (scalar s in ssol) {
              solprev[kk2]=s[];
              kk2+=1; 
            }
	  
            double Phiprev[Nequation*Nsubpoint];
            constPhi(Nsubpoint,Nequation,Nboundary0,1,Phiprev,meshprev,solprev,0,param);
            double normphiprev=0.;
            for (int pp=0;pp<Nequation*Nsubpoint;pp++) {
              normphiprev+=sq(Phiprev[pp]);
            }
            normphiprev=sqrt(normphiprev);
	  
            if ((normphiinit>=normphipart)&&(normphiprev>=normphipart)) {
              for (int kk1=0;kk1<Nsubpoint;kk1++) {
                mesh[kk1]=meshpart[kk1];
              }
              for (int kk2=0;kk2<Nequation*Nsubpoint;kk2++) {
                sol[kk2]=solpart[kk2];
              }
              initial=false;
            }
            else if ((normphipart>=normphiinit)&&(normphiprev>=normphiinit)) {
              for (int kk1=0;kk1<Nsubpoint;kk1++) {
                mesh[kk1]=meshinit[kk1];
              }
              for (int kk2=0;kk2<Nequation*Nsubpoint;kk2++) {
                sol[kk2]=solinit[kk2];
              }
              initial=true;
            }
            else if ((normphiinit>=normphiprev)&&(normphipart>=normphiprev)) {
              for (int kk1=0;kk1<Nsubpoint;kk1++) {
                mesh[kk1]=meshprev[kk1];
              }
              for (int kk2=0;kk2<Nequation*Nsubpoint;kk2++) {
                sol[kk2]=solprev[kk2];
              }
              initial=true;
            }  
          }
statssolverpart=ODE_solver(sol,mesh,Nsubpoint,Nequation,Nboundary0,initial ? 1:10,param);
          if (statssolverpart.norm2>1e-5) {  //if we didn't converge we seek a neighbor which as converged to use its solution
            fprintf(stderr,"echec : %g,%g,%g,%d\n",x,y,statssolverpart.norm2,statssolverpart.niter);
            fflush(stderr);
            findneighbor=false;
            foreach_neighbor(1) {
              if (updatedbvp[]>=0.5) {
                findneighbor=true;
                int kk1=0;
                for (scalar s in smesh) {
                  mesh[kk1]=s[];
                  kk1+=1;
                }
                int kk2=0;
                for (scalar s in ssol) {
                  sol[kk2]=s[];
                  kk2+=1; 
                }
              }
            }
            if (!findneighbor) {
              foreach_neighbor() {
                if (updatedbvp[]>=0.5) {
                  findneighbor=true;
                  int kk1=0;
                  for (scalar s in smesh) {
                    mesh[kk1]=s[];
                    kk1+=1;
                  }
                  int kk2=0;
                  for (scalar s in ssol) {
                    sol[kk2]=s[];
                    kk2+=1; 
                  }
		
                }
              }
            }
            if (!findneighbor) {
              fprintf(stderr,"no neighbor : %g,%g\n",x,y);
              fflush(stderr);
              updatedbvp[]=0.;
            }
          }
          else {
            updatedbvp[]=1.;
          }
          if (!initial) {
            initial=true;
          }
          if ((statssolverpart.norm2<1e-5)||(findneighbor)) {
            int kk3=0;
            for (scalar s in smesh) {
              s[]=mesh[kk3];
              kk3+=1;
            }
            int kk4=0;
            for (scalar s in ssol) {
              s[]=sol[kk4];
              kk4+=1;
            }
          }
	
          //we take the information about the solver
          statssolver.npoint+=1.;
          if (statssolverpart.niter<statssolver.nitermin) {
            statssolver.nitermin=statssolverpart.niter;
          }
          if (statssolverpart.niter>statssolver.nitermax) {
            statssolver.nitermax=statssolverpart.niter;
          }
          statssolver.niteravg+=statssolverpart.niter;
          if (statssolverpart.norm2>statssolver.norm2max) {
            statssolver.norm2max=statssolverpart.norm2;
          }
          statssolver.norm2avg+=statssolverpart.norm2;
	
          dpprec=param.dp;
          advprec=param.adv;

        }
        else {
          findneighbor=false;
          statssolverpart.norm2=1.;
        }
	
        //quantity that are needed for the wall function
        bool correctduF=false;
        if ((i>1000)&&((statssolverpart.norm2<1e-5)||(findneighbor))) {
          double mutotF=molvis[]+1.*sol[Nequation*(Nsubpoint-1)+1]*fv1(sol[Nequation*(Nsubpoint-1)+1]/molvis[]);//mu[]+rho[]*sol[Nequation*(Nsubpoint-1)+1]*fv1(sol[Nequation*(Nsubpoint-1)+1]/(mu[]/rho[]));
          if (i<2000) {
            muhat.duF[]=(2000.-i)/1000.*sq(utau)/molvis[]*fwallSAprime(d_IP*utau/molvis[])+(i-1000.)/1000.*sol[Nequation*(Nsubpoint-1)+2]/mutotF;
            correctduF=true;
          }
          else {
            correctduF=true;
            muhat.duF[]=sol[Nequation*(Nsubpoint-1)+2]/mutotF;
          }
	
          if (i<2000) {//(correctduF)&&) {
            muhat.utau[]=(2000.-i)/1000.*utau+(i-1000.)/1000.*sqrt(sol[2])/1.;//rho[];
          }
          else {
            muhat.utau[]=sqrt(sol[2]/1.);//utau;
          }
        }
        else {
          muhat.utau[]=utau;
          muhat.duF[]=sq(utau)/molvis[]*fwallSAprime(d_IP*utau/molvis[]);
        }
      
#else
      
        muhat.utau[]=utau;
        muhat.duF[]=sq(utau)/molvis[]*fwallSAprime(d_IP*utau/molvis[]);
      
#endif

#if BVP4
        if ((i>1000)&&((statssolverpart.norm2<1e-5)||(findneighbor))) {
          muhat.nuhatIP[]=muhat_IP/rho_IP;
          if (i<2000) {
            muhat.dnuhatF[]=(2000.-i)/1000.*kappa*utau+(i-1000.)/1000.*sol[Nequation*(Nsubpoint-1)+3]/(molvis[]+muhat_IP/rho_IP);
            viscSAdiss[]=((2000.-i)/1000.*(molvis[]+kappa*utau*d_IP)*kappa*utau+(i-1000.)/1000.*sol[Nequation*(Nsubpoint-1)+3])/muhat.dnuhatF[]-molvis[];//muhat_IP/rho_IP;
            muhat.viscSAdiss[]=viscSAdiss[];//muhat_IP/rho_IP;
          }
          else {
            muhat.dnuhatF[]=sol[Nequation*(Nsubpoint-1)+3]/(molvis[]+muhat_IP/rho_IP);
            viscSAdiss[]=muhat_IP/rho_IP;
            muhat.viscSAdiss[]=muhat_IP/rho_IP;
          }
        }
        else {
          muhat.nuhatIP[]=kappa*utau*d_IP;
          muhat.dnuhatF[]=kappa*utau;
          viscSAdiss[]=kappa*utau*d_IP;//muhat_IP/rho_IP;
          muhat.viscSAdiss[]=kappa*utau*d_IP;//muhat_IP/rho_IP;
        }
 #else
        muhat.nuhatIP[]=kappa*utau*d_IP;
        muhat.dnuhatF[]=kappa*utau;
        viscSAdiss[]=kappa*utau*d_IP;//muhat_IP/rho_IP;
        muhat.viscSAdiss[]=kappa*utau*d_IP;//muhat_IP/rho_IP;
 #endif      
      
      
      
#if LINEARISATION
      
        double ut=ut_IP+muhat.duF[]*(dstar-d_IP);
        double un=-un_IP*sq(dstar/d_IP)+2.*un_IP*dstar/d_IP;
      
        double cplus=0.;
#if BVP4
        if ((i>1000)&&((statssolverpart.norm2<1e-5)||(findneighbor))) {//&&(correctduF)) {
          if ((correctduF)) {
            if (i<2000) {
              double mutotF=molvis[]+1.*sol[Nequation*(Nsubpoint-1)+1]*fv1(sol[Nequation*(Nsubpoint-1)+1]/molvis[]);
              cplus=((2000.-i)/1000.*sq(utau)+(i-1000.)/1000.*sol[Nequation*(Nsubpoint-1)+2])/muhat.duF[]-molvis[];
            }
            else {
              cplus=sol[Nequation*(Nsubpoint-1)+2]/muhat.duF[]/1.-molvis[]*1.;
            }
          }
          else {
            cplus=molvis[]*(1./fwallSAprime(d_IP*muhat.utau[]/molvis[])-1.);
          }
        }
        else {
          cplus=molvis[]*(1./fwallSAprime(d_IP*muhat.utau[]/molvis[])-1.);//((sq(muhat.utau[])/molvis[])/muhat.duF[]-1.);
        }
#else
        cplus=molvis[]*(1./fwallSAprime(d_IP*utau/molvis[])-1.);//((sq(muhat.utau[])/molvis[])/muhat.duF[]-1.);
#endif
        double sigmaplus=cplus*pow(cv1*molvis[],3.);
        double Rea=0.;
        double Reb=cplus/4.;
        double Rem=(Rea+Reb)/2.;
        double resm;
        while (fabs(Rea-Reb)/Reb>1e-5) {
          Rem=(Rea+Reb)/2.;
          resm=function_Relambda(Rem,cplus,sigmaplus);
          if (resm>0.) {
            Reb=Rem;
          }
          else {
            Rea=Rem;
          }
        }
        double discriminant=sq(cplus-2.*Rem)*(1.+8.*Rem/(cplus-4.*Rem));
        double viscC=0.;
#if BVP4
        viscC=1./2.*(cplus-2.*Rem+sqrt(discriminant));//sol[Nequation*(Nsubpoint-1)+1];
	    
#else
        viscC=1./2.*(cplus-2.*Rem+sqrt(discriminant));
#endif
      
#else

#if BVP4
        double ut=1.;//interp(muhat.sol,muhat.mesh);  /fixme : need to be done
#else
        double ut=muhat.utau[]*fwall(dstar*muhat.utau[]/molvis[]);
#endif
        double un=un_IP*dstar/d_IP;
        double cplus=1.;
        double viscC=max(0.,kappa*muhat.utau[]*dstar);
      
#endif
        muhat.Cplus[]=cplus;
        muhat.uboundary[]=ut_IP+(-d_IP)*muhat.duF[];
        muhat.uyboundary[]=(u.y[]+u.y[0,-1])/2.;
        muhat.dnuhatF[]=0.;
   
        visc[]=viscC;
#if BVP4
        viscSA[]=muhat.nuhatIP[]+(dstar-d_IP)*muhat.dnuhatF[];
#else
        viscSA[]=kappa*utau*dstar;
#endif
        omega[]=sq(utau)/(molvis[]+1.*kappa*utau*dstar*fv1(kappa*utau*dstar/molvis[]));
        ux[]=un*nx+ut*tx;
        uy[]=un*ny+ut*ty;
        pavg2[]=p_IP;
        pfavg2[]=pf_IP;
        muhat.dunIPdt[]=0.;
        muhat.unIPprevious[]=un_IP;
        muhat.pIP[]=p_IP;
        muhat.pfIP[]=pf_IP;
      
        if (reverse) {  //we take utau<0 to show that we need to reverse the tangential vector
          muhat.utau[]=-muhat.utau[];
        }
      }
    }

#if BVP4
    statssolver.niteravg/=(statssolver.npoint+1e-30);
    statssolver.norm2avg/=(statssolver.npoint+1e-30);
#endif
  

    foreach() {
#if BVP4
      updatedbvp[]=0.;
#endif
      double dstar=distance_to_wall(x,y);
      if ((dstar>-Deltamin)&&(dstar<d_IP+Deltamin)&&(x>leading)&&(x<trailing)) {
        muhat[]=max(visc[],0.);
        if (i<100) {
          muhatSA[]=visc[];
          muhatSAdiss[]=visc[];
        }
        else {
          muhatSA[]=viscSA[];
          muhatSAdiss[]=viscSAdiss[];
        }
        u.x[]=ux[];
        u.y[]=uy[];
        p[]=pavg2[];
        pf[]=pfavg2[];
      }
    }
    foreach_face() {
      double chi=face_value(muhat,0)/face_value(molvis,0);
      mut.x[]=fm.x[]*face_value(molvis,0)*(1.+chi*fv1(chi));
    }

    centered_gradient (p,g);
  }
}

/**
   Boundary conditions need to be change. This event allows to hook up other boundary conditions that need to be modified (for example if a solid is a frontier of the domain) */
event boundary_condition(i++) {
#if EMBED
  if (i>100) {
    muhat[embed]=dirichlet(wall_condition_viscosity(x,y,muhat));
    muhatSA[embed]=dirichlet(wall_condition_viscosity_SA(x,y,muhat));
    muhatSAdiss[embed]=dirichlet(wall_condition_viscosity_SA_diss(x,y,muhat));
    u.n[embed]=dirichlet(wall_condition_velocity(x,y,true,muhat));
    u.t[embed]=dirichlet(wall_condition_velocity(x,y,false,muhat));
  }
#endif
  boundary({muhat,muhatSA});
  foreach_face() {
    double chi=face_value(muhat,0)/face_value(molvis,0);
    mut.x[]=fm.x[]*face_value(molvis,0)*(1.+chi*fv1(chi));
  }
  boundary({u.x,u.y});
  boundary({p,pf});
}

#if BVP4
event end (t=end) {
  free(smesh);
  free(ssol);
}
#endif