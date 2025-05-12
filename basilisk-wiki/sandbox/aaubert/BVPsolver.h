/**
   This file implement a Boundary value problem solver 
   
   It solves problem of the form
   $$
   \frac{\partial w}{\partial y}=f(w,y)
   $$
   with $w(0)=w_0$ and $w(L)=w_L$
   
   The function f needs to be defined in the intial file as rhs, the derivative as drhs, the boundary condition at $y=0$ as boundary0 and the boundary condition at $y=L$ as boundaryF
   */


/**
  Defintion of the coefficient for the different method */

#if EULER
//coefficient for the first order Euler scheme
#define sMIRK 1
double c[sMIRK]={0.};//{1./2.};
double v[sMIRK]={0.};//{1./2.};
double X[sMIRK][sMIRK]={0.};
double b[sMIRK]={1.};
#endif

#if MIRK221L
//coefficient for the MIRK221L scheme
//int s=2;   //number of stage for the MIRK scheme
#define sMIRK 2
double c[sMIRK]={1.,1./3.};
double v[sMIRK]={1.,332./825.};
double X[sMIRK][sMIRK]={{0.,0.},{-19./275.,0.}};
double b[sMIRK]={1./4.,3./4.};
#endif

#if MIRK6
//coefficient for the 6th order MIRK scheme
#define sMIRK 9
double c[sMIRK]={0.,1.,1./4.,3./4.,1./2.,7./16.,1./8.,9./16.,3./8.};
double v[sMIRK]={0.,1.,5./32.,27./32.,1./2.,7./16.,1./8.,9./16.,3./8.};
double X[sMIRK][sMIRK]={{0.,0.,0.,0.,0.,0.,0.,0.,0.},
			{0.,0.,0.,0.,0.,0.,0.,0.,0.},
			{9./64.,-3./64.,0.,0.,0.,0.,0.,0.,0.},
			{3./64.,-9./64.,0.,0.,0.,0.,0.,0.,0.},
			{5./24.,-5./24.,2./3.,-2./3.,0.,0.,0.,0.,0.},
			{1547./32768.,-1225./32768.,749./4096.,-287./2048.,-861./16384.,0.,0.,0.,0.},
			{83./1536.,-13./384.,283./1536.,-167./1536.,-49./512.,0.,0.,0.,0.},
			{1225./32768.,-1547./32768.,287./2048.,-749./4096.,861./16384.,0.,0.,0.,0.},
			{233./3456.,-19./1152.,0.,0.,0.,-5./72.,7./72.,-17./216.,0.}};
double b[sMIRK]={7./90.,7./90.,32./90.,32./90.,12./90.,0.,0.,0.,0.};
#endif

#if RK4
//coefficient for the 4th order Explicit RK scheme
#define sMIRK 4
double c[sMIRK]={0,1./2.,1./2.,1.};
double v[sMIRK]={0.,0.,0.,0.};
double X[sMIRK][sMIRK]={{0.,0.,0.,0.},
			{1./2.,0.,0.,0.},
			{0.,1./2.,0.,0.},
			{0.,0.,1.,0.}};
double b[sMIRK]={1./6.,1./3.,1./3.,1./6.};
#endif


//Parameter for the Newton method
int Niteration; //maximal number of iteration
double tol; //tolerance on the solution


typedef struct {
  int niter;
  double norm2;
} BVPstats;



event defaults(i=0) {
  Niteration=50;
  tol=1e-8;
}



/**
   Computation of the array */

//Allows to build the right hand side between the mesh point $y_k$ and $y_{k+1}$
void constpartialPhi(double *res, double *wk, double *wkp1, double yk, double dy,int Nequation,int Nboundary0, int Nitermax, int niter, Param param) {
  double fYj[Nequation][sMIRK];
  for (int jj=0;jj<sMIRK;jj++) {
    double Yj[Nequation];
    for (int kk=0;kk<Nequation;kk++) {
      Yj[kk]=(1-v[jj])*wk[kk]+v[jj]*wkp1[kk];
      for (int ll=0;ll<jj;ll++) {
	Yj[kk]+=dy*X[jj][ll]*fYj[kk][ll];
      }
    }
    double resf[Nequation];
    rhs(resf,yk+c[jj]*dy,Yj,niter,Nitermax,param);

    for (int kk=0;kk<Nequation;kk++) {
      fYj[kk][jj]=resf[kk];
    }
  }
  for (int kk=0;kk<Nequation;kk++) {
    res[kk]=wkp1[kk]-wk[kk];
    for (int jj=0;jj<sMIRK;jj++) {
      res[kk]-=dy*b[jj]*fYj[kk][jj];
    }
  }
}


// This function construct the array that we want equals to zero
void constPhi(int Nsubpoint, int Nequation, int Nboundary0, int Nitermax, double *phi, double *Y, double *W, int niter, Param param) {

  double w0[Nequation];
  for (int kk=0;kk<Nequation;kk++) {
    w0[kk]=W[kk];
  }
  double Boundary0[Nboundary0];
  boundary0(Boundary0,w0,param);
  for (int kk=0;kk<Nboundary0;kk++) {
    phi[kk]=Boundary0[kk];
  }
  for (int kk=0;kk<Nsubpoint-1;kk++) {
    double wk[Nequation];
    double wkp1[Nequation];
    for (int ll=0;ll<Nequation;ll++) {
      wk[ll]=W[kk*Nequation+ll];
      wkp1[ll]=W[(kk+1)*Nequation+ll];
    }
    double result[Nequation];
    constpartialPhi(result,wk,wkp1,Y[kk],Y[kk+1]-Y[kk],Nequation,Nboundary0,Nitermax,niter,param);
    for (int ll=0;ll<Nequation;ll++) {
      phi[Nboundary0+kk*Nequation+ll]=result[ll];
    }
  }
  double wF[Nequation];
  for (int kk=0;kk<Nequation;kk++) {
    wF[kk]=W[(Nsubpoint-1)*Nequation+kk];
  }
  double BoundaryF[Nequation-Nboundary0];
  boundaryF(BoundaryF,wF,param);
  for (int kk=0;kk<Nequation-Nboundary0;kk++) {
    phi[Nboundary0+(Nsubpoint-1)*Nequation+kk]=BoundaryF[kk];
  }
}

/**
   Computation of the Jacobian */

//This function build the Jacobian between two mesh points $y_k$ and $y_{k+1}$
void constpartialJacobian(int Nequation, int Nitermax, double Lk[Nequation][Nequation], double Rk[Nequation][Nequation], double *wk, double *wkp1, double yk, double dy, int niter, Param param) {
  double dfYjdk[Nequation][Nequation][sMIRK];
  double dfYjdkp1[Nequation][Nequation][sMIRK];
  double fYj[Nequation][sMIRK];
  for (int jj=0;jj<sMIRK;jj++) {
    double Yj[Nequation];
    double dYjdk[Nequation][Nequation];
    double dYjdkp1[Nequation][Nequation];
    for (int pp=0;pp<Nequation;pp++) {
      Yj[pp]=(1-v[jj])*wk[pp]+v[jj]*wkp1[pp];
      for (int ll=0;ll<jj;ll++) {
	Yj[pp]+=dy*X[jj][ll]*fYj[pp][ll];
      }
      for (int qq=0;qq<Nequation;qq++) {
	dYjdk[pp][qq]=(pp==qq) ? 1.-v[jj] : 0.;
	dYjdkp1[pp][qq]=(pp==qq) ? v[jj] : 0.;
	for (int ll=0;ll<jj;ll++) {
	  dYjdk[pp][qq]+=dy*X[jj][ll]*dfYjdk[pp][qq][ll];
	  dYjdkp1[pp][qq]+=dy*X[jj][ll]*dfYjdkp1[pp][qq][ll];
	}
      }
    }

    double dresf[Nequation][Nequation];
    drhs(Nequation,dresf,yk+c[jj]*dy,Yj,niter,Nitermax,param);

    double resf[Nequation];
    rhs(resf,yk+c[jj]*dy,Yj,niter,Nitermax,param);

    for (int pp=0;pp<Nequation;pp++) {
      fYj[pp][jj]=resf[pp];
      for (int qq=0;qq<Nequation;qq++) {
	dfYjdk[pp][qq][jj]=0.;
	dfYjdkp1[pp][qq][jj]=0.;
	for (int kk=0;kk<Nequation;kk++) { //changement produit matriciel
	  dfYjdk[pp][qq][jj]+=dresf[pp][kk]*dYjdk[kk][qq];
	  dfYjdkp1[pp][qq][jj]+=dresf[pp][kk]*dYjdkp1[kk][qq];
	}
      }
    }
  }
  for (int pp=0;pp<Nequation;pp++) {
    for (int qq=0;qq<Nequation;qq++) {
      Lk[pp][qq]=(pp==qq) ? -1. : 0.;
      Rk[pp][qq]=(pp==qq) ? 1. : 0.;
      for (int jj=0;jj<sMIRK;jj++) {
	Lk[pp][qq]-=dy*b[jj]*dfYjdk[pp][qq][jj];
	Rk[pp][qq]-=dy*b[jj]*dfYjdkp1[pp][qq][jj];
      }
    }
  }
}

//This function build the Jacobian for the Newton iterations. This matrix will be almost block diagonal with left diagonal L and right diagonal R. This block will be square matrix of size Nequation.
void constJacobian(int Nsubpoint, int Nequation, int Nboundary0, int Nitermax, double L[Nsubpoint-1][Nequation][Nequation], double R[Nsubpoint-1][Nequation][Nequation], double dbound0[Nboundary0][Nequation], double dboundF[Nequation-Nboundary0][Nequation], double *Phi, double *mesh, double *W, int niter, Param param) {

  //boundary conditions
  double w0[Nequation];
  for (int ll=0; ll<Nequation;ll++) {
    w0[ll]=W[ll];
  }
  dboundary0(Nequation,Nboundary0,dbound0,w0,param);
  
  double wF[Nequation];
  for (int ll=0; ll<Nequation;ll++) {
    wF[ll]=W[(Nsubpoint-1)*Nequation+ll];
  }
  dboundaryF(Nequation,Nboundary0,dboundF,wF,param);

  for (int kk=0;kk<Nsubpoint-1;kk++) {
    double wk[Nequation];
    double wkp1[Nequation];
    for (int ll=0; ll<Nequation;ll++) {
      wk[ll]=W[kk*Nequation+ll];
      wkp1[ll]=W[(kk+1)*Nequation+ll];
    }
    constpartialJacobian(Nequation,Nitermax,L[kk],R[kk],wk,wkp1,mesh[kk],mesh[kk+1]-mesh[kk],niter,param);
  }
}

/**
   Inversion of the Jacobian */

//This function solve the linear system Jacobian*deltaW=Phi by using Gaussian elimination. The special structure of the Jacobian allows some simplifications. We assume without loss of generality that dbound0[0][0] (the upper left coefficient of the derivative of the boundary condition at 0) is not zero (a permutation of the variable can be needed depending on what are the boundary conditions).
void invJacobian(int Nsubpoint, int Nequation, int Nboundary0, double *deltaW, double L[Nsubpoint-1][Nequation][Nequation], double R[Nsubpoint-1][Nequation][Nequation], double dbound0[Nboundary0][Nequation], double dboundF[Nequation-Nboundary0][Nequation], double *Phi,int Perm[Nequation*Nsubpoint]) {

  //zero step : we conditionned the system to avoid numerical error : we divide each lign by the absolute maximum of this lign (we do nothing for the boundary condition are they are already well-conditionned)
  for (int ll=0;ll<Nsubpoint-1;ll++) {
    for (int ii=0;ii<Nequation;ii++) {
      double max=0.;
      for (int jj=0;jj<2*Nequation;jj++) {
	double value=0.;
	if (jj<Nequation) {
	  value=L[ll][ii][jj];
	}
	else {
	  value=R[ll][ii][jj-Nequation];
	}
	if (fabs(value)>fabs(max)) {
	  max=value;
	}
      }
      assert(max);
      if (fabs(max)>1e-5) {
	for (int jj=0;jj<2*Nequation;jj++) {
	  if (jj<Nequation) {
	    L[ll][ii][jj]/=fabs(max);
	  }
	  else {
	    R[ll][ii][jj-Nequation]/=fabs(max);
	  }
	}
	Phi[Nboundary0+ll*Nequation+ii]/=fabs(max);
      }
    }
  }
  

  //first step : triangulation of the dBound0 matrix in the upper left corner
  for (int ii=1;ii<Nboundary0;ii++) {  //nothing to do on the first line because dbound0[0][0] is not null
    for (int kk=ii;kk<Nboundary0;kk++) {
      double pivot=dbound0[kk][ii-1]/dbound0[ii-1][ii-1];
      for (int jj=ii-1;jj<Nequation;jj++) {
	dbound0[kk][jj]-=pivot*dbound0[ii-1][jj];
      }
      Phi[kk]-=pivot*Phi[ii-1];
    }
  }

  //second step : triangulation of the left and right diagonal matrix
  for (int ll=0;ll<Nsubpoint-1;ll++) {
    if (ll==0) {
      for (int kk=0;kk<Nboundary0;kk++) {
	for (int ii=0;ii<Nequation;ii++) {
	  double pivot=L[ll][ii][kk]/dbound0[kk][kk];
	  for (int jj=kk;jj<Nequation;jj++) {
	    L[ll][ii][jj]-=pivot*dbound0[kk][jj];
	  }
	  Phi[Nboundary0+ii]-=pivot*Phi[kk];
	}
      }
    }
    else {
      for (int kk=0;kk<Nboundary0;kk++) {
	for (int ii=0;ii<Nequation;ii++) {
	  double pivot=L[ll][ii][kk]/R[ll-1][Nequation-Nboundary0+kk][kk];
	  for (int jj=kk;jj<Nequation;jj++) {
	    L[ll][ii][jj]-=pivot*R[ll-1][Nequation-Nboundary0+kk][jj];
	  }
	  Phi[Nboundary0+ll*Nequation+ii]-=pivot*Phi[Nboundary0+ll*Nequation-Nboundary0+kk];
	}
      }
    }
    for (int kk=0;kk<Nequation-Nboundary0;kk++) {
      //we seek the best pivot
      double maxp=L[ll][kk][Nboundary0+kk];
      int intm=kk;
      for (int pp=kk+1;pp<Nequation;pp++) {
	if (fabs(L[ll][pp][Nboundary0+kk])>fabs(maxp)) {
	  maxp=L[ll][pp][Nboundary0+kk];
	  intm=pp;
	}
      }
      if (intm!=kk) {
	for (int jj=Nboundary0+kk;jj<Nequation;jj++) {
	  double value=L[ll][kk][jj];
	  L[ll][kk][jj]=L[ll][intm][jj];
	  L[ll][intm][jj]=value;
	}
	for (int jj=0;jj<Nequation;jj++) {
	  double value=R[ll][kk][jj];
	  R[ll][kk][jj]=R[ll][intm][jj];
	  R[ll][intm][jj]=value;
	}
	double valuep=Phi[Nboundary0+ll*Nequation+kk];
	Phi[Nboundary0+ll*Nequation+kk]=Phi[Nboundary0+ll*Nequation+intm];
	Phi[Nboundary0+ll*Nequation+intm]=valuep;
        int intp=Perm[Nboundary0+ll*Nequation+kk];  //we keep track of the permutation to recover the solution
	Perm[Nboundary0+ll*Nequation+kk]=Perm[Nboundary0+ll*Nequation+intm];
	Perm[Nboundary0+ll*Nequation+intm]=intp;
      }
      for (int ii=kk+1;ii<Nequation;ii++) {
	assert(L[ll][kk][Nboundary0+kk]);
	double pivot=L[ll][ii][Nboundary0+kk]/L[ll][kk][Nboundary0+kk];
	for (int jj=Nboundary0+kk;jj<Nequation;jj++) {
	  L[ll][ii][jj]-=pivot*L[ll][kk][jj];
	}
	for (int jj=0;jj<Nequation;jj++) {
	  R[ll][ii][jj]-=pivot*R[ll][kk][jj];
	}
	Phi[Nboundary0+ll*Nequation+ii]-=pivot*Phi[Nboundary0+ll*Nequation+kk];
      }
    }
    for (int kk=Nequation-Nboundary0;kk<Nequation;kk++) {
      //we seek the best pivot
      double maxp=R[ll][kk][kk-Nequation+Nboundary0];
      int intm=kk;
      for (int pp=kk+1;pp<Nequation;pp++) {
	if (fabs(R[ll][pp][kk-Nequation+Nboundary0])>fabs(maxp)) {
	  maxp=R[ll][pp][kk-Nequation+Nboundary0];
	  intm=pp;
	}
      }
      if (intm!=kk) {
	for (int jj=kk-Nequation+Nboundary0;jj<Nequation;jj++) {
	  double value=R[ll][kk][jj];
	  R[ll][kk][jj]=R[ll][intm][jj];
	  R[ll][intm][jj]=value;
	}
	double valuep=Phi[Nboundary0+ll*Nequation+kk];
	Phi[Nboundary0+ll*Nequation+kk]=Phi[Nboundary0+ll*Nequation+intm];
	Phi[Nboundary0+ll*Nequation+intm]=valuep;
	int intp=Perm[Nboundary0+ll*Nequation+kk];  //we keep track of the permutation to recover the solution
	Perm[Nboundary0+ll*Nequation+kk]=Perm[Nboundary0+ll*Nequation+intm];
	Perm[Nboundary0+ll*Nequation+intm]=intp;
      }
      for (int ii=kk+1;ii<Nequation;ii++) {
	double pivot=R[ll][ii][kk-Nequation+Nboundary0]/R[ll][kk][kk-Nequation+Nboundary0];
	for (int jj=kk-Nequation+Nboundary0;jj<Nequation;jj++) {
	  R[ll][ii][jj]-=pivot*R[ll][kk][jj];
	}
	Phi[Nboundary0+ll*Nequation+ii]-=pivot*Phi[Nboundary0+ll*Nequation+kk];
      }
    }
  }

  //third step : triangulation of the lower right block
  for (int kk=0;kk<Nboundary0;kk++) {
    for (int ii=0;ii<Nequation-Nboundary0;ii++) {
      double pivot=dboundF[ii][kk]/R[Nsubpoint-2][Nequation-Nboundary0+kk][kk];
      for (int jj=kk;jj<Nequation;jj++) {
	dboundF[ii][jj]-=pivot*R[Nsubpoint-2][Nequation-Nboundary0+kk][jj];
      }
      Phi[Nboundary0+(Nsubpoint-1)*Nequation+ii]-=pivot*Phi[Nboundary0+(Nsubpoint-1)*Nequation-Nboundary0+kk];
    }
  }
  for (int kk=0;kk<Nequation-Nboundary0;kk++) {
    for (int ii=kk+1;ii<Nequation-Nboundary0;ii++) {
      double pivot=dboundF[ii][Nboundary0+kk]/dboundF[kk][Nboundary0+kk];
      for (int jj=Nboundary0+kk;jj<Nequation;jj++) {
	dboundF[ii][jj]-=pivot*dboundF[kk][jj];
      }
      Phi[Nboundary0+(Nsubpoint-1)*Nequation+ii]-=pivot*Phi[Nboundary0+(Nsubpoint-1)*Nequation+kk];
    }
  }

  //fourth step : the Jacobian is now upper triangular so we can solve backward to obtain the solution
  for (int ii=Nequation-Nboundary0-1;ii>=0;ii--) {
    deltaW[Nboundary0+(Nsubpoint-1)*Nequation+ii]=Phi[Nboundary0+(Nsubpoint-1)*Nequation+ii];
    for (int jj=Nequation-1;jj>=ii+Nboundary0+1;jj--) {
      deltaW[Nboundary0+(Nsubpoint-1)*Nequation+ii]-=deltaW[Nboundary0+(Nsubpoint-1)*Nequation-Nboundary0+jj]*dboundF[ii][jj];
    }
    deltaW[Nboundary0+(Nsubpoint-1)*Nequation+ii]/=dboundF[ii][Nboundary0+ii];
  }

  for (int ll=Nsubpoint-2;ll>=0;ll--) {
    for (int ii=Nequation-1;ii>=0;ii--) {
      deltaW[Nboundary0+Nequation*(ll)+ii]=Phi[Nboundary0+Nequation*(ll)+ii];
      for (int jj=Nequation-1;jj>=ii-Nequation+Nboundary0+1;jj--) {
	if (jj>=0) {
	  deltaW[Nboundary0+Nequation*ll+ii]-=deltaW[Nboundary0+Nequation*(ll)+jj+Nequation-Nboundary0]*R[ll][ii][jj];
	}
	else {
	  deltaW[Nboundary0+Nequation*ll+ii]-=deltaW[Nboundary0+Nequation*(ll)+jj+Nequation-Nboundary0]*L[ll][ii][jj+Nequation];
	}
      }
      if (ii-Nequation+Nboundary0>=0) {
	deltaW[Nboundary0+Nequation*ll+ii]/=R[ll][ii][ii-Nequation+Nboundary0];
      }
      else {
	deltaW[Nboundary0+Nequation*ll+ii]/=L[ll][ii][ii+Nboundary0];
      }
    }
  }

  for (int ii=Nboundary0-1;ii>=0;ii--) {
    deltaW[ii]=Phi[ii];
    for (int jj=Nequation-1;jj>=ii+1;jj--) {
      deltaW[ii]-=deltaW[jj]*dbound0[ii][jj];
    }
    deltaW[ii]/=dbound0[ii][ii];
  }
}



BVPstats ODE_solver(double *sol, double *mesh, int Nsubpoint, int Nequation,int Nboundary0, int Nitermax, Param param) {

  double Phi[Nequation*Nsubpoint];   //function that we want equals to zero
  
  //right and left diagonal for the Jacobian
  double R[Nsubpoint-1][Nequation][Nequation];
  double L[Nsubpoint-1][Nequation][Nequation];

  //jacobian of the boundary condition
  double dbound0[Nboundary0][Nequation];
  double dboundF[Nequation-Nboundary0][Nequation];

  double norm2=1.;
  double normphi=1.;
  int compteurt=0;
  for (int nn=0;nn<Nitermax;nn++) {
    norm2=1.;
    normphi=1.;
    int compteur=0;
    bool change=false;
    double normphiinit=0.;
    while ((compteur<Niteration)&&((norm2>tol)||(normphi>min(tol,max(1e-10,1e-5*normphiinit))))&&(normphi<1e10)&&(norm2<1e10)) {
      constPhi(Nsubpoint,Nequation,Nboundary0,Nitermax,Phi,mesh,sol,nn,param);
      constJacobian(Nsubpoint,Nequation,Nboundary0,Nitermax,L,R,dbound0,dboundF,Phi,mesh,sol,nn,param);
      
      double normphiprec=normphi;
      normphi=0.;
      for (int pp=0;pp<Nequation*Nsubpoint;pp++) {
	normphi+=sq(Phi[pp]);
      }
      normphi=sqrt(normphi);
      if (compteur==0) {
	normphiinit=normphi;
      }

      double deltaW[Nequation*Nsubpoint];
      int Perm[Nequation*Nsubpoint];
      for (int pp=0;pp<Nequation*Nsubpoint;pp++) {
	Perm[pp]=pp;
      }
      invJacobian(Nsubpoint,Nequation,Nboundary0,deltaW,L,R,dbound0,dboundF,Phi,Perm);
      
      norm2=0.;
      for (int pp=0;pp<Nequation*Nsubpoint;pp++) {
	norm2+=sq(deltaW[pp]/(sol[pp]+1e-20));
      }
      norm2=sqrt(norm2);
      double alpha=1.;
      if ((compteur>2)&&(!change)&&(normphi>0.2)&&(Nitermax<2)) {
	Nitermax=1;
	change=true;
      }
      else {
	change=false;
	double soltemp[Nequation*Nsubpoint];
	for (int pp=0;pp<Nequation*Nsubpoint;pp++) {
	  soltemp[pp]=sol[pp]-deltaW[pp];
	}
	if (compteur>0) {
	  constPhi(Nsubpoint,Nequation,Nboundary0,Nitermax,Phi,mesh,soltemp,nn,param);
	  double normphitemp=0.;
	  for (int pp=0;pp<Nequation*Nsubpoint;pp++) {
	    normphitemp+=sq(Phi[pp]);
	  }
	  normphitemp=sqrt(normphitemp);
	  int compteurtemp=0;
	  double minnormphitemp=normphitemp;
	  while ((compteurtemp<100)&&(normphitemp>normphi)) {
	    normphitemp=0.;
	    for (int pp=0;pp<Nequation*Nsubpoint;pp++) {
	      soltemp[pp]=sol[pp]+0.1*noise()*sol[pp];
	    }
	    constPhi(Nsubpoint,Nequation,Nboundary0,Nitermax,Phi,mesh,soltemp,nn,param);
	    for (int pp=0;pp<Nequation*Nsubpoint;pp++) {
	      normphitemp+=sq(Phi[pp]);
	    }
	    normphitemp=sqrt(normphitemp);
	    if (normphitemp<minnormphitemp) {
	      minnormphitemp=normphitemp;
	    }
	    compteurtemp+=1;
	  }
	  if (compteurtemp>0) {
	  }
	  if ((compteurtemp>0)&&(compteurtemp<100)) {
	    for (int pp=0;pp<Nequation*Nsubpoint;pp++) {
	      sol[pp]=soltemp[pp];
	    }
	    normphi=normphitemp;
	  }
	  else {
	    for (int pp=0;pp<Nequation*Nsubpoint;pp++) {
	      sol[pp]-=deltaW[pp];
	    }
	  }
	}
	else {
	  for (int pp=0;pp<Nequation*Nsubpoint;pp++) {
	    sol[pp]-=deltaW[pp];
	  }
	}
      }
      compteur+=1;

    }
    compteurt+=compteur;
  }

  
  BVPstats st;
  st.niter=compteurt;
  st.norm2=norm2;

  return st;
}