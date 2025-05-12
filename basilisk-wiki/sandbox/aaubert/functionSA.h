/**
   Function for the Spalart-Allmaras model */


/**
   Definition of all the constant for the model */
double cb1,sigma,cb2,kappa,cw2,cw3,cv1,ct3,ct4,cw1,c2,c3;

/**
   Definition of the constant for the Spalart-Allmaras wall function */
double Bbar,c1,c2,c3,c4;
double a1,a2,b1,b2;

/**
   Number of iteration for the Newton method */
int N_itermax;


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
  Bbar=5.0333908790505579;
  a1=8.148221580024245;
  a2=-6.9287093849022945;
  b1=7.4600876082527945;
  b2=7.468145790401841;
  c1=2.5496773539754747;
  c2=1.3301651588535228;
  c3=3.599459109332379;
  c4=3.6397531868684494;
  N_itermax=50;
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

double fwallSA(double yplus) {
  return Bbar+c1*log(pow(yplus+a1,2.)+pow(b1,2.))-c2*log(pow(yplus+a2,2.)+pow(b2,2.))-c3*atan2(b1,yplus+a1)-c4*atan2(b2,yplus+a2);
}

double fwallSAprime(double yplus) {
  return (pow(cv1,3.)+pow(kappa*yplus,3.))/(pow(cv1,3.)+pow(kappa*yplus,3.)*(1.+kappa*yplus));//2*c1*(yplus+a1)/(pow(yplus+a1,2.)+pow(b1,2.))-2*c2*(yplus+a2)/(pow(yplus+a2,2.)+pow(b2,2.))+c3*b1/(pow(yplus+a1,2.)+pow(b1,2.))+c4*b2/(pow(yplus+a2,2.)+pow(b2,2.));
}


double gwallSA(double utau,double y1,double ut1, double molvis) {
  double yplus=y1*utau/molvis;
  return utau*fwallSA(yplus)-ut1;
}

double gwallSAprime(double utau, double y1, double ut1, double molvis) {
  double yplus=y1*utau/molvis;
  return fwallSA(yplus)+y1*utau/molvis*fwallSAprime(yplus);
}

//This function does the Newton method to find $u_\tau$ given the tangential velocity and the distance to the wall for the image point and an intial guess for $u_\tau$ using the Spalart-Allmaras wall function.
double newton_utau(double utau, double ut1, double d, double molvis) {
  double utauprec;
  double err=1.;
  int n_iter=0;
  while ((n_iter<N_itermax)&&(err>1e-5)) {
    utauprec=utau;
    utau=utau-gwallSA(utau,d,ut1,molvis)/gwallSAprime(utau,d,ut1,molvis);
    err=fabs(utau-utauprec)/(utau+1e-10);
    n_iter+=1;
  }
  return utau;
}

/**
   Used for the computation of $\tilde{\nu}$ when $\nu_t$ is known */
double function_Relambda(double Rel,double Cplus,double sigmaplus) {
  return 8.*pow(Rel,3.)*pow(Cplus-2.*Rel,3.)-sigmaplus*sq(Cplus-4.*Rel);
}