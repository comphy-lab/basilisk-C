/**
   Function for the k-epsilon model */



/**
   Definition of the constant for the Spalart-Allmaras wall function */
double Bbar,c1,c2,c3,c4;
double a1,a2,b1,b2;

double cv1;

double kappa;   //von Karman constant


/**
   Number of iteration for the Newton method */
int N_itermax;


event defaults(i=0) {
  kappa=0.41;
  cv1=7.1;
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

double viscosity_profile(double utau,double y, double ut1) {
  double uplus=ut1/utau;
  return kappa*exp(-kappa*Bbar)*(exp(kappa*uplus)-1.-kappa*uplus-sq(kappa*uplus)/2.);
}