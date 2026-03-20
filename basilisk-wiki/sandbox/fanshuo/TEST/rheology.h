/** 
# Properties and functions for non newtonian flows */

/*For Bingham*/
double tauy = 0.0;
double mu = 1;

/*For cohesive Bagnold and dry Bagnold*/
double I0 = 0.3 [0];
double mu0 = 0.1 [0];
double deltamu = 0.26 [0];
double dg = 0.4 [1];
double rho = 1.0;
double tauc = 0.0 ; // cohesive stress, =0 if BAGNOLDDRY


/* For turbulent viscous */
double kappa = 0.4;
double rugo = 0.01;
double ll = 0.0;

/* up to now constante slope, to be changed */
double slope = 0.25;

/* new parameter for test case with different viscosity*/
int bilayer = 0;


double shear(Point point, scalar s, scalar h, int layer, int layercoef){

  double shear=0;
  if (layercoef>0&&layercoef<(nl-1)){ 
    shear = (s[0,0,layer+1]-s[0,0,layer-1])/(h[0,0,layer]+0.5*(h[0,0,layer+1]+h[0,0,layer-1]));
  }
  else if(layercoef==0){
     shear = (s[0,0,layer+1]-s[0,0,layer])/(0.5*(h[0,0,layer]+h[0,0,layer+1]));
  }
  else if(layercoef==nl-1){
    shear = (s[0,0,layer]-s[0,0,layer-1])/(0.5*(h[0,0,layer]+h[0,0,layer-1]));
  }
  else if(layercoef==-1){
    shear = 2.*s[0,0,0]/h[0,0,0];
  }
  return fabs(shear);
  //return sqrt(sq(shear));   // Why?
}

/*
Functions reserved for dry  Bagnold flow
*/
double pressionHydro(Point point, scalar h,int layer){

  double H = 0.;
  double zc = 0.;
    for (int l = 0; l < layer; l++) {
    H+=h[0,0,l];
  }
  zc = H + 0.5*h[0,0,layer];
  return rho*G*cos(slope)*(eta[]-zb[]-zc);
}


double nombreInertie(Point point, scalar s, scalar h, int layer){

  double ans;
  ans = dg*shear(point,s,h,layer,layer)/sqrt(pressionHydro(point,h,layer)/rho);
  return ans;
}

double coeffFrotte(Point point, scalar s, scalar h, int layer){

  double _rapport;
  _rapport = I0/nombreInertie(point, s, h, layer);
  return mu0 + deltamu/(_rapport + 1);

}

double Nuturbulent(Point point, scalar s, scalar h, int layer){

  double nu_eq = 0.0;
  double _y = 0.0;

    for (int l = 0; l < layer; l++) {
    _y+=h[0,0,l];
  }
  _y = _y + 0.5*h[0,0,layer];

  ll = kappa*(_y+rugo)*sqrt(1-(_y));
  nu_eq = shear(point,s,h,layer,layer)*ll*ll;
  return nu_eq;
}



double Nueq(Point point, scalar s, scalar h, int layer){

  double ans=0;
#if BAGNOLDDRY
  ans =  coeffFrotte(point,s,h,layer)*pressionHydro(point,h,layer)/shear(point,s,h,layer,layer);

#elif BAGNOLDCOH
  double term1 = 0.;
  double term2 = 0.;
  term1 = deltamu*dg*pow(rho*pressionHydro(point,h,layer),0.5)/I0;
  term2 = (tauc+mu0*pressionHydro(point,h,layer))/shear(point,s,h,layer,layer);
  
  ans = term1 + term2;
  
#elif BINGHAM
    ans = mu + tauy/shear(point,s,h,layer,layer);
  
#elif TURBULENT
  ans = Nuturbulent(point,  s, h, layer);
#else 
  ans = D;
#endif

  return ans;
}

//with bilayer
// double Nueq(Point point, scalar s, scalar h, int layer, double D, int bilayer){
//   double ans=0;
//   if (bilayer == 0) ans = D;
//   else if (bilayer == 1) ans = coeffFrotte(point,s,h,layer)*pressionHydro(point,h,layer)/(shear(point,s,h,layer,layer)); // other regularization
//   else  ans = mu + tauy/(shear(point,s,h,layer,layer)+dry);
//   return ans;
// }

/** Régularisation for $nu_{eq}$ */
void regularization(Point point, scalar s, scalar h, double nueq[], double D)
{
  int nlc = nl, l;

  for(l=0 ; l<nl;l++){ // Find the limiting layer shear
    if (shear(point,s,h,l,l) <=1e-3) {
      nlc = l-1;
    break;
    }
  }

  for( l=0 ; l<nl;l++){
    if ( l<=nlc ) nueq[l] = Nueq(point,s,h,l);
    else nueq[l] = nueq[l-1];
   // nueq[l] = Nueq(point,s,h,l);
  }

  // for(l=0 ; l<nl;l++){ // Potential former bilayer
  //   if ( l<=(nl/2.) )  nueq[l] = Nueq(point,s,h,l, D, bilayer);
  //   else nueq[l] = Nueq(point,s,h,l, D, bilayer+2);
  // }

// Nueq with bilayer need to be tested and compared to analytical solutions
  // for( l=0 ; l<nl;l++){ // Apply regularization
  //   if ( l<=nlc && l <=nl/2 ) nueq[l] = Nueq( point,s,h,l, D, bilayer+1 );
  //   else if ( l<=nlc && l > nl/2 ) nueq[l] = Nueq( point,s,h,l, D, bilayer  );
  //   else nueq[l] = nueq[l-1];
  // }
}



