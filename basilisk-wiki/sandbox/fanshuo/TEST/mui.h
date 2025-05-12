/** 
 # Properties and functions for non newtonian flows */

/*For Bingham*/
double tauy = 0.0;
double mu = 0.1;

/*For cohesive Bagnold*/
double I0 = 0.3 [0];
double mu0 = 0.1 [0];
double deltamu = 0.26 [0];
double dg = 0.4 [1];
double rho = 1.0;
double tauc = 0.0 ;

double slope = 0.25;


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
   
  return shear;
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


double Nueq(Point point, scalar s, scalar h, int layer){

  double ans=0;
#if BAGNOLDDRY
  	ans =  coeffFrotte(point,s,h,layer)*pressionHydro(point,h,layer)/shear(point,s,h,layer,layer);
#elif BINGHAM
  	ans = mu + tauy/shear(point,s,h,layer,layer);
#else 
  	ans = D;
#endif

  return ans;
}


/** Régularisation for nu_{eq} */
void regularization(Point point, scalar s, scalar h, double nueq[])
{
	int nlc = nl, l;

	for(l=0 ; l<nl;l++){
 		if (shear(point,s,h,l,l) <=1e-3) {
 			nlc = l-1;
 		break;
		}
	}

  for( l=0 ; l<nl;l++){
  	if ( l<=nlc ) nueq[l] = Nueq(point,s,h,l);
	  else nueq[l] = nueq[l-1];
	}
}

//to do: add other properties for non newtonian flows









