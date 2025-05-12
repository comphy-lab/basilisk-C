//Some functions are defined 

#include "./hydroF.h"
#include "./diffusionF.h"



/*
Functions reserved for Poiseuille flow
*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - */
double Up(Point point,double HR)
{
  double zc = 0.;
  for (int l = - point.l; l < nl - point.l; l++) {
    if (l < 0)
      zc += h[0,0,l];
  }
  zc += h[]/2.;
  return HR*G*sin(slope)/nu*zc*(1-0.5*zc);
}

double Up2(double z,double HR)
{
    return HR*G*sin(slope)/nu*z*(1-0.5*z);
}

double Sp(double z, double HR)
{
  
  return HR*G*sin(slope)/nu*(1-z);

}
/*- - - - - - - - - - - - - - - - - - - - - - - - - */



/*
Functions reserved for Bingham flow
*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - */
double Ubgm(Point point, double HR, double tauy, double mu)
{
  double Yc = 0.0;
  double zc = 0.0;
  double Umax;
  for (int l = - point.l; l < nl - point.l; l++) {
    if (l < 0)
      zc += h[0,0,l];
  }
  zc += h[]/2.;

  Yc = HR - tauy/(G*sin(slope));
  if (Yc>0) Umax = pow(tauy-HR*G*sin(slope),2)/(2*mu*G*sin(slope));
  else {
    printf("No flow, change the rheology parameters.\n");
	return 1;
  }


  if(zc>=Yc) return Umax;
  else return Umax*(1-pow(zc/Yc-1,2));

}

double Ubgm2(double z,double HR, double tauy, double mu)
{


  double Yc = 0.0;
  double Umax;
  
  Yc = HR - tauy/(G*sin(slope));
  if (Yc>0) Umax = pow(tauy-HR*G*sin(slope),2)/(2*mu*G*sin(slope));
  else {
    printf("No flow, change the rheology parameters.\n");
	return 1;
  }

  if(z>=Yc) return Umax;
  else return Umax*(1-pow(z/Yc-1,2));

}

double Sbgm(double z,double HR, double tauy, double mu)
{


  double Yc = 0.0;
  double Umax;
  
  Yc = HR - tauy/(G*sin(slope));
  if (Yc>0) Umax = pow(tauy-HR*G*sin(slope),2)/(2*mu*G*sin(slope));
  else {
    printf("No flow, change the rheology parameters.\n");
	return 1;
  }

  if(z>=Yc) return 0;
  else return -2*Umax/Yc*(z/Yc-1);

}

double muBingham(Point point, int nl, scalar s, scalar h,int layer, double tauy, double mu){

  double shear;

  if (layer>0&&layer<(nl-1)){ 
	shear = (s[0,0,layer+1]-s[0,0,layer])/(0.5*(h[0,0,layer]+h[0,0,layer+1]));}
	else if(layer==0){
	  shear = (s[0,0,layer]-0)/(0.5*(h[0,0,layer]));
	}
	else if(layer==nl-1){
	  shear = (s[0,0,layer]-s[0,0,layer-1])/(0.5*(h[0,0,layer]+h[0,0,layer-1]));
	}
  return mu+tauy/shear;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - */



/*
Functions reserved for cohesive Bagnold flow 
*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - */


double Kb(double I0, double mu0, double deltamu, double dg){

  return I0*(tan(slope)-mu0)*pow(G*cos(slope),0.5)/(deltamu*dg);

}

double tc(double I0, double tauc, double deltamu, double rho, double dg)
{
  double ans;
  ans = I0*tauc*pow(G/cos(slope),0.5)/(deltamu*rho*G*dg);
  return ans;
}

double Ubc(Point point,double HR, double I0, double tauc, double deltamu, double rho, double dg)
{
  double yc;
  yc = HR-tc(I0,tauc,deltamu,rho,dg)/Kb(I0,mu0,deltamu,dg);
  if (yc<=0){
    printf("No flow, change the rheology parameters.\n");
	return 1;
  }

  double zc = 0.;
  for (int l = - point.l; l < nl - point.l; l++) {
    if (l < 0)
      zc += h[0,0,l];
  }
  zc += h[]/2.;

  if(zc<yc) {

    return 2./3.*Kb(I0,mu0,deltamu,dg)*(1.-pow(1-zc/HR,1.5))-2.*tc(I0,tauc,deltamu,rho,dg)*
      pow(HR,0.5)*(1-pow(1-zc/HR,0.5));
  }
  else{

    return 2./3.*Kb(I0,mu0,deltamu,dg)*pow(HR,1.5)-2.*tc(I0,tauc,deltamu,rho,dg)*pow(HR,0.5)+
	  (4./3.)*pow(tc(I0,tauc,deltamu,rho,dg),1.5)/pow(Kb(I0,mu0,deltamu,dg),0.5);
	 
  }
}


double Ubc2(double zc,double HR, double I0, double tauc, double deltamu, double rho, double dg)
{
  double yc;
  yc = HR-tc(I0,tauc,deltamu,rho,dg)/Kb(I0,mu0,deltamu,dg);
  if(zc<yc) {

    return 2./3.*Kb(I0,mu0,deltamu,dg)*(1.-pow(1-zc/HR,1.5))-2.*tc(I0,tauc,deltamu,rho,dg)*
      pow(HR,0.5)*(1-pow(1-zc/HR,0.5));
  }
  else{

    return 2./3.*Kb(I0,mu0,deltamu,dg)*pow(HR,1.5)-2.*tc(I0,tauc,deltamu,rho,dg)*pow(HR,0.5)+
	  (4./3.)*pow(tc(I0,tauc,deltamu,rho,dg),1.5)/pow(Kb(I0,mu0,deltamu,dg),0.5);
	 
  }
}


double Sbc(double zc,double HR, double I0, double tauc, double deltamu, double rho, double dg)
{
  double yc;
  yc = HR-tc(I0,tauc,deltamu,rho,dg)/Kb(I0,mu0,deltamu,dg);
  if(zc<yc) {

    return Kb(I0,mu0,deltamu,dg)*pow(1-zc/HR,0.5)/HR-tc(I0,tauc,deltamu,rho,dg)*
      pow(HR,-0.5)*pow(1-zc/HR,-0.5);
  }
  else{

    return 0;
	 
  }
}



double pressionHydro(Point point, scalar h,int layer, double rho){

  double H = 0.;
  double zc = 0.;
  for (int l = 0; l < layer; l++){
    H+=h[0,0,l];
  }
  zc = H + 0.5*h[0,0,layer];
  return rho*G*cos(slope)*(eta[]-zb[]-zc);
  }


double muBagnoldC(Point point, double theta, int nl, scalar s, scalar h,int layer,double deltamu, double I0, double mu0, double dg, double rho, double tauc){
  double shear;
  double term1;
  double term2;
  if (layer>0&&layer<(nl-1)){ 
    shear = (s[0,0,layer+1]-s[0,0,layer])/(0.5*(h[0,0,layer]+h[0,0,layer+1]));}
  else if(layer==0){
    shear = (s[0,0,layer]-0)/(0.5*(h[0,0,layer]));
  }
  else if(layer==nl-1){
    shear = (s[0,0,layer]-s[0,0,layer-1])/(0.5*(h[0,0,layer]+h[0,0,layer-1]));
  }

  double H = 0.;
  double zc = 0.;
  for (int l = 0; l < layer; l++) {
    H+=h[0,0,l];
  }
  term1 = deltamu*dg*pow(rho*pressionHydro(point,h,layer,rho),0.5)/I0;
  term2 = (tauc+mu0*pressionHydro(point,h,layer,rho))/shear;
  return term1+term2;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - */
/*
Functions reserved for dry Bagnold flow 
*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - */


double nombreInertie(Point point, int nl, scalar s, scalar h, int layer, double dg, double rho){

  double shear;
  double ans;
  if (layer>0&&layer<(nl-1)){ 
      shear = (s[0,0,layer+1]-s[0,0,layer])/(0.5*(h[0,0,layer]+h[0,0,layer+1]));}
  else if(layer==0){
    shear = (s[0,0,layer]-0)/(0.5*(h[0,0,layer]));
  }
  else if(layer==nl-1){
    shear = (s[0,0,layer]-s[0,0,layer-1])/(0.5*(h[0,0,layer]+h[0,0,layer-1]));
  }
	
  ans = dg*shear/sqrt(pressionHydro(point, h, layer, rho)/rho);
  return ans;

}

double coeffFrotte(Point point,int nl, scalar s, scalar h, int layer, double dg, double rho, double mu0, double deltamu, double I0){

  double _rapport;
  _rapport = I0/nombreInertie(point, nl, s, h, layer, dg, rho);
  return mu0 + deltamu/(_rapport + 1);

}





double muBagnold(Point point, double theta, int nl, scalar s, scalar h,int layer,double deltamu, double I0, double mu0, double dg, double rho){
  double shear;
  double ans;
  if (layer>0&&layer<(nl-1)){ 
    shear = (s[0,0,layer+1]-s[0,0,layer])/(0.5*(h[0,0,layer]+h[0,0,layer+1]));}
  else if(layer==0){
    shear = (s[0,0,layer]-0)/(0.5*(h[0,0,layer]));
  }
  else if(layer==nl-1){
    shear = (s[0,0,layer]-s[0,0,layer-1])/(0.5*(h[0,0,layer]+h[0,0,layer-1]));
  }

  double H = 0.;
  double zc = 0.;
  for (int l = 0; l < layer; l++) {
    H+=h[0,0,l];
  }
    ans=coeffFrotte(point,nl,s,h,layer,dg,rho,mu0,deltamu,I0)*pressionHydro(point,h,layer,rho)/shear;
  return ans;
}

double Ub(Point point,double HR, double I0, double mu0, double deltamu, double rho, double dg)
{
  double Ia = 0.; 
  double zc = 0.;
  double ans;
  for (int l = - point.l; l < nl - point.l; l++) {
    if (l < 0)
      zc += h[0,0,l];
  }
  zc += h[]/2.;

  Ia = (tan(slope)-mu0)*I0/(mu0+deltamu-tan(slope));
 
  ans = 2./3.*Ia/dg*sqrt(G*pow(HR,3)*cos(slope))*(1-pow((1-zc/HR),1.5));

  return ans;

}


double Ub2(double zc,double HR, double I0, double mu0, double deltamu, double rho, double dg)
{
  double Ia = 0.; 
  double ans;
  
  Ia = (tan(slope)-mu0)*I0/(mu0+deltamu-tan(slope));
 
  ans = 2./3.*Ia/dg*sqrt(G*pow(HR,3)*cos(slope))*(1-pow((1-zc/HR),1.5));

  return ans;


}

double Sb(double zc,double HR, double I0, double mu0, double deltamu, double rho, double dg)
{
  double Ia = 0.; 
  double ans;
  
  Ia = (tan(slope)-mu0)*I0/(mu0+deltamu-tan(slope));
 
  ans = Ia/dg*sqrt(G*pow(HR,1)*cos(slope))*pow((1-zc/HR),0.5);

  return ans;


}



/*- - - - - - - - - - - - - - - - - - - - - - - - - - - */
