/**
#Subroutines to obtain the theoretical solution of the Sod's tube problem
*/

extern double rhoL, rhoR, pL, pR, gammaL, gammaR;

double get_u4(double p1, double p2, double g, double r) {

  double m = (g - 1)/(g + 1);
  return (p1 - p2)*sqrt((1.-m)/(r*(p1 + m*p2))); 
}

double get_du4dp(double p1, double p2, double g, double r) {

  double m = (g - 1)/(g + 1);
  return sqrt((1.-m)/(r*(p1 + m*p2)))*(1. - 1./2.*(p1-p2)/(p1 + m*p2)); 
}


double get_u2(double p1, double p2, double g, double r) {

  double m = (g - 1)/(g + 1);
  return (pow(p1,(g-1)/(2*g)) - pow(p2,(g-1)/(2*g)))*sqrt((1.-m*m)*pow(p1,1/g)/(r*m*m)); 
}

double get_du2dp(double p1, double p2, double g, double r) {

  double m = (g - 1)/(g + 1);
  return - (g-1)/(2*g)*pow(p2,-(g+1)/(2*g))*sqrt((1.-m*m)*pow(p1,1/g)/(r*m*m)); 
}


double get_p3 ()
{

  double a,b;
  a = pR;
  b = (pR+pL)/2;
  double p3 = (a + b)/2;

  int Imax = 1000, i = 1;
  double error = 1.e10, tol = 1.e-3;

  //Newton-Raphson algorithm
  while (error > tol && i < Imax) {
    double u4 = get_u4(p3,pR,gammaR,rhoR);
    double u2 = get_u2(pL,p3,gammaL, rhoL);
    error = fabs(u4 - u2);
    p3 -= (u4 - u2)/(get_du4dp(p3,pR,gammaR,rhoR)- get_du2dp(pL,p3,gammaL, rhoL));
    if (p3 < a)
      p3 = a;
    else if (p3 > b)
      p3 = b;

    i++;
  }

  return p3;
}

struct SodSol{
    double pe, rhoe, ue; 
};

struct SodSol Sod_theory ( double chi ) 
{
    struct SodSol ss; 
    double pe, rhoe, ue; 
  
    double p3 = get_p3();
    double p4 = p3;
    double u4 = get_u4(p4,pR,gammaR,rhoR);
    double u3 = u4;
    double mR = (gammaR - 1)/(gammaR + 1);
    double mL = (gammaL - 1)/(gammaL + 1);
    double rho4 = rhoR*(p4 + mR*pR)/(pR + mR*p4);
//  double rho3 = rhoL*pow(p3/pL,1./gammaL);
    double ushock =  u4*rho4/(rho4-rhoR);

    double csonL = sqrt(gammaL*pL/rhoL);
    double chi1 = -csonL;
    double chi2 = (u3/(1. - mL) - csonL);


    if (chi < chi1) {

      rhoe = rhoL;
      pe = pL;
      ue = 0.;

    }
    else if (chi < chi2) {

      ue = (1. - mL)*(chi + csonL);
      rhoe = pow(pow(rhoL,gammaL)/(gammaL*pL)*pow(ue - chi,2),1./(gammaL-1.));
      pe = pL*pow(rhoe/rhoL, gammaL);

    }
    else if (chi < u4) {

      ue = (1. - mL)*(chi2 + csonL);
      rhoe = pow(pow(rhoL,gammaL)/(gammaL*pL)*pow(ue - chi2,2),1./(gammaL-1.));
      pe = pL*pow(rhoe/rhoL, gammaL);

    }
    else if (chi < ushock) {

      ue = u4;
      pe = p4;
      rhoe = rho4;

    } else {

      ue = 0.;
      pe = pR;
      rhoe = rhoR;

    }

    ss.pe = pe;
    ss.rhoe = rhoe;
    ss.ue = ue;
    return ss;

}