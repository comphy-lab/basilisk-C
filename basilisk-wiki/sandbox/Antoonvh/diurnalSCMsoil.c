/**
# An atmospheric single-column model over a condictive soil 

This page is derived from [the scm with a lumped parameter description](lumpedscm.c), yet here the term $G$ is based on a resolved soil profile. 

More comments will follow...
*/

#include "grid/bitree.h"
#include "run.h"
#include "diffusion.h"

//#define fris(Ri) (sq((1 - (Ri/0.20)))*(Ri < 0.20)) // Critical Ri
#define fris(x) (exp(-10.*x))                        // Exponential
//#define fris(x) (1./(1. + (10*x*(1.+8.*x))))       // Long Tail
#define friu(Ri) (sqrt(1. - (18.*Ri)))               // Holtslag en Boville 1992
                                                     // Note: no surface f(Ri)
#define BSURF ((b[1] - b[]*c[level])/(1. - c[level]))
#define GFLX (-Lambda*(BSURF - (1.5*bsoil[0] - 0.5*bsoil[1])))
#define Qn (max(B0*sin(2.*M_PI*t/T), B1))

scalar u[], v[], b[];
double Lc = 1030.;
double T = 24.*3600.;
double Lambda, B0, B1, f, Nv, Ugeo, zom, zoh;
double Pi1 = -6;
double Pi2 = 2160.;
double Pi3 = 10.;
double Pi4 = 5366./3.;
double Pi5 = 4;
double Pi6 = 5150.;
double Pi7 = 1.;

double bo = 0., k = 0.4;
double c[20], lut[20], Lmix;

int ni, maxlevel = 9;
mgstats mgb;
double Umax, inv, u40;
int ns = 40;
double dzs = .05; // 5cm soil layers
double kappas = 0.5E-6; // Soil conductivity
double rat = 1000;
double *bsoil; 

b[left] = dirichlet(BSURF);
b[right] = dirichlet(bo + sq(Nv)*x);
u[left] = dirichlet(0.);
v[left] = dirichlet(0.);

int main(){
  bsoil = malloc(ns*sizeof(double));
  memset(bsoil, 0, ns*sizeof(double));
  Nv = Pi2/T;
  f = Pi3/T;
  B0 = sq(Lc*Nv)*M_PI/(2.*T);
  B1 = B0/Pi1;
  Lambda = sqrt(B0*T)/Pi4;
  //Lambda = 1.;
  zom = Lc / Pi6;
  zoh = zom * Pi7;
  L0 = 3.*Lc;
  Ugeo = Pi5*pow(B0*Lc, 1./3.);
  init_grid(128);
  run();
  
}

event init(t = 0){
  ni = 0;
  inv = Umax = u40 = 0.;
  DT = T/(3600.*24.);
  TOLERANCE = 10E-4;
  refine(x < (Lc/10.) && level < maxlevel);
  foreach(){
    u[] = Ugeo;
    v[] = 0.;
    b[] = bo + sq(Nv)*x;
  }
  for (int j = 1; j <= maxlevel; j++){
    double d = (L0/((double)(1 << j)))/zoh;
    double d2 = (L0/((double)(1 << j)))/zom;
    c[j] = (log(4.*d) - 1.)/(log(d) - 1.);
    lut[j] = sq(k/(log(d2) - 1.));
  }
}
/**
This event time integrated the diffusion equation for the soil heat. 
*/
event soil(i++){
  double G = 0;
  foreach_boundary(left)
    G = GFLX;
  for (int j = 0; j < ns; j++){
    if (j == 0)
      bsoil[j] += dt/dzs*(-G/rat - kappas*(bsoil[j]  - bsoil[j + 1])/dzs);
    else if (j < ns)
      bsoil[j] += kappas*(bsoil[j-1] - 2*bsoil[j] + bsoil[j+1])*dt/sq(dzs);
    else // (j == ns)
      bsoil[j] += kappas*(bsoil[j-1] - 2*bsoil[j] + 0)*dt/(sq(dzs));
  }
}

event diff(i++){
  scalar rx[],ry[],rb[];
  face vector kh[];
  double B = 0;
  double ws = 0;
  Lmix = 0.;
  foreach()
    Lmix += (b[] - x*sq(Nv)) * Delta;
  if (Lmix > 0.)
    Lmix = sqrt (Lmix*2./sq(Nv));
  printf("%g Lmix = %g\n",t,Lmix);
  foreach(){
    rx[] = f*v[];
    ry[] = f*(Ugeo - u[]);
    rb[] = 0.;
    if (x < Delta){
      B = (Qn + GFLX);
      rb[] += B / Delta;
      rx[] += -sign(u[])*lut[level]*sq(u[])/Delta;
      ry[] += -sign(v[])*lut[level]*sq(v[])/Delta;
    }
  }
  printf("%g B=%g\n",t,B);
  if (B > 0 && Lmix > 0)
    ws = pow(B*Lmix, 1./3.);
  printf("%g ws=%g\n",t,ws);
  foreach_face(){
    double sqd = (sq((u[] - u[-1])/(Delta)) + sq((v[] - v[-1])/(Delta)));
    double Ri = ((b[] - b[-1])/(Delta))/(sqd + 0.00001);
    double fRi;
    if (Ri < 0)
      fRi = friu(Ri);
    else
      fRi = fris(Ri);
    double l = min(k*x, (70./1030.)*Lc);
    double fraction = 0;
    if (Lmix > 0)
      fraction = 3.*(x/Lmix * sq(1. - x/Lmix))*(x < Lmix);
    double Vs = sqrt(fraction*sq(ws) + sq(l)*sqd*sq(fRi));
    kh.x[] = l*Vs;
  }
  boundary(all);
  dt = dtnext(DT);
  int n = 0;
  mgb = diffusion(u, dt, kh, r = rx);
  n += mgb.i;
  mgb = diffusion(v, dt, kh, r = ry);
  n += mgb.i;
  mgb = diffusion(b, dt, kh, r = rb);
  n += mgb.i;
  if (n>10)  //Quickly reduce the timestep if things get rough
    DT=max(DT/(1+((double)n/10.)), T/(24.*3600.));
  if (n<5)   //Slowly increase the timestep when time integration is easy.
    DT=min(DT*(1+((double)n/100.)), T/(24.*3600.));
}

event adapt(i++){
  double ue = Ugeo/20.;
  double be =  sq(Nv)*Lc/50.;
  adapt_wavelet({u,v,b}, (double[]){ue, ue, be}, maxlevel);
}

event profile(t += T/(24*15)){
  char fname[99];
  sprintf(fname, "ProfsPi5%gatmos",Pi5);
  static FILE * fp = fopen(fname, "w");
  sprintf(fname, "ProfsPi5%gsoil",Pi5);
  static FILE * fp2 = fopen(fname, "w");
  
  foreach(){
    if (x < Delta){
      fprintf(fp, "0 0 0 %g \n", BSURF);
      fprintf(fp2, "%g\t", BSURF);
    }
    fprintf(fp, "%g %g %g %g \n", x, u[], v[], b[]);
  }
  
  for (int j = 0; j < ns; j++)
    fprintf(fp2, "%g\t", bsoil[j]);
  fprintf(fp2, "\n");
  fflush(fp2);
  
}

double xuvb[4][1000] = {0};
double trec = 0;
event llj (t = T/2; i += 5){
  bool new_record = false;
  foreach()
    if (sq(u[]) + sq(v[]) > sq(Umax))
      new_record = true;
  if (new_record){
    trec = t;
    int j = 0;
    foreach(){
      xuvb[0][j] = x;
      xuvb[1][j] = u[];
      xuvb[2][j] = v[];
      xuvb[3][j++] = b[];
    }
    xuvb[0][j] = -1.;
  }
}   


event timeseries(i += 10){
  char fname[99];
  sprintf(fname, "timeseriesPi5%gsoil",Pi5);
  static FILE * fp = fopen(fname, "w");
  double B, G, bs, bmix;
  Lmix = 0;
  foreach()
    Lmix += (b[] - x*sq(Nv)) * Delta;
  if (Lmix > 0.)
    Lmix = sqrt (Lmix*2./sq(Nv));
  bmix = Lmix * sq(Nv);
  double Lambdae;
  foreach_boundary(left){
    G = GFLX;
    Lambdae = G/(BSURF);
    B = Qn + G;
    bs = BSURF; 
  }
  if (i == 0)
    fprintf(fp, "t\tQn\tG\tB\tbs\tbmix\tLmix\n");
  fprintf(fp, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", t, Qn, G, B, bs, bmix, Lmix, Lambdae);
  
}

event stop (t = T){
  return 1;
}


