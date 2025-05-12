#include "grid/multigrid1D.h"
#include "spherisym.h"
#include "compressible/thermal.h"
#include "tension.h"
#include "compressible/tension.h"
#include "compressible/NASG.h"

#define LEVEL 18

double rhoL = 1., rhoR = 1.787929636728e-3;
double Rbub = 4.5e-4;
double tend = 1.5;
double tr;

double pdim = 101325.;
double rhodim = 996.56;
double Tdim = 300.;
double Rdim = 0.01;

double FREQ = 26.326163109183227;
double OMEGA = 165.4121612420333;
double CLL = 155.3326;
double KWL = 1.064890185589073;

uf.n[right] = neumann(0.);
p[right]    = dirichlet(1. + 0.98969110902632*sin(OMEGA*t));
q.n[right]  = neumann(0.);

uf.n[left] = 0.;

int main() {
  L0 = 1.;
  CFLac = 0.5;
  DT = HUGE [0];
  
  double Reynolds = 100487.035;
  mu1 = 1./Reynolds;
  mu2 = 0.01*mu1;

  double Weber = 14072.91667;
  f.sigma = 1./Weber;
  f.gradient = zero;

  gamma1 = 1.187;
  PI1 = 7028e5/pdim;
  b1 = 6.61e-4*rhodim;
  q1 = -1177788*rhodim/pdim;

  gamma2 = 5./3.;

  cv1 = 3610*rhodim*Tdim/pdim; cv2 = 313.17*rhodim*Tdim/pdim;
  cp1 = 4285*rhodim*Tdim/pdim; cp2 = 523.*rhodim*Tdim/pdim;

  kappa1 = 0.61032227156167/(Rdim/Tdim*sqrt(cube(pdim)/rhodim));
  kappa2 = 20.25e-3/(Rdim/Tdim*sqrt(cube(pdim)/rhodim));

  tr = 1./FREQ;
  
  tend *= tr;

  init_grid(1 << LEVEL);

  TOLERANCE = 1e-6;
  
  run();
}

event init (t = 0) {
  if (!restore (file = "restart")) {
  
    foreach() {
      f[] = (Rbub - x) > Delta/2. ? 0. :
	fabs(Rbub - x) < Delta/2. ? 1. - 0.5 - (Rbub - x)/Delta :
	(Rbub - x) == Delta/2. ? 0. :
	(Rbub - x) == -Delta/2. ? 1. :
	1.;
      
      frho1[]  = f[]*rhoL;
      frho2[]  = (1. - f[])*rhoR;

      double pL = 1.;
      p[] = f[]*pL + (1. - f[])*(pL + 2.*f.sigma/Rbub);

      double fc = clamp (f[],0.,1.);
      double rhocpmcvavg = (cp1 - cv1)*frho1[] + (cp2 - cv2)*frho2[];
      double const1 = (fc - frho1[]*b1) + (1. - fc - frho2[]*b2);
      double const2 = (fc - frho1[]*b1)*PI1 + (1. - fc - frho2[]*b2)*PI2;
      T[] = (const1*p[] + const2)/rhocpmcvavg;

      q.x[] = frho1[]*0.98969110902632*(cos(KWL*x) - sin(KWL*x)/KWL/x)/CLL/sin(KWL)/x;
    
      fE1[]   = f[] ? sq(q.x[])/frho1[]/2. + (pL + gamma1*PI1)/(gamma1 - 1.)*(f[] - frho1[]*b1) + frho1[]*q1 : 0.;
      fE2[]   = (1. - f[])*(pL/(gamma2 - 1.));
    }
  }
  fprintf(stderr,"OMEGA = %g\n",OMEGA);
}

event centroid (i++) {
  scalar ff[];
  double radius = 0.;
  foreach(reduction(+:radius)) {
    ff[] = 1. - f[];
    radius += ff[]*Delta;
  }

  if (pid() == 0) {
    char name[80];
    sprintf(name,"radius.txt");
    FILE * fp = fopen(name,"a");
    char str[80];
    sprintf(str,"%g %g\n",t/tr,radius/Rbub);
    fputs(str,fp);
    fclose(fp);
  }
}

event temperature (i++) {
  double temperature = 0.;
  for (scalar s in {T})
    temperature = interpolate (s, 0.);

  if (pid() == 0) {
    FILE * fp = fopen("signal.txt","a");
    char str[80];
    sprintf(str,"%g %g\n",t/tr,temperature);
    fputs(str,fp);
    fclose(fp);
  }
}

event logfile (i++) {
  stats sp = statsf (p);
  stats su = statsf (q.x);
  stats sT = statsf (T);
  fprintf (stderr,"t = %g, i = %d, dt = %g, min(p) = %g, max(p) = %g, min(T) = %g, max(T) = %g, min(u) = %g, max(u) = %g\n", t/tr, i, dt/tr, sp.min, sp.max, sT.min, sT.max, su.min, su.max);
}

event output (t = 0.5/FREQ; t += 0.001/FREQ) {
  scalar pdump[];
  
  foreach() {
    double Ek = 0.;
    foreach_dimension()
      Ek += sq(q.x[]);

    double fc = clamp (f[],0.,1.);
    double invgammaavg = (fc - frho1[]*b1)/(gamma1 - 1.) +
      (1. - fc - frho2[]*b2)/(gamma2 - 1.);
    double PIGAMMAavg = PI1*gamma1*(fc - frho1[]*b1)/(gamma1 - 1.) + frho1[]*q1 +
      PI2*gamma2*(1. - fc - frho2[]*b2)/(gamma2 - 1.) + frho2[]*q2;

    pdump[] = (fE1[] + fE2[] - Ek/(frho1[] + frho2[])/2. - PIGAMMAavg)/invgammaavg;
  }
  
  double data[2*N];
  for (int i = 0; i < N; i++) {
    double xx = L0*(1./2. + i)/N;
    Point point = locate(xx);
    for (scalar s in {pdump})
      data[i] = point.level >= 0 ? s[] : nodata;
    for (scalar ss in {T})
      data[i + N] = point.level >= 0 ? ss[] : nodata;
  }
  
  if (pid() == 0) {
    @if _MPI
      MPI_Reduce (MPI_IN_PLACE, &data, 2*N, MPI_DOUBLE, MPI_MIN, 0,
  		  MPI_COMM_WORLD);
    @endif
    char name[80];
    sprintf(name,"data-%g",t*FREQ);
    FILE * fpp = fopen (name,"w");
    for (int i = 0; i < N; i++) {
      double xx = L0*(1./2. + i)/N;
      for (scalar s in {pdump})
  	fprintf(fpp,"%g %g ",xx/Rbub,data[i]);
      for (scalar ss in {T})
	fprintf(fpp,"%g",data[i + N]);
      fputc ('\n',fpp);
    }
    fflush(fpp);
  }
  @if _MPI
  else
    MPI_Reduce (&data, NULL, 2*N, MPI_DOUBLE, MPI_MIN, 0,
  		  MPI_COMM_WORLD);
  @endif
}

event ending (t = tend) {
  return 1.;
}
