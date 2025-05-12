/** # Linear forcing method

This routine aims at forcing turbulence in a domain 
based on a target $k$ and $\epsilon$.
The integral lenght scale $\ell$ is fixed by the forcing 
such that $\ell = 0.19\mathcal{L}$ with $\mathcal{L}$ the
domain lenght.
The implementation is based on this work (see [Rosales
& Meneveau, 2005](/src/references.bib#rosales2005)).

This implementation provides a better control of turbulent 
quantities if method 4,5,6 or 7 is used and it is also
adapted for two-phase flows with the computation of 
capillary forces.

*/

FILE * fStatF;

/**
First, we define macros in the case of single-phase flow where
rho and mu are not defined.
*/

#ifndef rho
# define rho(f) (rho1)
#endif
#ifndef mu
# define mu(f)  (mu1)
#endif

/**
Then, we define the turbulent quantities if they are not 
already defined elsewhere in the main code. Those are the
target quantities of the forcing
*/

#ifndef A0
#define A0 0.1 
#endif
#ifndef k0
#define k0 (27./2.*sq(0.38*pi*A0))
#endif
#ifndef eps0
#define eps0 (27.*sq(0.38*pi)*pow(A0,3))
#endif
#ifndef teddy
#define teddy (2./3.*k0/eps0)
#endif

#ifndef TFMETH
  #define TFMETH 1
#endif
#ifndef TFRHSNUM
  #define TFRHSNUM 2
#endif
#ifndef TFTAU
  #define TFTAU 67.
#endif
#ifndef TFC1
  #define TFC1 0.25
#endif
#ifndef TFC2
  #define TFC2 0.25
#endif

/**
If ABC forcing is on, override all other forcing parameters
and define $u_f$ and $\kappa_f$ the forcing amplitude and 
forcing wavenumber respectively.
*/
#ifndef TFABC
  #define TFABC 0
#endif
#ifndef uABC
  #define uABC sqrt(2/3*k0)
#endif
#ifndef LtABC
  #define LtABC (0.38*pi)
#endif

/**
For two-phase flows, TFALL is a parameter to force in both phase
or only in the carrier phase ($f=1$) and UMEAN
*/
#ifndef TFALL
  #define TFALL 1
#endif
#ifndef UMEAN
  #define UMEAN 2
#endif

/**
Some quantities need to be stored for all timestep.
*/
double ken = 0., An = 0., vdn = 0., wn = 0.;
double kenm1 = 0., Anm1 = 0., vdnm1 = 0., wnm1 = 0.;
double dkdt = 0., dvddt = 0.;
double dkdtr = 0., dvddtr = 0.;
double rhsk = 0., rhsvd = 0.;
double Pacc = 0., Phiacc = 0.;
double w1 = 0., w2 = 0;
coord ubar = {0.};
coord Facc = {0.};
coord AdimTmp = coord_null;
face vector Fturb[];

double c1 = TFC1/(TFC1 + TFC2), c2 = TFC2/(TFC1 + TFC2);

// Lundgren (2003)
#if TFMETH == 0 
coord compute_A(double ke, double vd, double w)
{
  double Alundgren = A0;
#if TFRHSNUM == 2
  double vdtot = eps0 - vd - dkdt + 2.*Anm1*kenm1;
  Alundgren = pow(vdtot/27/sq(Lt),1./3.);
#endif
  coord Adim;
  foreach_dimension() Adim.x = Alundgren;

  return Adim;
}
// Carroll (2013)
#elif TFMETH == 1 
coord compute_A(double ke, double vd, double w)
{
  double vdtot = eps0;
#if TFRHSNUM == 2
  vdtot = eps0 - vd - dkdt + 2.*Anm1*kenm1;
#endif
  coord Adim;
  foreach_dimension() Adim.x = 0.5*vdtot/ke;

  return Adim;
}
// Ketterl (2018)
#elif TFMETH == 2
coord compute_A(double ke, double vd, double w)
{
  coord Adim;
  foreach_dimension() Adim.x = max(0.,(k0 - ke)/(dt*k0));
  
  return Adim;
}
// Duret (2012)
#elif TFMETH == 3
coord compute_A(double ke, double vd, double w)
{
  double CkTmp = dkdt - 2.*Anm1*kenm1;
  double A = max(0.5/ke*((k0 - ke)/(3.*dt) - CkTmp),0.);
  coord Adim;
  foreach_dimension() Adim.x = A;
  // Store for next step
  keTmp = ke;
  ATmp = A;

  return Adim;
}
// Constant k - Bassenne (2016)
#elif TFMETH == 4
coord compute_A(double ke, double vd, double w)
{
  double rhsk = -vd + Pacc;
#if TFRHSNUM == 0
  rhsk = -vd;
#elif TFRHSNUM == 2
  rhsk = dkdt - 2.*Anm1*kenm1;
#endif
  double A = max(0.5/ke*(dkdtr - rhsk),0.);
  coord Adim;
  foreach_dimension() Adim.x = A;

  return Adim;
}
// Constant epsilon - Bassenne (2016)
#elif TFMETH == 5
coord compute_A(double ke, double vd, double w)
{
  double rhsvd = -w + Phiacc;
#if TFRHSNUM == 0
  rhsvd = -w;
#elif TFRHSNUM == 2
  rhsvd = dvddt - 2.*Anm1*vdnm1;
#endif
  double A = max(0.5/vd*(dvddtr - rhsvd),0.);
  coord Adim;
  foreach_dimension() Adim.x = A;

  return Adim;
}
// Hybrid method of Bassene (2016)
#elif TFMETH == 6
coord compute_A(double ke, double vd, double w)
{
  double rhsk = -vd + Pacc;
  double rhsvd = -w + Phiacc;
#if TFRHSNUM == 0
  rhsk = -vd;
  rhsvd = -w;
#elif TFRHSNUM == 2
  rhsk = dkdt - 2.*Anm1*kenm1;
  rhsvd = dvddt - 2.*Anm1*vdnm1;
#endif
  c1 = TFC1/(TFC1 + TFC2)*8*sq(ke)/(4*sq(ke) + 9*sq(vd*teddy));
  c2 = TFC1/(TFC1 + TFC2)*18*sq(vd*teddy)/(4*sq(ke) + 9*sq(vd*teddy));
  double A = max(0.5*c1/ke*(dkdtr - rhsk)
               + 0.5*c2/vd*(dvddtr - rhsvd),0.);
  coord Adim;
  foreach_dimension() Adim.x = A;

  return Adim;
}
// Constant general turbulence quantity (TFC1 and TFC2 to prescribe)
#elif TFMETH == 7
coord compute_A(double ke, double vd, double w)
{
  double rhsk = -vd + Pacc;
  double rhsvd = -w + Phiacc;
#if TFRHSNUM == 0
  rhsk = -vd;
  rhsvd = -w;
#elif TFRHSNUM == 2
  rhsk = dkdt - 2.*Anm1*kenm1;
  rhsvd = dvddt - 2.*Anm1*vdnm1;
#endif
  double A = max(0.5*c1/ke*(dkdtr - rhsk)
               + 0.5*c2/vd*(dvddtr - rhsvd),0.);
  coord Adim;
  foreach_dimension() Adim.x = A;

  return Adim;
}
#endif

/** Initialize all turbulent quantities and print header of the output file*/
event init (t=0)
{
  if (pid()==0) {
    fStatF = fopen("sf.csv","w");
#if dimension == 2
    fprintf (fStatF, "t,ken,vdn,w1,w2,wn,dkdt,dkdtr,dvddt,");
    fprintf (fStatF, "dvddtr,An,c1,c2,Pacc,Phiacc,Faccx,Faccy,ux,uy\n");
#elif dimension == 3
    fprintf (fStatF, "t,ken,vdn,w1,w2,wn,dkdt,dkdtr,dvddt,");
    fprintf (fStatF, "dvddtr,An,c1,c2,Pacc,Phiacc,Faccx,Faccy,Faccz,ux,uy,uz\n");
#endif
    fclose(fStatF);
  }

  // Compute volume fraction
  scalar ff[];
  scalar fff[];
#if TFALL
  foreach() ff[] = 0.;
  if (f.i)
    foreach() fff[] = f[];
  else {
    foreach() fff[] = 0.;
  }
#else // Apply only in carrier phase!
  foreach() fff[] = 1.;
  if (f.i)
    foreach() ff[] = f[];
  else {
    foreach() ff[] = 0.;
  }
#endif
  boundary ({ff});
  boundary ({fff});

  // Compute ubar
  ubar = coord_null;
  double vol = 0.;
  foreach(reduction(+:ubar) reduction(+:vol)) {
    vol += (1. - ff[])*dv();
    foreach_dimension() {
      // mean fluctuating kinetic energy
      ubar.x += (1. - ff[])*dv()*u.x[];
    }
  }
  ubar = div_coord(ubar,vol);

  scalar vareps[];
  vector Sx[], Sy[], Sz[];
  vector nuSx[], nuSy[], nuSz[];
  foreach(){
    mat3 GU, S;
    vec_grad(u,GU);
    double nu = mu(fff[])/rho(fff[]);
    vareps[] = 0.;
    einstein_sum(i,j){
      S.i.j = (GU.i.j  +  GU.j.i);
      vareps[] += 0.5*nu*(GU.i.j  +  GU.j.i)*(GU.i.j  +  GU.j.i);
    }
    Sx.x[] = S.x.x;
    Sx.y[] = S.x.y;
    Sx.z[] = S.x.z;
    Sy.x[] = S.y.x;
    Sy.y[] = S.y.y;
    Sy.z[] = S.y.z;
    Sz.x[] = S.z.x;
    Sz.y[] = S.z.y;
    Sz.z[] = S.z.z;
    nuSx.x[] = nu*S.x.x;
    nuSx.y[] = nu*S.x.y;
    nuSx.z[] = nu*S.x.z;
    nuSy.x[] = nu*S.y.x;
    nuSy.y[] = nu*S.y.y;
    nuSy.z[] = nu*S.y.z;
    nuSz.x[] = nu*S.z.x;
    nuSz.y[] = nu*S.z.y;
    nuSz.z[] = nu*S.z.z;
  }

  // Compute energetics
  ken = 0.;
  vdn = 0.;
  w1 = 0.;
  w2 = 0.;
  foreach(reduction(+:ken) reduction(+:vdn) 
          reduction(+:w1) reduction(+:w2)) {
    double dvc = (1. - ff[])*dv();
    double nu = mu(fff[])/rho(fff[]);
    mat3 GU, GSx, GSy, GSz, nuGSx, nuGSy, nuGSz;
    vec_grad(u,GU);
    vec_grad(Sx,GSx);
    vec_grad(Sy,GSy);
    vec_grad(Sz,GSz);
    vec_grad(nuSx,nuGSx);
    vec_grad(nuSy,nuGSy);
    vec_grad(nuSz,nuGSz);
    // Intermediate tensors because of overflows using triple-index repetition
    mat3 t1, t2;
    einstein_sum(i,j,k){
      ken += dvc*(u.i[] - ubar.i)*(u.i[] - ubar.i);
      t1.i.j = GU.i.k*GU.k.j;
      t2.i.j = (nuGSx.i.j*GSx.i.j + nuGSy.i.j*GSy.i.j + nuGSz.i.j*GSz.i.j);
    }
    einstein_sum(i,j){
      w1 += dvc*nu*(t1.i.j + t1.j.i)*(GU.i.j + GU.j.i);
      w2 += dvc*nu*t2.i.j;
    }
    vdn += dvc*vareps[];
  }
  ken = 0.5*ken/vol;
  vdn /= vol;
  w1 /= vol;
  w2 /= vol;
  wn = w1 + w2;
  double dtrel = teddy/TFTAU;
  kenm1 = ken;
  dkdt = 0;
  dkdtr = (k0 - ken)/dtrel;
  vdnm1 = vdn;
  dvddt = 0.;
  dvddtr = (eps0 - vdn)/dtrel;
  wnm1 = wn;
  An = sqrt(2./27.*kenm1)/Lt;
  Anm1 = An;
}

/** Compute the forcing term based on the energetics.*/
event compute_forcing(i++) {
/** First, store the energetics from the previous timestep.*/
  Anm1 = An;
  kenm1 = ken;
  vdnm1 = vdn;
  wnm1 = wn;

/** Compute volume fractions usefull for forcing in only one phase
ff and fff are masks which ensure that forcing is applied at the
right location (ff) with the right fluid properties (fff)
*/
  scalar ff[];
  scalar fff[];
#if TFALL
  foreach() ff[] = 0.;
  if (f.i)
    foreach() fff[] = f[];
  else {
    foreach() fff[] = 0.;
  }
#else // Apply only in carrier phase!
  foreach() fff[] = 1.;
  if (f.i)
    foreach() ff[] = f[];
  else {
    foreach() ff[] = 0.;
  }
#endif
  boundary ({ff});
  boundary ({fff});

/** Compute all energetics
*/
  ubar = coord_null;
  double vol = 0.;
  foreach(reduction(+:ubar) reduction(+:vol)) {
    vol += (1. - ff[])*dv();
    foreach_dimension() {
      // mean fluctuating kinetic energy
      ubar.x += (1. - ff[])*dv()*u.x[];
    }
  }
  ubar = div_coord(ubar,vol);

  scalar vareps[];
  vector Sx[], Sy[], Sz[];
  vector nuSx[], nuSy[], nuSz[];
  foreach(){
    mat3 GU, S;
    vec_grad(u,GU);
    double nu = mu(fff[])/rho(fff[]);
    vareps[] = 0.;
    einstein_sum(i,j){
      S.i.j = (GU.i.j  +  GU.j.i);
      vareps[] += 0.5*nu*(GU.i.j  +  GU.j.i)*(GU.i.j  +  GU.j.i);
    }
    Sx.x[] = S.x.x;
    Sx.y[] = S.x.y;
    Sx.z[] = S.x.z;
    Sy.x[] = S.y.x;
    Sy.y[] = S.y.y;
    Sy.z[] = S.y.z;
    Sz.x[] = S.z.x;
    Sz.y[] = S.z.y;
    Sz.z[] = S.z.z;
    nuSx.x[] = nu*S.x.x;
    nuSx.y[] = nu*S.x.y;
    nuSx.z[] = nu*S.x.z;
    nuSy.x[] = nu*S.y.x;
    nuSy.y[] = nu*S.y.y;
    nuSy.z[] = nu*S.y.z;
    nuSz.x[] = nu*S.z.x;
    nuSz.y[] = nu*S.z.y;
    nuSz.z[] = nu*S.z.z;
  }

  // Compute energetics
  ken = 0.;
  vdn = 0.;
  w1 = 0.;
  w2 = 0.;
  foreach(reduction(+:ken) reduction(+:vdn) 
          reduction(+:w1) reduction(+:w2)) {
    double dvc = (1. - ff[])*dv();
    double nu = mu(fff[])/rho(fff[]);
    mat3 GU, GSx, GSy, GSz, nuGSx, nuGSy, nuGSz;
    vec_grad(u,GU);
    vec_grad(Sx,GSx);
    vec_grad(Sy,GSy);
    vec_grad(Sz,GSz);
    vec_grad(nuSx,nuGSx);
    vec_grad(nuSy,nuGSy);
    vec_grad(nuSz,nuGSz);
    // Intermediate tensors because of overflows using triple-index repetition
    mat3 t1, t2;
    einstein_sum(i,j,k){
      ken += dvc*(u.i[] - ubar.i)*(u.i[] - ubar.i);
      t1.i.j = GU.i.k*GU.k.j;
      t2.i.j = (nuGSx.i.j*GSx.i.j + nuGSy.i.j*GSy.i.j + nuGSz.i.j*GSz.i.j);
    }
    einstein_sum(i,j){
      w1 += dvc*nu*(t1.i.j + t1.j.i)*(GU.i.j + GU.j.i);
      w2 += dvc*nu*t2.i.j;
    }
    vdn += dvc*vareps[];
  }
  ken = 0.5*ken/vol;
  vdn /= vol;
  w1 /= vol;
  w2 /= vol;
  wn = w1 + w2;
  dkdt = (ken - kenm1)/dt;
  dkdtr = TFTAU*(k0 - ken)/teddy;
  dvddt = (vdn - vdnm1)/dt;
  dvddtr = TFTAU*(eps0 - vdn)/teddy;

#if TFABC
  vector fABC[];
  double TABC = 1.5*LtABC/uABC;
  foreach() {
    fABC.x[] = uABC/TABC*(cos(y*2*pi/LtABC) + sin(z*2*pi/LtABC));
    fABC.y[] = uABC/TABC*(cos(z*2*pi/LtABC) + sin(x*2*pi/LtABC));
    fABC.z[] = uABC/TABC*(cos(x*2*pi/LtABC) + sin(y*2*pi/LtABC));
  }
  boundary({fABC});
  coord fmABC = {0.};
  foreach(reduction(+:fmABC)) {
    foreach_dimension() 
      fmABC.x += dv()*0.5*(ff[]*fABC.x[] + 0.5*(ff[-1]*fABC.x[-1] + ff[1]*fABC.x[1]));
  }
  fmABC = div_coord(fmABC,vol);
  foreach_face(){
    Fturb.x[] = (1. - 0.5*(ff[] + ff[-1]))
               *(0.5*(fABC.x[] + fABC.x[1]) - fmABC.x);
  }
#else
  // Compute forcing term
  AdimTmp = compute_A(ken,vdn,wn);
  An = AdimTmp.x;

#if UMEAN==0
  ubar = coord_null;
#endif

  // Update acceleration
  foreach_face(){
    Fturb.x[] = (1. - 0.5*(ff[] + ff[-1]))*AdimTmp.x
               *(0.5*(u.x[] + u.x[-1]) - ubar.x);
  }
#endif

}

event acceleration (i++) {
  // Compute mean acceleration due to capillary forces/gravity
  Pacc = 0.;
  Facc = coord_null;
  Phiacc = 0.;
  double vol = 0.;

  vector Fsigc[];
  foreach(){
    foreach_dimension(){
      Fsigc.x[] = 0.5*(a.x[] + a.x[1]);
    }
  }
  boundary ({Fsigc});
  
  foreach(reduction(+:Facc) reduction(+:Pacc) reduction(+:Phiacc) reduction(+:vol)) {
    vol += dv();
    mat3 GU, gradFsig;
    vec_grad(u,GU);
    vec_grad(Fsigc,gradFsig);
    foreach_dimension() {
      Facc.x += dv()*Fsigc.x[];
      Pacc += dv()*Fsigc.x[]*u.x[];
      // better approximation when grad is aligned with staggered acceleration
      gradFsig.x.x = (a.x[1] - a.x[])/Delta; 
    }
    einstein_sum(i,j){
      Phiacc += dv()*mu(f[])/rho(f[])*2*GU.i.j*gradFsig.i.j;
    }
  }
  Facc = div_coord(Facc,vol);
  Pacc /= vol;
  Phiacc /= vol;

  // Update acceleration
  foreach_face(){
    a.x[] += Fturb.x[];
#if UMEAN==2
    a.x[] -= Facc.x;
#endif
  }
}

event print_forcing (t=0; t += 0.1*teddy)
{
  if (pid()==0) {
    fStatF = fopen("sf.csv","a");
    fprintf(fStatF,"%g",t);
    fprintf(fStatF,",%g",ken);
    fprintf(fStatF,",%g",vdn);
    fprintf(fStatF,",%g",w1);
    fprintf(fStatF,",%g",w2);
    fprintf(fStatF,",%g",wn);
    fprintf(fStatF,",%g",dkdt);
    fprintf(fStatF,",%g",dkdtr);
    fprintf(fStatF,",%g",dvddt);
    fprintf(fStatF,",%g",dvddtr);
    fprintf(fStatF,",%g",An);
    fprintf(fStatF,",%g",c1);
    fprintf(fStatF,",%g",c2);
    fprintf(fStatF,",%g",Pacc);
    fprintf(fStatF,",%g",Phiacc);
    einstein_sum(i,j){
      fprintf(fStatF,",%g",Facc.i);
      fprintf(fStatF,",%g",ubar.i);
    }
    fprintf(fStatF,"\n");
    fclose(fStatF);
  }
}

/**
## References

~~~bib
@article{lundgren2003linearly,
  title={Linearly forced isotropic turbulence},
  author={Lundgren, Thomas S},
  journal={Center for Turbulence Research Annual Research Briefs 2003},
  year={2003}
}
@article{carroll2013proposed,
  title={A proposed modification to Lundgren's physical space velocity forcing method for isotropic turbulence},
  author={Carroll, Phares L and Blanquart, Guillaume},
  journal={Physics of Fluids},
  volume={25},
  number={10},
  year={2013},
  publisher={AIP Publishing}
}
@article{bassenne2016constant,
  title={Constant-energetics physical-space forcing methods for improved convergence to homogeneous-isotropic turbulence with application to particle-laden flows},
  author={Bassenne, Maxime and Urzay, Javier and Park, George I and Moin, Parviz},
  journal={Physics of Fluids},
  volume={28},
  number={3},
  year={2016},
  publisher={AIP Publishing}
}
@article{naso2010interaction,
  title={The interaction between a solid particle and a turbulent flow},
  author={Naso, Aurore and Prosperetti, Andrea},
  journal={New Journal of Physics},
  volume={12},
  number={3},
  pages={033040},
  year={2010},
  publisher={IOP Publishing}
}
@article{yao2021deagglomeration,
  title={Deagglomeration of cohesive particles by turbulence},
  author={Yao, Yuan and Capecelatro, Jesse},
  journal={Journal of Fluid Mechanics},
  volume={911},
  pages={A10},
  year={2021},
  publisher={Cambridge University Press}
}
~~~
*/