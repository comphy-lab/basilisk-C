/** ## Nanofluid convection "two-phase"
*/
/** ## Model equations
*/
/**
Continuity equation for nanofluid :
$$
        \widetilde{\nabla} . \widetilde{u} = 0
$$
The momentum equation for nanofluid : 
$$
        \partial_{\tau}\widetilde{u} + \widetilde{u}.\widetilde{\nabla}\widetilde{u} = - \widetilde{\nabla} \widetilde{p} + \frac{Pr_{nf}}{\sqrt{Ra_{nf}}}\widetilde{\nabla}^{2}\widetilde{u} + Pr_{nf} \widetilde{\theta}
$$
The energy equation for the nanofluid :
$$
        \partial_{t} \widetilde{\theta} + \widetilde{u}.\widetilde{\nabla} \widetilde{\theta} = \frac{\widetilde{\nabla}^{2}\widetilde{\theta}}{\sqrt{Ra_{nf}}} + \frac{\tau \widetilde{\nabla} \Phi . \widetilde{\nabla} \widetilde{\theta}}{Le_{nf}\sqrt{Ra_{nf}}} + \frac{\tau S_{t}\widetilde{\nabla} \widetilde{\theta} . \widetilde{ \nabla} \widetilde{\theta}}{\sqrt{Ra_{nf}}}
$$
Continuity equation for nanoparticles :
$$
        \partial_{t}\Phi + \widetilde{u}. \widetilde{\nabla} \Phi = \frac{1}{\sqrt{Ra_{nf}}}[\frac{\widetilde{\nabla}^{2}\Phi}{Le_{nf}} + S_{t} \widetilde{\nabla}^{2} \widetilde{\theta}]
$$
*/

/** ## Dimensionless Parameters
*/
/**
$$
        Ra_{nf} = \frac{g\beta_{nf} \Delta T H^{3}}{\nu_{nf}\kappa_{nf}}, \quad Pr_{nf} = \frac{\nu_{nf}}{\kappa_{nf}}, \quad Le_{nf} = \frac{k_{nf}}{\rho_{nf}Cp_{nf}D_{B}}, \quad S_{T} = \frac{D_{t} \Delta T}{\kappa_{nf} T_{0} \phi}, \quad \tau = \frac{\rho_{p} Cp_{p} \phi}{\rho_{nf} Cp_{nf}}
$$

*/

#define T0 (45.+273.15)
/** Characteristic scales (water + aluminia particles) */
#define Rhof (2446. -20.674*(T0) + 0.11576*sq(T0) -3.12895e-4*cube(T0) + 4.0505e-7*pow(T0,4) -2.0546e-10*pow(T0,5))
#define Cpf (exp((8.29041-0.012557*(T0))/(1.-1.52373e-3*(T0))))
#define kf (-0.76761 + 7.535211e-3*(T0) -0.98249e-5*sq(T0))
#define Betaf (21.e-5)
#define muf (2.414e-5*pow(10,247.8/((T0)-140.)))
#define Nuf (muf/Rhof)
#define Kappaf (kf/(Rhof*Cpf))


#if 0
/** Characteristic scales (temperature depend) */
#define Rhof (983.37)//998.30
#define Cpf (4193.92)//4183.
#define kf (0.652)//0.604
#define Betaf (21.e-5)
#define muf(T) (2.414e-5*pow(10,247.8/((T)-140.)))
#define Nuf (muf/Rhof)
#define Kappaf (kf/(Rhof*Cpf))
#endif

#define Rhop (3880.)//3600
#define Cpp (765.)
#define kp (36.)//(36.)
#define Betap (8.1e-6)
//#define dp (25.e-9) 
#define kb (1.38066e-23)
#define Kappap (kp/(Rhop*Cpp))

/** nanofluid mass density */
#define RHOnf ((1-Phi)*Rhof + Phi*Rhop)

/** nanofluid Heat capacity */
#define CPnf (((1-Phi)*Rhof*Cpf + Phi*Rhop*Cpp)/RHOnf)

/** nanofluid thermal expansion */
#define BETAnf (((1-Phi)*Betaf*Rhof + Phi*Betap*Rhop)/RHOnf)

/** nanofluid effective thermal conductivity */
#define knf ((1+4.4*pow(Rep,0.4)*pow(Prf,0.66)*pow(T0/Tfr,10)*pow(kp/kf,0.03)*pow(Phi,0.66))*kf)
#define Rep ((2*Rhof*kb*T0)/(M_PI*sq(muf)*dp))
#define Prf (Nuf/Kappaf)
#define Tfr (273.15)

/** nanofluid diffusivity */
#define KAPPAnf (knf/(CPnf*RHOnf))

/** nanofluid effective dynamic viscosity */
#define MUnf (muf + (Rhop*VB*sq(dp))/(72*C*deltaa))//attention dp en nm pour C
#define VB (sqrt((18.*kb*T0)/(M_PI*Rhop*dp))*1./dp)
#define deltaa (pow(M_PI/(6*Phi),1./3)*dp)
#define C (((-0.000001133*dp*1e9-0.000002771)*Phi+(0.00000009*dp*1e9-0.000000393))/muf)

/** nanofluid cinematic viscosity */
#define NUnf (MUnf/RHOnf)

/** Einstein + Thermophoretic Diffusion coefficient*/
#define DB (((kb*T0)/(3*M_PI*muf*dp)))
#define Dt ((0.26*kf)/(2*kf+kp)*(muf/Rhof)*Phi)

/** Temperature boundary + Temperature average*/
#define Thot (T0+5.)
#define Tcold (T0-5.)
#define deltaTemp (Thot - Tcold)
#define Tbot (((Thot)-T0)/deltaTemp)
#define Ttop (((Tcold)-T0)/deltaTemp)
#define Tav ((Tbot+Ttop)/2.)


/** This is the maximum and minimum level of refinement*/
#define MINLEVEL 5
#define err 1e-3

/** End time simulation */
#define EndTime 600.
#define nanofluid 1

//#include "display.h"
#include "convection_nanofluid.h"
#include "navier-stokes/perfs.h"
#include "view.h"
#include "output_vtu_foreach.h"
#include "parameters.h"
#include "output_gnuplot.h"
#include "profile6.h"

void streamfunction (vector u, scalar psi){
  scalar omega[];
  vorticity (u, omega);
  boundary ({psi, omega});
  poisson (psi, omega);
}

int main(int argc, char *argv[]) {
  Ranf = 1e4;
  MAXLEVEL = 8;
  if(argc > 1)
    Ranf = atof(argv[1]);
  if(argc > 2)
    MAXLEVEL = atoi(argv[2]);
#if nanofluid
  Phi = 0.01;
  if(argc > 3)
    Phi = atof(argv[3]);
  dp = 13.e-9;
  if(argc > 4)
    dp = atof(argv[4]);
  Prnf = NUnf/KAPPAnf;
  Lenf = KAPPAnf/DB;
  St = (Dt*deltaTemp)/(KAPPAnf*T0*Phi);
  tau = (Rhop*Cpp*Phi)/(RHOnf*CPnf);
  Np.gradient = minmod;
#endif
  L0 = 1.;
  X0 = Y0 = -0.5;
  TOLERANCE = 1e-5;
  init_grid(1<<(MAXLEVEL));
  mu = muc; //viscosity
  a = av; // acceleration
  run();
}

/** 
    ## Boundaries Conditions 
*/

T[top] = dirichlet(Ttop);
T[bottom] = dirichlet(Tbot);
T[left] = neumann(0.);
T[right] = neumann(0.);

u.n[top] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.n[right] = dirichlet(0.);
u.n[left] = dirichlet(0.);

u.t[top] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.t[left] = dirichlet(0.);

#if nanofluid
Np[left] = neumann(0.);
Np[right] = neumann(0.);
Np[top] = neumann(-D4ORDTy*St*Lenf);
Np[bottom] = neumann(-D4ORDTy*St*Lenf);
#endif

event restoring(i=0){
  if (!restore ("dump")){
    restore("dump");
  }
}

/**
   ## Initial conditions
   Initial conditions correspond to a linear temperature and no motion
*/

event init (t=0) {
#if nanofluid
  refine (y < -(L0/2.-2*Delta) && y>(L0/2.-2*(Delta)) && level <= MAXLEVEL);
  foreach(){
    foreach_dimension()
      u.x[] = noise()*0.001;
    T[] = Tbot -(y+0.5)*(Tbot-Ttop);
    Np[] = 1.;
  }
  boundary ({T,u,Np});
#else
  foreach(){
    T[] = (Tbot -(y+0.5)*(Tbot-Ttop));
    foreach_dimension()
      u.x[] = noise()*0.01;
  }
  boundary ({T,u});
#endif
  system("mkdir -p paraview");
  //DT=0.1;
}

/** 
    ## Adaptative mesh
*/
#define dA Delta
event couche(i++){
  foreach()
    CLT[] = 0.;
  foreach_boundary(bottom)
    CLT[] = (L0*L0*(Tbot-Ttop)*Delta)/(dA*(T[0,-1] - T[]));
  foreach_boundary(top)
    CLT[] = (L0*L0*(Tbot-Ttop)*Delta)/(dA*(T[] - T[0,1]));
}
event adapt (i++){
#if nanofluid
  adapt_wavelet ((scalar*){Np, T, CLT, u},
		 (double[]){err*0.1, err, err, err, err}, MAXLEVEL, MINLEVEL);
#else
  adapt_wavelet ((scalar*){T, CLT, u},
		 (double[]){err, err, err, err}, MAXLEVEL, MINLEVEL);
#endif
}

event statistics(t+=1.; t <= EndTime){
#if nanofluid
  stats s = statsf(Np);
#endif
  stats d = statsf(u.x);
  stats f = statsf(u.y);
  char name[80];
  sprintf(name, "values%2.2e.dat", Ranf);
  static FILE * fc = fopen (name, "w");
  if (i == 0)
#if nanofluid
    fprintf(fc, "#t Np.min Np.max Np.sum ux.max uy.max\n");
  fprintf (fc, "%f %g %g %g %g %g\n", t, s.min, s.max, s.sum, d.max, f.max);
#else
    fprintf(fc, "#t ux.max uy.max\n");
  fprintf (fc, "%f %g %g\n", t, d.max, f.max);
#endif
  fflush(fc);
}

#if 1
scalar yref[];
yref[top] = dirichlet(-0.5);
yref[bottom] = dirichlet(0.5);
#include "available_potential.h"
event reference (i++,last){
	reference_height(yref,T,-0.5,0.5,N);
	boundary ({yref});
}

event time_series (t += 0.05; t <= EndTime){

	vector dT[], dyref[];
	gradients ({T}, {dT});
	gradients ({yref}, {dyref});

	tensor du[];
	foreach(){
		du.x.y[] = (u.x[0,1] - u.x[0,-1])/2./Delta;
		du.y.x[] = (u.y[1,0] - u.y[-1,0])/2./Delta;
		du.x.x[] = (u.x[1,0] - u.x[-1,0])/2./Delta;
		du.y.y[] = (u.y[0,1] - u.y[0,-1])/2./Delta;
	}

	double nu_top = 0., nu_bot = 0.;
	foreach_boundary (top,reduction(+:nu_top))
		nu_top += dA*(T[] - T[0,1])/Delta;

	foreach_boundary (bottom,reduction(+:nu_bot))
		nu_bot += dA*(T[0,-1] - T[])/Delta;

	double ekin = 0., epot = 0., bpot = 0.;
	foreach(reduction(+:ekin), reduction(+:epot), reduction(+:bpot)){
		epot += -dv() * Prnf * T[] * y;
		bpot += -dv() * Prnf * T[] * yref[];
		foreach_dimension()
			ekin += dv() * 0.5 * sq(u.x[]);
	}

	double nu_vol=0., nu_eps=0., nu_tmp=0., nu_mix=0.;
	foreach(reduction(+:nu_vol), reduction(+:nu_eps), reduction(+:nu_tmp), reduction(+:nu_mix)){
		nu_vol += dv() * (sqrt(Ranf)*u.y[]*T[] - dT.y[]);
		nu_eps += dv() * (sq(du.x.x[]) + sq(du.x.y[]) +  sq(du.y.x[]) + sq(du.y.y[]));
		foreach_dimension(){
			nu_tmp += dv() * (dT.x[]*dT.x[]);
			nu_mix += dv() * (dT.x[]*dyref.x[]);
		}
	}

	FILE * fp = fopen("rb2d.dat", "a");
	fprintf (fp, "%f %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g \n", t, ekin, epot, bpot, nu_top, nu_bot, nu_vol, nu_eps, nu_tmp, nu_mix);
	fclose (fp);
}
#endif

/** 
    ## Nusselts
*/

#include "nusselts.h"
event Nusselts_time_series (t += 1.; t <= EndTime){

  double nu_top = nusselt_top(T, Ranf);
  double nu_bot = nusselt_bot(T, Ranf);
  double nu_right = nusselt_right(T);
  double nu_left = nusselt_left(T);
  double nu_pheno_top = nusselt_pheno_top(T);
  double nu_pheno_bot = nusselt_pheno_bot(T);
  double nu_buon_top = nusselts_buon_top(T, Ranf, Np);
  double nu_buon_bot = nusselts_buon_bot(T, Ranf, Np);

/**
  nu_bot=nu_bot/(L0*(Tbot-Ttop));
  nu_top=nu_top/(L0*(Tbot-Ttop));
  nu_right=nu_right/(L0*(Tbot-Ttop));
  nu_left=nu_left/(L0*(Tbot-Ttop));
*/ 
  
  char name[80];
  sprintf(name, "Nusselts_%2.2e.dat", Ranf);
  static FILE * fp = fopen(name, "w");
  if (i == 0)
    fprintf(fp, "#t nu_{top} nu_{bot} nu_{left} nu_{right} nu_phe_{top} nu_phe_{bot} nu_buon_{top} nu_buon_{bot}\n");
  fprintf (fp, "%f %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g\n", t,nu_top, nu_bot, nu_left, nu_right, nu_pheno_top, nu_pheno_bot, nu_buon_top, nu_buon_bot);
  fflush(fp);
  
  //Nusselts local
  scalar Nuloc[];
  foreach_boundary(bottom)
    Nuloc[] = dA*(T[0,-1] - T[])/Delta;
    
  char names[256];
  sprintf(names, "Nulocbot%2.2e-%d.dat", Ranf, pid());
  FILE * fp1 = fopen(names, "w");
  foreach_boundary(bottom)
    fprintf (fp1, "%f %.9g %.9g %.9g %.9g %.9g\n", t, x, y, Nuloc[], T[0,-1], T[]);
  fflush(fp1);
  fclose(fp1);
  char command[256];
  sprintf(command, "LC_ALL=C cat Nulocbot%2.2e*.dat >> Nulocbot.dat && rm Nulocbot%2.2e*.dat", Ranf, Ranf);
  system(command);
  
  scalar Nuloct[];
  foreach_boundary(top)
    Nuloct[] = dA*(T[] - T[0,1])/Delta;
  char namess[256];
  sprintf(namess, "Nuloctop%2.2e-%d.dat", Ranf, pid());
  FILE * fp2 = fopen(namess, "w");
  foreach_boundary(top)
    fprintf (fp2, "%f %.9g %.9g %.9g %.9g %.9g\n", t, x, y, Nuloct[], T[0,1], T[]);
  fflush(fp2);
  fclose(fp2);
  char command2[256];
  sprintf(command2, "LC_ALL=C cat Nuloctop%2.2e*.dat >> Nuloctop.dat && rm Nuloctop%2.2e*.dat", Ranf, Ranf);
  system(command2);
}

/** 
    ## Output movies
*/

event movie (t += 5.; t <= EndTime){
#if 1
  clear();
  box();
  view (quat = {0.000, 0.000, 0.000, 1.000},
      fov = 30, near = 0.01, far = 1000,
      tx = -0.023, ty = 0.034, tz = -2.387);
  //view(width = 1920, height=1080);
  char time1[80];
  sprintf(time1, "t=%.2gTauD",t*sq(L0)/(Kappaf*sqrt(Ranf)));
  draw_string(time1,pos=0);
  #if nanofluid
  char time2[80];
  sprintf(time2, "t=%.2gTaubrow",t*sq(L0)/(Kappap*sqrt(Ranf)*Lenf));
  draw_string(time2,pos=1);
  char time3[80];
  sprintf(time3, "t=%.2gTauthermo",t*(sq(L0)*St)/(Kappap*sqrt(Ranf)*Lenf));
  draw_string(time3,pos=2);
  #endif
  squares ("T", min = Ttop, max = Tbot);
  char name[80];
  sprintf(name, "temperature_%2.2e.mp4", Ranf);
  save(name);
  
#if nanofluid
  clear();
  box();
  //view(width = 1920, height=1080);
  draw_string(time1,pos=0);
  draw_string(time2,pos=1);
  draw_string(time3,pos=2);
  stats s = statsf(Np);
  squares ("Np", min = (s.max - s.stddev/10.), max = s.max);
  char name2[80];
  sprintf(name2, "nanoparticles_%2.2e.mp4", Ranf);
  save(name2);
#endif

  clear();
  box();
  //view(width = 1920, height=1080);
  draw_string(time1,pos=0);
  #if nanofluid
  draw_string(time2,pos=1);
  draw_string(time3,pos=2);
  #endif
  cells();
  scalar l[];
  foreach()
    l[] = level;
  squares("l", min = MINLEVEL, max = MAXLEVEL);
  char name3[80];
  sprintf(name3, "cells%2.2e.mp4", Ranf);
  save(name3);
  
  clear();
  box();
  //view(width = 1920, height=1080);
  draw_string(time1,pos=0);
  #if nanofluid
  draw_string(time2,pos=1);
  draw_string(time3,pos=2);
  #endif
  squares ("T", min = Ttop, max = Tbot, linear = true, map = cool_warm);
  for (float val = -0.4; val < 0.4; val=val+0.1) {
    isoline ("T", val, min = -0.5, max = 0.5,  linear = true);
  }
  save ("isoT.mp4");
  
  clear();
  box();
  //view(width = 1920, height=1080);
  draw_string(time1,pos=0);
  #if nanofluid
  draw_string(time2,pos=1);
  draw_string(time3,pos=2);
  squares ("Np", min = 0., max = s.max*1.2, linear = true, map = cool_warm);
  for (float val = 0.; val < 1.; val=val+0.1) {
    isoline ("Np", val, min = 0, max = s.max*1.2,  linear = true);
  }
  save ("phi.mp4");
  #endif

  clear();
  box();
  //view(width = 1920, height=1080);
  draw_string(time1,pos=0);
  #if nanofluid
  draw_string(time2,pos=1);
  draw_string(time3,pos=2);
  #endif
  squares("l", min = MINLEVEL, max = MAXLEVEL);
  char name4[80];
  sprintf(name4, "level%2.2e.mp4", Ranf);
  save(name4);
  #endif
  
  scalar omega[];
  vorticity (u, omega);
  
  stats ss = statsf(omega);
  char min[80];
  sprintf(min, "min=%g", ss.min);
  char max[80];
  sprintf(max, "max=%g", ss.max);
  
  clear();
  box();
  //view(width = 1920, height=1080);
  draw_string(min,pos=0);
  draw_string(max,pos=1);
  isoline ("omega", n = 21);
  squares("omega", min = ss.min, max = ss.max);
  
  char video[80];
  sprintf(video, "omega%2.2e.mp4", Ranf);
  save(video);
}

#if 1
event dumping(t=0; t+=50; t<=EndTime){
    dump();
}

event backup(t+=10.; t<=EndTime){
  scalar omega[];
  vorticity (u, omega);
  char name[80];
  //sprintf(name, "./paraview/field_time%g_cpu%i",t, pid());
  sprintf(name, "./paraview/field_time%g_cpu",t);
  //output_vtu((scalar *) {T,Np}, (vector *) {u}, name);
#if nanofluid
  //output_vtu2((scalar *) {T,Np}, (vector *) {u}, name, t);
  output_vtu((scalar *) {T,Np,rbrow,rthermo,rc,omega}, (vector *) {u,D,Dc}, name);
#else
  output_vtu((scalar *) {T,omega}, (vector *) {u}, name);
#endif
}

void printingend(){
  clear();
  view (quat = {0.000, 0.000, 0.000, 1.000},
      fov = 30, near = 0.01, far = 1000,
      tx = -0.023, ty = 0.034, tz = -2.387,
      width = 1920, height=1080);
  box();
  squares("T", linear = true, map = cool_warm, min = -0.45, max = 0.45);
  isoline("T", n = 21, min = -0.45, max = 0.45);
  save ("isoT.png");
  
  clear();
  view (quat = {0.000, 0.000, 0.000, 1.000},
      fov = 30, near = 0.01, far = 1000,
      tx = -0.023, ty = 0.034, tz = -2.387,
      width = 1920, height=1080);
  scalar psi[];
  streamfunction (u, psi);
  boundary ({psi});
  squares ("psi", linear = true, map = cool_warm);
  isoline ("psi", n = 21);
  box();
  save ("isoPsi.png");
}

void printdata(){
  char namess[256];
  sprintf(namess, "FIELD_%2.2e_%g.dat", Ranf, t);
  FILE * fg = fopen (namess, "w");
  scalar psi[];
  streamfunction (u, psi);
  scalar omega[];
  vorticity (u, omega);
  #if nanofluid
  output_field({T, u.x, u.y, omega, psi, Np},fg, linear=true);   
  #else
  output_field({T, u.x, u.y, omega, psi},fg, linear=true);
  #endif
}

void printbindata(){
  char nameT[128];
  sprintf(nameT, "Tfield_%2.2e.bin", Ranf);
  FILE * fT = fopen (nameT, "w");
  output_matrix_mpi(T,fT, linear=true);
  fclose(fT);
  char nameux[128];
  sprintf(nameux, "Uxfield_%2.2e.bin", Ranf);
  FILE * fux = fopen (nameux, "w");
  output_matrix_mpi(u.x,fux, linear=true);
  fclose(fux);
  char nameuy[128];
  sprintf(nameuy, "Uyfield_%2.2e.bin", Ranf);
  FILE * fuy = fopen (nameuy, "w");
  output_matrix_mpi(u.y,fuy, linear=true);
  scalar psi[];
  streamfunction (u, psi);
  scalar omega[];
  vorticity (u, omega);
  char nameomg[128];
  sprintf(nameomg, "Omgfield_%2.2e.bin", Ranf);
  FILE * fomg = fopen (nameomg, "w");
  output_matrix_mpi(omega,fomg, linear=true);
  fclose(fomg);
  char namepsi[128];
  sprintf(namepsi, "Psifield_%2.2e.bin", Ranf);
  FILE * fpsi = fopen (namepsi, "w");
  output_matrix_mpi(psi,fpsi, linear=true);
  fclose(fpsi);
  #if nanofluid
  char nameNp[128];
  sprintf(nameNp, "Npfield_%2.2e.bin", Ranf);
  FILE * fNp = fopen (nameNp, "w");
  output_matrix_mpi(Np,fNp, linear=true);
  fclose(fNp);
  #endif
}


event pictures(t=EndTime){
  printingend();
  printbindata();
  printdata();
  profile({u.y}, y, "uy.dat");
}

scalar un[];
event init_un (i = 0) {
  foreach()
    un[] = u.x[];
}

event logfile (t += 1.0; t <= EndTime) {
  double deltau = change (u.x, un);
  if (deltau <= 1e-5 && t > 200.){
      printingend();
      printbindata();
      printdata();
      dump("converged");
      profile({u.y}, y, "uy.dat");
      return 1;
  }
}

/** ## Results
![Temperature field.](nanofluid/Temperature*.mp4)
![nanoparticle field.](nanofluid/nanoparticle*.mp4)

~~~gnuplot
set xlabel 'times t'
set ylabel 'Nusselts'
plot '< cat Nusselts*.dat' u 1:2 w l title 'Nu-top',\
      '' u 1:3 w l title 'Nu-bot'
~~~
*/
/** ## References
Buongiorno, J.. (03.2006). Convective Transport in Nanofluids. Journal of Heat Transfer. (128)3. p.240 - 250.

Corcione, M., Cianfrini, M., & Quintino, A. (2013). Two-phase mixture modeling of natural convection of nanofluids with temperature-dependent properties. International Journal of Thermal Sciences, 71, 182-195.
*/