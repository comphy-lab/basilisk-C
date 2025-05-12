#include "grid/multigrid3D.h"
#include "momentum.h"
#include "tension.h"
#include "reduced.h"
#include "tag.h"
#include "view.h"

#define STEEPNESS 0.55
#define BO 200.
#define RE 40000.0
#define LEVEL 9

#define RATIO (1./850.)
#define MURATIO (17.4e-6/8.9e-4)

double ak = STEEPNESS ;
double k_ = 2.0*M_PI;
double h_ = 0.5;
double g_ = 1.;
double lam_ = 1.0;

double wave (double x, double y) {
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  //double eta3 = 0.0;
  return eta1 + ak*eta2 + sq(ak)*eta3 - y;
  //return eta1 - y;
  //return eta1 + eta2 + eta3 - y;
  //return eta1 + eta2 - y;
  //return eta1 + ak*eta2 - y;
}

double eta (double x) {
  double eta1 = ak/k_*cos(k_*x);
  double a_ = ak/k_;
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  //return eta1 + ak*eta2 + sq(ak)*eta3 ;
  return eta1;
  //return eta1 + eta2 + eta3;
}

/*double gaus (double y, double yc, double T) {
  double deltaw = (sqrt(2./RE)/k_);
  double deltaa = (sqrt(2./(RE*RATIO/MURATIO))/k_);
  double r = y - yc ;
  
  return 2./(sqrt(2.*M_PI*deltaw*deltaw)+sqrt(2.*M_PI*deltaa*deltaa))*
    (T*exp(-r*r/(2.*deltaw*deltaw))+(1.-T)*exp(-r*r/(2.*deltaa*deltaa))); 
}*/

int main()
{
  X0 = Y0 = -L0/2.;
  periodic (right);
  periodic (front);
  rho1 = 1.;
  rho2 = RATIO;
  mu1 = 1.0/RE; //using wavelength as length scale
  mu2 = 1.0/RE*MURATIO;
  //mu1 = 0.0;
  //mu2 = 0.0;
  f.sigma = 1./(BO*sq(k_));
  //f.sigma = 1.0/BO;
  //f.sigma = 0.0;
  G.y = -g_;
  N = 1 << LEVEL;
  run();
}

/*event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= g_;
}*/


event init (i = 0)
{
  //restore("dumpert1.6");
  fraction (f, wave(x,y));

  scalar Phi[];
  foreach() {
    double alpa = 1./tanh(k_*h_);
    double a_ = ak/k_;
    double sgma = sqrt(g_*k_*tanh(k_*h_)*
			(1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					   (sq(alpa) - 1.) + sq(alpa))));
    //double sgma = sqrt(g_*k_*tanh(k_*h_));
    double A_ = a_*g_/sgma;
   // double A_ = a_*sgma/k_;
    double phi1 = A_*cosh(k_*(y + h_))/cosh(k_*h_)*sin(k_*x);
    //double phi1 = A_*cosh(k_*(y+h_))/sinh(k_*h_)*sin(k_*x);
    double phi2 = 3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
      //cosh(2.*k_*(y + h_))*sin(2.*k_*x);
      cosh(2.0*k_*(y + h_))*sin(2.0*k_*x)/cosh(2.0*k_*h_);
    double phi3 = 1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
      (9.*sq(alpa) - 13.)*
      cosh(3.*k_*(y + h_))/cosh(3.*k_*h_)*a_*a_*k_*k_*A_*sin(3.*k_*x);
    //double phi3=0.0;
    Phi[] = phi1 + ak*phi2 + ak*ak*phi3;
    //Phi[] = phi1 + phi2 + phi3;
    //Phi[] = phi1;
    //Phi[] = phi1 + phi2;
    //u.x[] = A_*k_*cosh(k_*(y+h_))/cosh(k_*h_)*cos(k_*x)*f[];
    //u.y[] = A_*k_*sinh(k_*(y+h_))/cosh(k_*h_)*sin(k_*x)*f[];
  }
  boundary ({Phi});


  foreach() {
  //Define density field
  //rho[] = rho1*f[]+rho2*(1.0-f[]);
  //Define momentum field from potential flow
  q.x[] = (Phi[1,0] - Phi[-1,0])/(2.*Delta)*f[]*rho1;
  q.y[] = (Phi[0,1] - Phi[0,-1])/(2.*Delta)*f[]*rho1;
  q.z[] = 0.0;
  }
  //boundary ({u.x,u.y});
  //u.n[left] = dirichlet(-c_);
  //u.n[right] = dirichlet(-c_);
  //u.n[left] = neumann(0.0);
  //u.n[right] = neumann(0.0);
  //u.t[left] = neumann(0.0);
  //u.t[right]= neumann(0.0);
  //u.t[left] = dirichlet(0.0);
  //u.t[right] = dirichlet(0.0);

}

double getDissipationRate(vector u)
{
  //tensor SDeform[];
  scalar SDeformxx[], SDeformyy[];
  scalar SDeformxy[], SDeformyx[];
  double DissipationRate=0.0;
  scalar dudy[];
  scalar dvdx[];
#if dimension > 2
  scalar SDeformzz[];
  scalar SDeformxz[], SDeformzx[];
  scalar SDeformyz[], SDeformzy[];
  scalar dudz[];
  scalar dwdx[];
  scalar dvdz[];
  scalar dwdy[];
#endif
  //equivalent but simpler
  foreach() {
    dudy[] = (u.x[0,1]-u.x[0,-1])/(2.0*Delta);
    dvdx[] = (u.y[1]  -u.y[-1]  )/(2.0*Delta);
#if dimension > 2
    dudz[] = (u.x[0,0,1]-u.x[0,0,-1])/(2.0*Delta);
    dwdx[] = (u.z[1]    -u.z[-1]    )/(2.0*Delta);
    dvdz[] = (u.y[0,0,1]-u.y[0,0,-1])/(2.0*Delta);
    dwdy[] = (u.z[0,1]  -u.z[0,-1]  )/(2.0*Delta);
#endif
    SDeformxx[] = 0.0; //Incompressible fluid
    SDeformyy[] = 0.0;
    //SDeform.x.y[] = 0.5 *( (u.x[0,1]-u.x[0,-1])/(2.0*Delta) +
    //                       (u.y[1]  -u.y[-1]  )/(2.0*Delta) );
    SDeformxy[] = 0.5*(dudy[] + dvdx[]);
    SDeformyx[] = SDeformxy[];
#if dimension > 2
    SDeformzz[] = 0.0; //Incompressible fluid
    //SDeformxz[] = 0.5 *( (u.x[0,0,1]-u.x[0,0,-1])/(2.0*Delta) +
    //                       (u.z[1]    -u.z[-1]    )/(2.0*Delta) );
    //SDeformyz[] = 0.5 *( (u.y[0,0,1]-u.y[0,0,-1])/(2.0*Delta) +
    //                       (u.z[0,1]  -u.z[0,-1]  )/(2.0*Delta) );
    //Equivalent but simpler:
    SDeformxz[] = 0.5*(dudz[] + dwdx[]);
    SDeformzx[] = SDeformxz[];
    SDeformyz[] = 0.5*(dvdz[] + dwdy[]);
    SDeformzy[] = SDeformyz[];
#endif
  }
  //Now get tensor product Sij.Sij
  double TensorProduct = 0.0;
  foreach(reduction(+:DissipationRate)){
    //TensorProduct = 2.0* (sq(SDeformxy[]) + sq(SDeformxz[]) + sq(SDeformyz[]));
    //DissipationRate += 2.0*mu1/rho[]*TensorProduct*cube(Delta) * f[];
#if dimension < 3
    DissipationRate += 2.0*mu1/rho[]*f[] *
      (sq(SDeformxy[])+sq(SDeformyx[]))*sq(Delta);
#else
    DissipationRate += 2.0*mu1/rho[]*f[] *
      (sq(SDeformxy[])+sq(SDeformyx[]) +
       sq(SDeformxz[])+sq(SDeformzx[]) +
       sq(SDeformyz[])+sq(SDeformzy[]))*cube(Delta);
#endif
  }
  return DissipationRate;
}

event countDropsBubble(i++)
{
  //Count the droplets and bubbles...
  scalar m1[]; //droplets
  scalar m2[]; //bubbles
  foreach(){
    m1[] = f[] > 0.5; //i.e. set m true if f[] is locally nonzero (droplets)
    m2[] = f[] < 0.5; //m true if f[] locally (close to) zero (bubbles)
  }
  int n1 = tag(m1);
  int n2 = tag(m2);
  //and find their size. This example is due to the jet atomization problem.
  //Determine volumes v and positions b of each droplet/bubble
  double v1[n1]; //droplet
  coord b1[n1];  //droplet
  double v2[n2]; //bubble
  coord b2[n2];  //bubble
  
  //Initialize... on droplets
  for (int j=0; j<n1; j++)
    {
      v1[j] = b1[j].x = b1[j].y = b1[j].z = 0.0;
    }
  //Initialize... on bubbles
  for (int j=0; j<n2; j++)
    {
      v2[j] = b2[j].x = b2[j].y = b2[j].z = 0.0;
    }
  //... and calculate
  foreach_leaf() //droplets
    {
      if (m1[] > 0) { //i.e. if the stencil-centred value of f[] is true,
      int j = m1[] - 1; //i.e. set j to zero? Or the local tag number?
      v1[j] += dv()*f[]; //increment the volume of the droplet
      coord p = {x,y,z};
      foreach_dimension()
	b1[j].x += dv()*f[]*p.x;
      }
    }
  foreach_leaf() //bubbles
    {
      if (m2[] > 0) {
	int j = m2[] - 1;
	v2[j] += dv()*f[];
	coord p = {x,y,z};
	foreach_dimension()
	  b2[j].x += dv()*f[]*p.x;
      }
    }
  //Reduce for purposes of MPI
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v1, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b1, 3*n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, v2, n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b2, 3*n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  //and output the volume and position of each droplet to file,
  static FILE * fdrop = fopen("droplets.dat","w");
  static FILE * fbubb = fopen("bubbles.dat","w");
  for (int j=0; j<n1; j++)
    {
      fprintf (fdrop, "%d %g %d %g %g %g\n", i, t,
	       j, v1[j], b1[j].x/v1[j], b1[j].y/v1[j]);
      fprintf (ferr, "Droplets: %d %g %d %g %g %g\n", i, t,
	       j, v1[j], b1[j].x/v1[j], b1[j].y/v1[j]);
    }
  for (int j=0; j<n2; j++)
    {
      fprintf (fbubb, "%d %g %d %g %g %g\n", i, t,
	       j, v2[j], b2[j].x/v2[j], b2[j].y/v2[j]);
      fprintf (ferr, "Bubbles: %d %g %d %g %g %g\n", i, t,
	       j, v2[j], b2[j].x/v2[j], b2[j].y/v2[j]);
    }
}


/* Calculate total kinetic energy in the water*/
double calculate_ke()
{ 
  double kesum = 0.0;
  foreach(reduction(+:kesum))
  {
    kesum += (q.x[]*q.x[]/rho[] + q.y[]*q.y[]/rho[] + q.z[]*q.z[]/rho[])*f[] * cube(Delta);
  }
  return 0.5*kesum;
}
/* Calculate the gravitational potential energy*/
double calculate_gpe()
{
  double gpesum = 0.0;
  foreach(reduction(+:gpesum))
  {
    gpesum += rho1 * g_ * y * f[] * cube(Delta);
  }
  return gpesum + 0.125;
}

/* Calculate the surface tension potential energy*/
/*double calculate_stpe(scalar LL)
{
  double gam = (rho1-rho2)*g_/(BO*sq(k_));
  double xo = 0.0;
  double yo = 0.0;
  scalar * list = NULL;
  foreach()
    {
      if (f[] < 1. && f[] > 0.)
        {
          xo = x;
	  yo = y;
        }
    }   
  return 0.0;
}*/

int get_umag(scalar umag){
  foreach()
    {
      umag[] = sqrt(sq(q.x[])/rho[]+sq(q.y[])/rho[]+sq(q.z[])/rho[]);
    } 
  return 1;
}

int get_u(vector u){
  foreach()
    {
      u.x[] = q.x[]/rho[];
      u.y[] = q.y[]/rho[];
      u.z[] = q.z[]/rho[];
    }
  return 1;
}
  

/*Total energy, not counting losses*/

/*Graphs*/
event graphs (i++){
  double ke = calculate_ke();
  double gpe = calculate_gpe();
  vector u[];
  get_u(u);
  double dissipationRate = getDissipationRate(u);
  static FILE * fp = fopen ("budget.dat", "w");
  fprintf (fp, "%g %g %g %g %g\n", t/(k_/sqrt(g_*k_)), ke, gpe, dissipationRate, ke+gpe);
  fprintf (ferr, "%g %g %g %g %g\n", t/(k_/sqrt(g_*k_)), ke, gpe, dissipationRate, ke+gpe);
}



/*event imagesi(i++){
  char targetn[100];
  sprintf(targetn,"f%0.8d.ppm",i);
  double pp = 0.2;
  double th = 0.5;
  view (fov=30,tx=0.3,ty=0.2, theta=th,phi=pp);
  //box();
  draw_vof("f");
  save(targetn);
  clear();
}*/

event imagest(t+=0.005){
  char targetn[100];
  sprintf(targetn,"f%3.3g.ppm",t);
  double pp=0.2;
  double th = 0.5;
  view (fov=30,tx=0.3,ty=0.2, theta=th,phi=pp);
  box();
  draw_vof("f");
  save(targetn);
  clear();
}

/*event movie (t += 0.02) {
  static FILE * fp = popen ("ppm2mpeg > f.mpg", "w");
  output_ppm (f, fp, min =0, max =1, n=1024);
  vector u[];
  get_u(u);
  static FILE * fp2 = popen ("ppm2mpeg > umag.mpg", "w");
  scalar umag[];
  get_umag(umag);
  output_ppm (umag, fp2, n=1024);
  static FILE * fp3 = popen ("ppm2mpeg > ux.mpg", "w");
  output_ppm (u.x, fp3, n=1024);
  static FILE * fp4 = popen ("ppm2mpeg > uy.mpg", "w");
  output_ppm (u.y, fp4, n=1024);
  static FILE * fp7 = popen ("ppm2mpeg > uz.mpg", "w");
  output_ppm (u.z, fp7, n=1024);
  static FILE * fp5 = popen ("ppm2mpeg > rho.mpg", "w");
  output_ppm (rho, fp5, n=1024);
  static FILE * fp6 = popen ("ppm2mpeg > p.mpg", "w");
  output_ppm (p, fp6, n=1024);
}*/

/*event movie (t += 0.02) {
  double p = 1.0/5.0;
  double t= 0.5;
  view (fov=40,tx=-0.0,ty=-0.0, theta=t, phi=p);
  box();
  static FILE * fp = popen ("ppm2mpeg > test3d.mpg", "w");
  draw_vof("f");
  save (fp = fp);
  clear();
  
}*/

/*event images (t += 0.01) {
  static FILE * fp1 = fopen ("im.ppm", "w");
  scalar omega[];
  vorticity (u, omega);
  output_ppm (omega, fp1, n=1024);
}*/
/*
event images (t += 0.01) {
  static FILE * fp1 = fopen ("imf.ppm", "w");
  scalar f[];
  output_ppm (f, fp1, min = 0, max = 1, n=1024);
}*/

event dumpstep(t += 0.1) {
  char name[100];
  sprintf(name,"dumpert%g", t);
  dump(name);
}

event end (t = 4.0*k_/sqrt(g_*k_)) { /* k_/sqrt(g_*k_) is the wave period T - we want to run up to 4 non-dimensional time units */
  printf ("i = %d t = %g\n", i, t);
}

