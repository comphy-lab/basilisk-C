#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"
#include "output_hdf.h"
#include "signature.h"

scalar f[], cf[], phii[], sign[], M[];
scalar * interfaces = {f};

const int max_change = 10; // Maximum number of holes 
bool large = true;  //Make larger holes

#define rho1 900.
#define rho2 1000.
#define mu1  0.0195756
#define mu2  0.0195756
#define rho(f) (clamp(f,0,1)*(rho1 - rho2) + rho2)
#define mu(f) (1./(clamp(f,0,1)*(1./(mu1) - 1./(mu2)) + 1./(mu2)))
#define sigmaf 0.0153281

#define H 0.1
#define G 9.81
#define Ug sqrt((rho2 - rho1)/rho2*H*G/2.)
#define tc (H/(Ug))

#define Ar (rho2 - rho1)*rho2*G*cube(H)/sq(mu2)
#define Bo (rho2 - rho1)*G*sq(H)/sigmaf

face vector alphav[], muv[], av[];
scalar rhov[];

#define Nl 128

/**
  We set the grid dimension and position.
*/

double min_x = -H/2.;
double min_y = -H/2.;
double step = H/Nl;

timer tt;
int count_it = 0;

int main (int argc, char * argv[]) {
  fprintf(ferr, "# tc is %g \n", tc);
  fprintf(ferr, "# Ar^1/2 is %g \n", sqrt(Ar));
  fprintf(ferr, "# Bo is %g \n", Bo);
  size (H);
  origin (-H/2., -H/2., -H/2.);
  f.tracers = {cf};
  init_grid(Nl);
  a = av;
  mu = muv;
  alpha = alphav;
  rho = rhov;
  f.sigma = sigmaf;
  DT = 5e-5;
  tt = timer_start();
  run();
}

/**
  We compute the signature and perforate thin regions
*/

event calc_and_print (t = 0.015; t += 0.015){
  
  foreach(){
    phii[] = 2*f[] - 1;
    sign[] = 7;
  }
  
  int l_sign = 9;
  
  for (int ilev = depth() - 1; ilev >= l_sign; ilev--)  
    foreach_level(ilev){
      if(is_refined(cell))
      restriction_average(point, phii);
    }
  
  compute_signature_neigh_level (f, phii, sign, l_sign);
  
  /** 
   The signature `sign` is available only at the level `l_sign`. 
   We need to prolong it onto the finest grid. */
  
  printf("\n level used for moments %d and depth is %d \n", l_sign, depth()); 
  
  for (int ilev = l_sign; ilev < depth(); ilev++)  
    foreach_level(ilev){
      sign.prolongation = phii.prolongation = refine_injection;
      if(is_refined(cell)){
        sign.prolongation (point, sign);
        phii.prolongation (point, phii);
      }
    }

 change_topology (f, sign, M, l_sign, max_change, large);
  
  /**
     Finally we output the results in HDF5 format.  */
  
  scalar * list = {sign, f};
  vector * vlist = {u};
  
//   output_ppm (sign, linear = true, file = "f.mp4", n = 200);
  char buf[100];
  char time[100];
  char itime[100];
  count_it++;
  
  snprintf(buf, sizeof(buf), "out_%05d.xmf", count_it);
  snprintf(itime, sizeof(itime), "%05d", count_it);
  snprintf(time, sizeof(time), "%06g", t/tc);
  
  FILE * fp = fopen(buf, "w");
  output_xmf_h5_foreach(list, vlist, 64, fp, itime, time); 
  fclose(fp);
  
  char dumpname[80];
  sprintf (dumpname, "snapshot-%g", t/tc);
  dump (dumpname);
  
}


event init (i = 0) {
  #if TREE
  scalar f1[];
  foreach()
    f1[] = (x <= 0 && y <= 0 && z <= 0);
  boundary ({f1});
  astats s;
  do {
    s = adapt_wavelet ({f1}, (double[]){0.0}, maxlevel=9, minlevel=4, list = NULL);
    foreach()
      f1[] = (x <= 0 && y <= 0 && z <= 0);
    boundary ({f1});
  } while (s.nf);
  foreach()
    f[] = (x <= 0 && y <= 0 && z <= 0);
  boundary ({f});
  #else
  foreach()
    f[] = (x <= 0 && y <= 0 && z <= 0);
  boundary ({f});
  #endif
}

event acceleration (i++) {
  foreach_face(y)
    av.y[] -= G;
  boundary ((scalar *){av});
}

event properties (i++) {
  #if TREE
  f.prolongation = refine_bilinear;
  boundary ({f});
  #endif
  
  foreach_face() {
    double ff = (f[] + f[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff);
    muv.x[] = fm.x[]*mu(ff);
  }
  foreach()
    rhov[] = cm[]*rho(f[]);
  
  #if TREE
  f.prolongation = fraction_refine;
  boundary ({f});
  #endif
}

event logfile (i++; t <= 3) {
  double ke1 = 0., ke2 = 0.; //kinetic energy 
  double vd = 0., vol1 = 0.; //viscous dissipation and volume at the top
  double ep1 = 0., ep2 = 0.; //potentail energy
  double er1 = 0., er2 = 0.; //enstrophy
  double er_nb_1 = 0., er_nb_2 = 0.;  //enstrophy inside the domain 
  double area = 0.;
  int nc = 0;
  static long tnc = 0;
  
  foreach(reduction(+:ke1) reduction(+:ke2) reduction(+:vd)
	  reduction(+:vol1) reduction(+:ep1) reduction(+:ep2)
	  reduction(+:er1) reduction(+:er2) reduction(+:area)
	  reduction(+:nc) ) {
    if (y > H/2. - H/8.)
      vol1 += f[]*dv();
    ep1 += rho1*f[]*G*(y + H/2.)*dv();
    ep2 += rho2*(1. - f[])*G*(y + H/2.)*dv();
    // interfacial area
    if (f[] > 1e-4 && f[] < 1. - 1e-4) {
      coord m = mycs (point, f);
      double alpha = plane_alpha (f[], m);
      coord p;
      area += sq(Delta)*plane_area_center (m, alpha, &p);
    }
    double w2 = 0.;
    foreach_dimension() {
      // kinetic energy
      ke1 += dv()*f[]*rho1*sq(u.x[]);
      ke2 += dv()*(1. - f[])*rho2*sq(u.x[]);
      // viscous dissipation
      vd += dv()*(sq(u.x[1] - u.x[-1]) +
		  sq(u.x[0,1] - u.x[0,-1]) +
		  sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
      // enstrophy
      w2 += sq(u.x[0,1] - u.x[0,-1] - u.y[1,0] + u.y[-1,0]);
    }
    w2 /= sq(2.*Delta);
    er1 += dv()*f[]*w2;
    er2 += dv()*(1. - f[])*w2;
    nc++;
  }
  
  foreach(reduction(+:er_nb_1) reduction(+:er_nb_2)) {
    double w2 = 0.;
    foreach_dimension() {
      // enstrophy
      w2 += sq(u.x[0,1] - u.x[0,-1] - u.y[1,0] + u.y[-1,0]);
    }
    w2 /= sq(2.*Delta);
    if ((x < H/2. - H/10.) && (z < H/2. - H/10.)){
      er_nb_1 += dv()*f[]*w2;
      er_nb_2 += dv()*(1. - f[])*w2;
    }
  }
  
  ke1 /= 2.;
  ke2 /= 2.;
  er1 /= 2.;
  er2 /= 2.;
  er_nb_1 /= 2.;
  er_nb_2 /= 2.;
  //  vd *= MU/vol;

  if (i == 0)
    fprintf (ferr,
	     "# t ke1 ke2 ep1 ep2 er1 er2 R2 area mgp.i mgu.i nc time speed er_nb_1 er_nb_2\n");
  double elapsed = timer_elapsed (tt);
  tnc += nc;
  fprintf (ferr, "%g %g %g %g %g %g %g %g %g %d %d %d %g %g %g %g\n",
	   t/tc,
	   ke1/(1./16.*rho1*sq(Ug)*cube(H)),
	   ke2/(1./16.*rho2*sq(Ug)*cube(H)),
	   ep1/(rho1*G*15.*sq(H)*sq(H)/128.),
	   ep2/(rho2*G*49.*sq(H)*sq(H)/128.),
	   er1,
	   er2,
	   8.*vol1/cube(H),
	   area,
	   mgp.i, mgu.i, nc, elapsed,
     tnc/elapsed,
     er_nb_1,   
     er_nb_2
  );
}


#if TREE
event adapt (i++) {
  adapt_wavelet ({f, u}, (double[]){0.0005,0.005,0.005,0.005},
                   maxlevel = 9, minlevel = 4);
}
#endif
